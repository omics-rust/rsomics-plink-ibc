//! PLINK1 `--ibc`: the three GCTA method-of-moments inbreeding estimators
//! per sample, over non-missing autosomal markers.
//!
//! With `x` = A1-allele count (0/1/2), `p` = founder A1 frequency, `d = 2p(1-p)`,
//! the per-marker contributions (Yang et al. 2011 AJHG) are
//!   Fhat1 = (x - 2p)^2 / d - 1                  variance-standardized
//!   Fhat2 = 1 - x(2-x) / d                       excess homozygosity (= --het F)
//!   Fhat3 = (x^2 - (1+2p)x + 2p^2) / d           correlation of uniting gametes
//! each averaged over the sample's `NOMISS` markers.

#![allow(clippy::cast_precision_loss)]

use rayon::prelude::*;
use rsomics_pgen::Pgen;
use std::io::{self, Write};

pub struct IbcRecord {
    pub fid: String,
    pub iid: String,
    pub nomiss: u32,
    pub fhat1: f64,
    pub fhat2: f64,
    pub fhat3: f64,
}

fn founder_mask(pgen: &Pgen) -> Vec<bool> {
    pgen.samples
        .iter()
        .map(|s| s.pid == "0" && s.mid == "0")
        .collect()
}

fn is_autosome(chrom: &str) -> bool {
    matches!(chrom.parse::<u32>(), Ok(1..=22))
}

fn founder_allele_counts(row: &[u8], founders: &[bool], n_samp: usize) -> (u32, u32) {
    let (mut a1, mut a2) = (0u32, 0u32);
    for s in 0..n_samp {
        if !founders[s] {
            continue;
        }
        match (row[s / 4] >> ((s % 4) * 2)) & 0b11 {
            0b00 => a1 += 2,
            0b10 => {
                a1 += 1;
                a2 += 1;
            }
            0b11 => a2 += 2,
            _ => {}
        }
    }
    (a1, a2)
}

/// Per packed byte, the four lanes' A1 dosage (HomA1=2/Het=1/HomA2=0) for the
/// all-founder fast path.
struct AlleleLut {
    a1: [u32; 256],
    a2: [u32; 256],
}

impl AlleleLut {
    fn build() -> Self {
        let mut a1 = [0u32; 256];
        let mut a2 = [0u32; 256];
        for byte in 0usize..256 {
            let mut b = byte as u8;
            for _ in 0..4 {
                match b & 0b11 {
                    0b00 => a1[byte] += 2,
                    0b10 => {
                        a1[byte] += 1;
                        a2[byte] += 1;
                    }
                    0b11 => a2[byte] += 2,
                    _ => {}
                }
                b >>= 2;
            }
        }
        Self { a1, a2 }
    }
}

/// Per-sample sums of the three estimators plus the non-missing count, laid out
/// AoS so one sample's four running totals share a cache line — the inner
/// per-variant loop scatters by sample, so a contiguous 32-byte slot per sample
/// beats four separate strided arrays. Lane 3 holds NOMISS as an `f64`.
#[derive(Clone)]
struct Accum(Vec<[f64; 4]>);

impl Accum {
    fn zeros(n: usize) -> Self {
        Self(vec![[0.0; 4]; n])
    }

    fn merge(mut self, other: Self) -> Self {
        for (a, b) in self.0.iter_mut().zip(&other.0) {
            for k in 0..4 {
                a[k] += b[k];
            }
        }
        self
    }
}

/// Per-2-bit-code contribution slots for one variant: index by the raw PLINK
/// code (00 HomA1, 01 Missing, 10 Het, 11 HomA2). Missing maps to all-zeros so
/// the inner loop is branchless.
fn variant_table(p: f64, d: f64) -> [[f64; 4]; 4] {
    let term = |x: f64| {
        let f1 = (x - 2.0 * p).powi(2) / d - 1.0;
        let f2 = 1.0 - x * (2.0 - x) / d;
        let f3 = (x * x - (1.0 + 2.0 * p) * x + 2.0 * p * p) / d;
        [f1, f2, f3, 1.0]
    };
    let mut t = [[0.0; 4]; 4];
    t[0b00] = term(2.0);
    t[0b10] = term(1.0);
    t[0b11] = term(0.0);
    // t[0b01] (Missing) stays all-zeros.
    t
}

/// Compute the per-sample `--ibc` records.
#[must_use]
pub fn ibc(pgen: &Pgen) -> Vec<IbcRecord> {
    let n_samp = pgen.n_samples();
    let bpv = pgen.bytes_per_variant();
    let founders = founder_mask(pgen);
    let all_founders = founders.iter().all(|&f| f);
    let allele = AlleleLut::build();
    let full_bytes = n_samp / 4;
    let last_lanes = n_samp % 4;

    let usable: Vec<usize> = (0..pgen.n_variants())
        .filter(|&v| is_autosome(&pgen.variants[v].chrom))
        .collect();
    let gt = &pgen.gt_raw;

    let fold = |mut acc: Accum, &v: &usize| -> Accum {
        let row = &gt[v * bpv..v * bpv + bpv];

        let (a1, a2) = if all_founders {
            let (mut a1, mut a2) = (0u32, 0u32);
            for &byte in &row[..full_bytes] {
                a1 += allele.a1[byte as usize];
                a2 += allele.a2[byte as usize];
            }
            for lane in 0..last_lanes {
                match (row[full_bytes] >> (lane * 2)) & 0b11 {
                    0b00 => a1 += 2,
                    0b10 => {
                        a1 += 1;
                        a2 += 1;
                    }
                    0b11 => a2 += 2,
                    _ => {}
                }
            }
            (a1, a2)
        } else {
            founder_allele_counts(row, &founders, n_samp)
        };

        let n_obs = a1 + a2;
        if n_obs == 0 {
            return acc;
        }
        let p = f64::from(a1) / f64::from(n_obs);
        let d = 2.0 * p * (1.0 - p);
        if d == 0.0 {
            return acc;
        }
        let t = variant_table(p, d);
        let slots = &mut acc.0;

        for (byte_idx, &b) in row[..full_bytes].iter().enumerate() {
            let base = byte_idx * 4;
            let codes = [b & 0b11, (b >> 2) & 0b11, (b >> 4) & 0b11, (b >> 6) & 0b11];
            for (lane, &code) in codes.iter().enumerate() {
                let add = &t[code as usize];
                let slot = &mut slots[base + lane];
                for k in 0..4 {
                    slot[k] += add[k];
                }
            }
        }
        if last_lanes != 0 {
            let b = row[full_bytes];
            let base = full_bytes * 4;
            for lane in 0..last_lanes {
                let add = &t[((b >> (lane * 2)) & 0b11) as usize];
                let slot = &mut slots[base + lane];
                for k in 0..4 {
                    slot[k] += add[k];
                }
            }
        }
        acc
    };

    let acc = if rayon::current_num_threads() > 1 {
        usable
            .par_iter()
            .fold(|| Accum::zeros(n_samp), fold)
            .reduce(|| Accum::zeros(n_samp), Accum::merge)
    } else {
        usable.iter().fold(Accum::zeros(n_samp), fold)
    };

    (0..n_samp)
        .map(|s| {
            let slot = acc.0[s];
            let n = slot[3] as u32;
            let inv = if n == 0 { f64::NAN } else { 1.0 / f64::from(n) };
            IbcRecord {
                fid: pgen.samples[s].fid.clone(),
                iid: pgen.samples[s].iid.clone(),
                nomiss: n,
                fhat1: slot[0] * inv,
                fhat2: slot[1] * inv,
                fhat3: slot[2] * inv,
            }
        })
        .collect()
}

/// Write the records in plink's `.ibc` text layout: tab-delimited
/// `FID IID NOMISS Fhat1 Fhat2 Fhat3`, floats as C `%g` (6 significant digits).
pub fn write_ibc<W: Write>(records: &[IbcRecord], out: &mut W) -> io::Result<()> {
    writeln!(out, "FID\tIID\tNOMISS\tFhat1\tFhat2\tFhat3")?;
    for r in records {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            r.fid,
            r.iid,
            r.nomiss,
            fmt_g(r.fhat1),
            fmt_g(r.fhat2),
            fmt_g(r.fhat3),
        )?;
    }
    Ok(())
}

/// Faithful C `printf("%g")` — plink's numeric output format (6 significant
/// digits, trailing zeros stripped, scientific when exponent `< -4` or `>= 6`).
fn fmt_g(x: f64) -> String {
    const P: i32 = 6;
    if x.is_nan() {
        return "nan".to_string();
    }
    if x == 0.0 {
        return "0".to_string();
    }
    let rounded = format!("{:.*e}", (P - 1) as usize, x);
    let exp: i32 = rounded
        .split('e')
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap();

    if !(-4..P).contains(&exp) {
        let (mant, e) = rounded.split_once('e').unwrap();
        let mant = strip_trailing(mant);
        let e: i32 = e.parse().unwrap();
        format!("{mant}e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
    } else {
        let frac = (P - 1 - exp).max(0) as usize;
        strip_trailing(&format!("{x:.frac$}"))
    }
}

fn strip_trailing(s: &str) -> String {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g_formatting_matches_plink_samples() {
        // 6 significant digits, trailing zeros stripped, like C printf("%g").
        assert_eq!(fmt_g(-0.0972656), "-0.0972656");
        assert_eq!(fmt_g(-0.155059), "-0.155059");
        assert_eq!(fmt_g(0.319444), "0.319444");
        assert_eq!(fmt_g(1.0), "1");
        assert_eq!(fmt_g(0.0), "0");
        assert_eq!(fmt_g(-0.5), "-0.5");
        assert_eq!(fmt_g(1.23456e-7), "1.23456e-07");
        assert_eq!(fmt_g(1234567.0), "1.23457e+06");
    }
}
