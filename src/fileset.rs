//! PLINK1 `.bed/.bim/.fam` reader.
//!
//! Distinct from the foundation `rsomics_pgen::Pgen::load` in one way that
//! matters for `--ibc`: the `.fam` sex column is read leniently. PLINK treats
//! any code other than `1`/`2` (including the `-9` missing marker it writes
//! itself) as ambiguous sex and still runs the analysis; `--ibc` never consults
//! sex, so a strict numeric parse would reject valid filesets for no reason.

use anyhow::{Context, Result, bail};
use rsomics_pgen::{Pgen, Sample, Variant};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub fn load(prefix: &Path) -> Result<Pgen> {
    let variants = parse_bim(&prefix.with_extension("bim"))?;
    let samples = parse_fam(&prefix.with_extension("fam"))?;
    let gt_raw = read_bed_raw(&prefix.with_extension("bed"), variants.len(), samples.len())?;
    Ok(Pgen {
        variants,
        samples,
        gt_raw,
    })
}

fn parse_fam(path: &Path) -> Result<Vec<Sample>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out = Vec::new();
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 6 {
            bail!(
                "malformed {} line {}: {} fields, expected 6",
                path.display(),
                i + 1,
                fields.len()
            );
        }
        let sex = match fields[4] {
            "1" => 1,
            "2" => 2,
            _ => 0,
        };
        out.push(Sample {
            fid: fields[0].to_string(),
            iid: fields[1].to_string(),
            pid: fields[2].to_string(),
            mid: fields[3].to_string(),
            sex,
            phen: fields[5].to_string(),
        });
    }
    Ok(out)
}

fn parse_bim(path: &Path) -> Result<Vec<Variant>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out = Vec::new();
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 6 {
            bail!(
                "malformed {} line {}: {} fields, expected 6",
                path.display(),
                i + 1,
                fields.len()
            );
        }
        let cm: f64 = fields[2].parse().with_context(|| {
            format!("{} line {}: bad cM {:?}", path.display(), i + 1, fields[2])
        })?;
        let pos: u64 = fields[3].parse().with_context(|| {
            format!("{} line {}: bad pos {:?}", path.display(), i + 1, fields[3])
        })?;
        out.push(Variant {
            chrom: fields[0].to_string(),
            id: fields[1].to_string(),
            cm,
            pos,
            a1: fields[4].to_string(),
            a2: fields[5].to_string(),
        });
    }
    Ok(out)
}

fn read_bed_raw(path: &Path, n_variants: usize, n_samples: usize) -> Result<Vec<u8>> {
    let mut f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut header = [0u8; 3];
    f.read_exact(&mut header)?;
    if header[0] != 0x6c || header[1] != 0x1b {
        bail!(
            "{}: not a PLINK1 .bed (magic {:#04x} {:#04x})",
            path.display(),
            header[0],
            header[1]
        );
    }
    match header[2] {
        0x01 => {}
        0x00 => bail!("{}: sample-major .bed unsupported", path.display()),
        other => bail!("{}: unknown .bed mode {other:#04x}", path.display()),
    }
    let bytes_per_variant = n_samples.div_ceil(4);
    let expected = (bytes_per_variant * n_variants) as u64;
    let have = f.metadata()?.len() - 3;
    if have != expected {
        bail!(
            "{}: size mismatch, have {have} payload bytes, expected {expected} \
             ({n_variants} variants x {n_samples} samples)",
            path.display()
        );
    }
    let mut out = vec![0u8; bytes_per_variant * n_variants];
    f.read_exact(&mut out)?;
    Ok(out)
}
