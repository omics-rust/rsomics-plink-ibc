# rsomics-plink-ibc

The three GCTA method-of-moments inbreeding-coefficient estimators per sample
from a PLINK1 binary fileset — a Rust reimplementation of `plink --ibc`.

For each sample, over its non-missing autosomal markers (`NOMISS` of them), with
`x` the A1-allele count (0/1/2), `p` the founder A1 frequency and `d = 2p(1−p)`:

| column  | estimator                                                              |
|---------|------------------------------------------------------------------------|
| `Fhat1` | variance-standardized: mean of `(x − 2p)² / d − 1`                     |
| `Fhat2` | excess homozygosity (= `--het` F): mean of `1 − x(2 − x) / d`          |
| `Fhat3` | correlation of uniting gametes: mean of `(x² − (1+2p)x + 2p²) / d`     |

`p` is estimated from **founders only** (samples whose parental IDs are both
`0`), matching PLINK's default. Only autosomes (chromosomes 1–22) are counted;
markers with no founder genotypes or a monomorphic founder set are skipped.

## Usage

```sh
# write the .ibc table to stdout
rsomics-plink-ibc path/to/fileset

# write to <out>.ibc instead of stdout (matches plink --out)
rsomics-plink-ibc path/to/fileset --out result

# parallel over markers
rsomics-plink-ibc path/to/fileset -t 8
```

`path/to/fileset` is the prefix shared by `fileset.bed`, `fileset.bim`,
`fileset.fam` (no extension), exactly as PLINK's `--bfile`.

## Compatibility

The `NOMISS`, `Fhat1`, `Fhat2`, and `Fhat3` fields are byte-identical to
PLINK 1.9 (`%.4g` numeric formatting reproduced exactly, including the switch to
scientific notation). `tests/compat.rs` diffs every field against PLINK live
when the `plink` binary is on `PATH`, and against checked-in PLINK 1.9 golden
output otherwise.

The inter-column whitespace padding of PLINK's `.ibc` layout depends on
PLINK-internal column-width bookkeeping that is not recoverable from the file
format alone; the column **values** match field-for-field, so downstream
parsers (which split on whitespace) see identical data.

## Origin

This crate is an independent Rust reimplementation of `plink --ibc` based on:

- The published method: Yang et al. 2011 AJHG (GCTA,
  doi:10.1016/j.ajhg.2010.11.011), from which PLINK ported `--ibc`; and
  Chang et al. 2015 (PLINK 1.9, doi:10.1186/s13742-015-0047-8).
- The public PLINK 1.9 `--ibc` / basic-statistics documentation
  (<https://www.cog-genomics.org/plink/1.9/basic_stats>) and binary-fileset
  format spec (<https://www.cog-genomics.org/plink/1.9/formats>).
- Black-box behaviour testing against the `plink` 1.9 binary.

No source code from the GPL upstream was used as reference during
implementation. Test fixtures are independently generated.

License: MIT OR Apache-2.0.
Upstream credit: [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
(Christopher Chang et al., GPLv3); `--ibc` method from
[GCTA](https://yanglab.westlake.edu.cn/software/gcta/) (Yang et al.).
