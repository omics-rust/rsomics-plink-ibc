//! Differential compatibility tests against PLINK 1.9 `--ibc`.
//!
//! Each test runs our binary and compares its `.ibc` output to PLINK's
//! field-by-field: FID/IID strings exact, NOMISS integers exact, the three
//! Fhat numeric tokens exact (PLINK prints `%g`, which we reproduce). When the
//! `plink` binary is on PATH we diff live; otherwise we diff against checked-in
//! PLINK 1.9 golden output.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-plink-ibc"))
}

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn plink_available() -> bool {
    Command::new("plink")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn fields(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .collect()
}

fn run_ours_on(prefix: &Path) -> String {
    let out = Command::new(ours())
        .arg(prefix.to_string_lossy().into_owned())
        .output()
        .expect("run rsomics-plink-ibc");
    assert!(
        out.status.success(),
        "rsomics-plink-ibc failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8")
}

fn assert_fields_equal(ours: &str, ref_text: &str) {
    let a = fields(ours);
    let b = fields(ref_text);
    assert_eq!(a.len(), b.len(), "row count differs");
    for (i, (x, y)) in a.iter().zip(&b).enumerate() {
        assert_eq!(x, y, "row {i} differs:\n ours: {x:?}\n ref:  {y:?}");
    }
}

fn run_live_plink(prefix: &Path) -> String {
    let tmp = tempfile::Builder::new()
        .prefix("plink-ibc-compat-")
        .tempdir_in(std::env::var("TMPDIR").unwrap_or_else(|_| "/tmp".into()))
        .expect("tempdir");
    let out_prefix = tmp.path().join("ref");
    let status = Command::new("plink")
        .args([
            "--bfile",
            prefix.to_str().unwrap(),
            "--ibc",
            "--allow-no-sex",
            "--out",
            out_prefix.to_str().unwrap(),
        ])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("run plink");
    assert!(status.success(), "plink --ibc failed");
    std::fs::read_to_string(out_prefix.with_extension("ibc")).expect("read .ibc")
}

#[test]
fn default_matches_golden() {
    let ours = run_ours_on(&golden_dir().join("small"));
    let golden =
        std::fs::read_to_string(golden_dir().join("small.ibc.golden")).expect("read golden");
    assert_fields_equal(&ours, &golden);
}

#[test]
fn autosome_only_matches_golden() {
    // PLINK excludes the X-chromosome variants; only the chr-1 markers count.
    let ours = run_ours_on(&golden_dir().join("withx"));
    let golden =
        std::fs::read_to_string(golden_dir().join("withx.ibc.golden")).expect("read withx golden");
    assert_fields_equal(&ours, &golden);
}

#[test]
fn header_is_plink_shape() {
    let ours = run_ours_on(&golden_dir().join("small"));
    let header: Vec<&str> = ours.lines().next().unwrap().split_whitespace().collect();
    assert_eq!(header, ["FID", "IID", "NOMISS", "Fhat1", "Fhat2", "Fhat3"]);
}

#[test]
fn monomorphic_markers_counted_like_plink() {
    // A monomorphic marker is not skipped: PLINK adds Fhat1 += -1, Fhat2/Fhat3
    // += 0 and includes it in NOMISS. The `mono` fixture is one polymorphic plus
    // one monomorphic marker; the golden is PLINK 1.9 output.
    let ours = run_ours_on(&golden_dir().join("mono"));
    let golden =
        std::fs::read_to_string(golden_dir().join("mono.ibc.golden")).expect("read mono golden");
    assert_fields_equal(&ours, &golden);
}

#[test]
fn missing_sex_minus9_is_accepted() {
    // PLINK treats a `-9` sex code as ambiguous-but-loaded and still runs --ibc;
    // a strict sex parse must not reject the fileset.
    let ours = run_ours_on(&golden_dir().join("nosex"));
    let golden =
        std::fs::read_to_string(golden_dir().join("nosex.ibc.golden")).expect("read nosex golden");
    assert_fields_equal(&ours, &golden);
}

#[test]
fn single_sample_fails_loud_and_writes_no_ibc() {
    // PLINK: "Error: At least 2 people required for pairwise analysis." and no
    // .ibc is produced. We must fail non-zero and not emit a record.
    let tmp = tempfile::Builder::new()
        .prefix("plink-ibc-single-")
        .tempdir_in(std::env::var("TMPDIR").unwrap_or_else(|_| "/tmp".into()))
        .expect("tempdir");
    let out_prefix = tmp.path().join("out");
    let output = Command::new(ours())
        .arg(golden_dir().join("single").to_string_lossy().into_owned())
        .arg("--out")
        .arg(out_prefix.to_string_lossy().into_owned())
        .output()
        .expect("run rsomics-plink-ibc");
    assert!(
        !output.status.success(),
        "expected failure on single-sample fileset"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("At least 2 people required for pairwise analysis."),
        "stderr should carry PLINK's message, got: {stderr}"
    );
    assert!(
        !out_prefix.with_extension("ibc").exists(),
        "no .ibc must be written for a single-sample fileset"
    );
}

#[test]
fn default_matches_live_plink() {
    if !plink_available() {
        eprintln!("plink not on PATH; skipping live differential");
        return;
    }
    let prefix = golden_dir().join("small");
    assert_fields_equal(&run_ours_on(&prefix), &run_live_plink(&prefix));
}

#[test]
fn autosome_matches_live_plink() {
    if !plink_available() {
        eprintln!("plink not on PATH; skipping live differential");
        return;
    }
    let prefix = golden_dir().join("withx");
    assert_fields_equal(&run_ours_on(&prefix), &run_live_plink(&prefix));
}
