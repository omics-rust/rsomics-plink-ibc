use clap::Parser;
use rsomics_pgen::Pgen;
use rsomics_plink_ibc::{ibc, write_ibc};
use std::fs::File;
use std::io::{self, BufWriter};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "rsomics-plink-ibc",
    about = "PLINK1 method-of-moments inbreeding coefficients Fhat1/Fhat2/Fhat3 per sample (plink --ibc)",
    version
)]
struct Cli {
    /// Path prefix for the .bed/.bim/.fam fileset (without extension).
    bfile: PathBuf,

    /// rayon worker threads (0 = all cores).
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Write the report to <OUT>.ibc instead of stdout (plink --out).
    #[arg(long)]
    out: Option<PathBuf>,
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> anyhow::Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .ok();

    let pgen = Pgen::load(&cli.bfile)?;
    let records = ibc(&pgen);

    match cli.out {
        Some(prefix) => {
            let path = prefix.with_extension("ibc");
            let mut w = BufWriter::new(File::create(path)?);
            write_ibc(&records, &mut w)?;
        }
        None => {
            let stdout = io::stdout();
            let mut w = BufWriter::new(stdout.lock());
            write_ibc(&records, &mut w)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_definition_is_valid() {
        <Cli as clap::CommandFactory>::command().debug_assert();
    }
}
