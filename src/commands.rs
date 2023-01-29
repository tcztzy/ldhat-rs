use super::LDhatResult as Result;
use clap::Parser;
use std::path::PathBuf;

pub trait Executable {
    fn execute(&self) -> Result<()>;
}

/// Convert FASTA-style file to LDhat format.
#[derive(Parser, Debug)]
pub struct Convert {
    /// Input FASTA-style format file.
    #[arg(value_name = "FILE")]
    seq: PathBuf,
    /// SNP positions in seq file. Assumed contiguous if absent
    #[arg(short, long, value_name = "FILE")]
    loc: Option<PathBuf>,
    /// Only output sites with exactly two alleles
    #[arg(long, default_value_t = false)]
    only2: bool,
    /// Min Minor Allele Frequency (between 0 and 1)
    #[arg(long, default_value_t = 0., value_name = "FLOAT")]
    freqcut: f64,
    /// Max Missing data frequency (between 0 and 1)
    #[arg(long, default_value_t = 1., value_name = "FLOAT")]
    missfreqcut: f64,
    /// Only print sites between these two values
    #[arg(long, value_name = "INT", num_args = 2)]
    sites: Vec<usize>,
    /// Number of sequences to output: default=all
    #[arg(long, value_name = "INT", required = false)]
    nout: Option<usize>,
    /// Prefix of output files
    #[arg(long, value_name = "STRING", default_value = "")]
    prefix: String,
    /// Random seed
    #[arg(long, value_name = "INT")]
    seed: Option<u64>,
    #[clap(flatten)]
    verbose: clap_verbosity_flag::Verbosity,
}

impl Executable for Convert {
    fn execute(&self) -> Result<()> {
        Ok(())
    }
}
