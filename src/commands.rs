use crate::{
    io::{read_locs, read_sites, Base, Locs, Ploidy},
    LDhatResult as Result, SEQ_MAX,
};
use clap::Parser;
use ndarray_stats::QuantileExt;
use polars::prelude::UInt32Type;
use rand::{rngs::StdRng, SeedableRng};
use std::{fs::File, io::Write, path::PathBuf};

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
    #[arg(long, visible_alias = "2only", default_value_t = false)]
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

fn sc(base: Base) -> char {
    match base {
        Base::N => '?',
        Base::A => '2',
        Base::C => '1',
        Base::G => '0',
        Base::T => '0',
    }
}

impl Executable for Convert {
    fn execute(&self) -> Result<()> {
        // Original use Unix timestamp as seed. From entropy might be a better choice.
        let mut rng: StdRng = SeedableRng::from_entropy();
        let seqs = read_sites(&self.seq)?;
        let (lseq, mut nseq) = seqs.shape();
        let locs = if let Some(loc) = &self.loc {
            read_locs(&loc)?
        } else {
            Locs::new_from_length(lseq)
        };
        if nseq > SEQ_MAX {
            log::warn!(
                "More than max no. sequences: Using first {} for analysis",
                SEQ_MAX
            );
            nseq = SEQ_MAX;
        }
        log::info!(
            "Reading {} sequences of length {} bases .........",
            nseq,
            lseq
        );
        let (lower, upper) = if self.sites.len() == 2 {
            (self.sites[0], self.sites[1])
        } else {
            (0usize, lseq)
        };
        let output_sites_path = PathBuf::from(format!("{}sites.txt", self.prefix));
        let output_locs_path = PathBuf::from(format!("{}locs.txt", self.prefix));
        log::info!(
            "Segregating sites written to file	: {}",
            output_sites_path.to_str().unwrap()
        );
        log::info!(
            "Locations of segregating sites to file	: {}",
            output_locs_path.to_str().unwrap()
        );
        let nout = core::cmp::min(self.nout.unwrap_or(nseq), nseq);
        let mut index = rand::seq::index::sample(&mut rng, nseq, nout).into_vec();
        index.sort();
        let nall = seqs
            .allele_count(Some(self.prefix.as_str()))?
            .to_ndarray::<UInt32Type>()?;
        let mut minor: Vec<u8> = vec![];
        let fl = if self.only2 || self.freqcut > 0. {
            2
        } else {
            1
        };
        let mut output_site = vec![false; lseq];
        let mut fsnp = vec![true; lseq];
        for i in lower..upper {
            let row = nall.slice(ndarray::s![i, 1usize..]);
            let na = row.mapv(|x| (x > 0) as u32).sum();
            let jmin = row
                .mapv(|x| if x > 0 { x as f32 } else { f32::NAN })
                .argmin_skipnan()?;
            let nmin = row[jmin];
            minor.push(jmin as u8 + 1);
            if (fl == 1 && na > 1)
                || (fl == 2
                    && na == 2
                    && nmin != nseq as u32 * seqs.ploidy as u32
                    && nmin as f64 > (nseq * seqs.ploidy as usize) as f64 * self.freqcut)
                    && (nall.index_axis(ndarray::Axis(0), i)[0] as f64
                        <= (nseq * seqs.ploidy as usize) as f64 * self.missfreqcut)
            {
                output_site[i] = true;
                if na > 2 {
                    fsnp[i] = false;
                }
            }
        }
        let psite: usize = output_site.iter().map(|&x| x as usize).sum();
        if psite == 0 {
            return Err(anyhow::anyhow!("No data to output"));
        }
        let mut ofp = File::create(output_sites_path)?;
        writeln!(ofp, "{} {} {}", nout, psite, seqs.ploidy as usize)?;
        let mut loc = File::create(output_locs_path)?;
        write!(
            loc,
            "{} {} {}",
            psite,
            locs.data.last().unwrap(),
            locs.model
        )?;

        for i in index {
            let seq = seqs.data.select_at_idx(i).unwrap();
            ofp.write(format!(">{}\n", seq.name()).as_bytes())?;
            let seq: Vec<Option<u8>> = seq.u8().expect("").into_iter().collect();
            let mut na = 0;
            for (j, &o) in &mut output_site.iter().enumerate() {
                if o {
                    let a = Base::from(seq[j].unwrap());
                    let to_write = if seqs.ploidy == Ploidy::Diploid {
                        sc(a).to_string()
                    } else if !fsnp[j] {
                        a.to_string()
                    } else if fsnp[j] && a == Base::N {
                        "?".to_string()
                    } else {
                        (a as usize).to_string()
                    };
                    ofp.write(to_write.as_bytes())?;
                    na += 1;
                    if (na % 50) == 0 {
                        ofp.write(b"\n")?;
                    }
                }
            }
            if (na % 50) != 0 {
                ofp.write(b"\n")?;
            }
        }
        for (i, &o) in &mut output_site.iter().enumerate() {
            if o {
                write!(loc, "\n{:.3}", locs.data[i])?;
            }
        }
        Ok(())
    }
}
