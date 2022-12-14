use crate::{scan, Error};
use bio::io::fasta;
use polars::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

pub struct Locs {
    pub data: Vec<f64>,
    length: f64,
    c: char,
}
pub fn read_locs(path: &PathBuf) -> Result<Locs, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let mut buf = String::new();
    let mut reader = BufReader::new(file);
    reader.read_line(&mut buf)?;
    let (l, length, c) = scan!(buf, char::is_whitespace, usize, f64, char);

    let mut data = Vec::<f64>::new();
    for _ in 0..l {
        reader.read_line(&mut buf)?;
        let x = buf.parse::<f64>()?;
        if data.last().unwrap_or(&0.) >= &x {
            return Err(Box::new(Error::new(
                "loc file SNPs not monotonically increasing",
            )));
        }
        data.push(x);
    }
    assert_eq!(l, data.len());
    Ok(Locs { data, length, c })
}

pub struct Seqs {
    pub data: polars::frame::DataFrame,
}

pub fn read_fasta(path: &PathBuf) -> Result<Seqs, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut buf = String::new();
    reader.read_line(&mut buf)?;
    let (mut nseq, lseq, hd) = scan!(buf, char::is_whitespace, usize, usize, u8);
    if hd < 1 || hd > 2 {
        eprintln!("The first line of seqs file should have N_Seqs L_Seq Haploid(1)/Diploid(2)");
        return Err(Box::new(Error::new("have not read haploid/diploid status")));
    }
    let reader = fasta::Reader::new(reader);
    let data = DataFrame::from_iter(reader.records().into_iter().map(|res| {
        let rec = res.expect("Error during fasta record parsing");
        let ret = Series::new(rec.id(), rec.seq());
        ret
    }));
    Ok(Seqs { data })
}
