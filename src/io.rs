use crate::{scan, Error};
use bio::io::fasta;
use polars::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;

#[derive(Debug, PartialEq)]
pub struct Locs {
    pub data: Vec<f64>,
    length: f64,
    model: Model,
}

/// Model of crossing-over or gene conversion
#[derive(Debug, PartialEq)]
enum Model {
    CrossingOver,
    GeneConversion,
}

impl std::str::FromStr for Model {
    type Err = crate::Error;
    fn from_str(s: &str) -> Result<Model, Self::Err> {
        match s {
            "L" => Ok(Model::CrossingOver),
            "C" => Ok(Model::GeneConversion),
            _ => Err(Self::Err::new("Variant not found")),
        }
    }
}

#[test]
fn test_parse_model() {
    use std::str::FromStr;
    assert_eq!(Model::from_str("L").unwrap(), Model::CrossingOver);
    assert_eq!("L".parse::<Model>().unwrap(), Model::CrossingOver);
    assert_eq!(Model::from_str("C").unwrap(), Model::GeneConversion);
    assert_eq!("C".parse::<Model>().unwrap(), Model::GeneConversion);
}

pub fn read_locs(path: &PathBuf) -> Result<Locs, Box<dyn std::error::Error>> {
    let mut file = File::open(path)?;
    let mut content = String::new();
    file.read_to_string(&mut content)?;
    parse_locs(content.as_str())
}

fn parse_locs(content: &str) -> Result<Locs, Box<dyn std::error::Error>> {
    let mut lines = content.lines();
    let first_line = lines.next().ok_or(crate::Error::new("Cannot parse the first line"))?;
    let (l, length, model) = scan!(first_line, char::is_whitespace, usize, f64, Model);
    let content = lines.collect::<Vec<&str>>().join("\n");
    let mut data = Vec::<f64>::new();
    for loc in content.split_whitespace() {
        let x = loc.parse::<f64>()?;
        if data.last().unwrap_or(&0.) >= &x {
            return Err(Box::new(Error::new(
                "loc file SNPs not monotonically increasing",
            )));
        }
        data.push(x);
    }
    assert_eq!(l, data.len());
    Ok(Locs {
        data,
        length,
        model,
    })
}

#[test]
fn test_parse_locs() {
    let locs_content = r#"10 1200 L
1 57 180 187 223 250 438 509 878 1034"#;
    let locs = parse_locs(locs_content).unwrap();
    assert_eq!(
        locs,
        Locs {
            data: vec![1., 57., 180., 187., 223., 250., 438., 509., 878., 1034.],
            length: 1200.,
            model: Model::CrossingOver
        }
    );
}

pub struct Seqs {
    pub data: polars::frame::DataFrame,
}

pub fn read_fasta(path: &PathBuf) -> Result<Seqs, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut buf = String::new();
    reader.read_line(&mut buf)?;
    let (_nseq, _lseq, hd) = scan!(buf, char::is_whitespace, usize, usize, u8);
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
