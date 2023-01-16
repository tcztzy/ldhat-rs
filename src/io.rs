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

impl From<char> for Model {
    fn from(value: char) -> Self {
        match value {
            'L' => Model::CrossingOver,
            'C' => Model::GeneConversion,
            _ => Err(crate::Error::new("Variant not found")).unwrap(),
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

fn is_monotonic_increasing(arr: &[f64]) -> bool {
    let mut last = 0.;
    for &i in arr {
        if i < last {
            return false;
        }
        last = i
    }
    true
}

fn parse_locs(content: &str) -> Result<Locs, Box<dyn std::error::Error>> {
    use nom::{
        character::complete::{digit1, multispace1, newline, one_of, space1},
        combinator::map_res,
        multi::separated_list0,
        number::complete::double,
        sequence::{terminated, tuple},
    };
    let (content, (l, length, model, data)) = tuple((
        terminated(
            map_res(digit1::<_, (_, nom::error::ErrorKind)>, str::parse::<usize>),
            space1,
        ),
        terminated(double, space1),
        terminated(
            map_res(one_of("LC"), |c| Ok::<_, crate::Error>(Model::from(c))),
            newline,
        ),
        separated_list0(multispace1, double),
    ))(content)
    .unwrap();
    assert_eq!(content.len(), 0);
    assert_eq!(l, data.len());
    if !is_monotonic_increasing(&data) {
        return Err(Box::new(Error::new(
            "loc file SNPs not monotonically increasing",
        )));
    }
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
