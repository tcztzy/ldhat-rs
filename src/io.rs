use crate::LDhatResult as Result;
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
    fn from_str(s: &str) -> std::result::Result<Model, Self::Err> {
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

pub fn read_locs(path: &PathBuf) -> Result<Locs> {
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

fn parse_locs(content: &str) -> Result<Locs> {
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
        return Err(anyhow::anyhow!(
            "loc file SNPs not monotonically increasing"
        ));
    }
    Ok(Locs {
        data,
        length,
        model,
    })
}

#[test]
fn test_parse_locs() {
    let locs_content = r#"10  1200 L
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
    ploidy: Ploidy,
    data: polars::frame::DataFrame,
}

impl Seqs {
    pub fn names(&self) -> Vec<&str> {
        self.data.get_column_names()
    }
    pub fn shape(&self) -> (usize, usize) {
        self.data.shape()
    }
}

impl std::ops::Index<&str> for Seqs {
    type Output = polars::series::Series;

    fn index(&self, idx: &str) -> &Self::Output {
        &self.data[idx]
    }
}

/// Ploidy (/ˈplɔɪdi/) is the number of complete sets of chromosomes in a cell,
/// and hence the number of possible alleles for autosomal and pseudoautosomal genes.
#[derive(Debug, PartialEq)]
enum Ploidy {
    Haploid = 1,
    Diploid = 2,
}

impl From<char> for Ploidy {
    fn from(value: char) -> Self {
        match value {
            '1' => Ploidy::Haploid,
            '2' => Ploidy::Diploid,
            _ => Err(crate::Error::new("Variant not found")).unwrap(),
        }
    }
}

pub fn read_sites(path: &PathBuf) -> Result<Seqs> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    parse_sites(&mut reader)
}

/// Backward compatible for original code
enum Base {
    N = 1,
    T,
    C,
    A,
    G,
}

fn parse_sites(reader: &mut impl BufRead) -> Result<Seqs> {
    use nom::character::complete::{digit1, newline, one_of, space1};
    use nom::combinator::map_res;
    use nom::sequence::terminated;
    let mut buf = String::new();
    reader.read_line(&mut buf)?;
    let first_line = buf.as_str();
    let (_, (nseq, lseq, ploidy)) = nom::sequence::tuple((
        terminated(
            map_res(digit1::<_, (_, nom::error::ErrorKind)>, str::parse::<usize>),
            space1,
        ),
        terminated(
            map_res(digit1::<_, (_, nom::error::ErrorKind)>, str::parse::<usize>),
            space1,
        ),
        terminated(
            map_res(one_of("12"), |c| Ok::<_, crate::Error>(Ploidy::from(c))),
            newline,
        ),
    ))(first_line)
    .unwrap();
    let reader = fasta::Reader::new(reader);
    let data = DataFrame::from_iter(reader.records().into_iter().map(|res| {
        let rec = res.expect("Error during fasta record parsing");
        let ret = Series::new(
            rec.id(),
            rec.seq()
                .iter()
                .map(|&x| match x {
                    48 | 84 | 116 => Base::T,
                    49 | 67 | 99 => Base::C,
                    50 | 65 | 97 => Base::A,
                    51 | 71 | 103 => Base::G,
                    _ => Base::N,
                } as u8)
                .collect::<Vec<u8>>(),
        );
        ret
    }));
    assert_eq!(data.shape(), (lseq, nseq));
    Ok(Seqs { ploidy, data })
}

#[test]
fn test_parse_sites() {
    let content = r#"4 10 1
>SampleA
TCCGC??RTT
>SampleB
TACGC??GTA
>SampleC
TC?-CTTGTA
>SampleD
TCC-CTTGTT"#;
    let mut reader = std::io::BufReader::new(content.as_bytes());
    let seqs = parse_sites(&mut reader).unwrap();
    assert_eq!(seqs.ploidy, Ploidy::Haploid);
    let sample_a = &seqs["SampleA"];
    println!("{}", sample_a.to_string());
    assert!(sample_a.series_equal(&polars::series::Series::new(
        "SampleA",
        [2u8, 3, 3, 5, 3, 1, 1, 1, 2, 2]
    )));
}
