use crate::{scan, Error};
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
