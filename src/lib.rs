#![feature(str_split_as_str)]
pub mod error;
pub mod io;
pub mod tools;
pub use error::Error;
pub use io::read_locs;
pub const SEQ_MAX: usize = 1000;
