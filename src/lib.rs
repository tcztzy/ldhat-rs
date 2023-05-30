pub mod commands;
pub mod error;
pub mod io;
pub use error::Error;
pub use io::read_locs;
pub const SEQ_MAX: usize = 1000;

pub type LDhatResult<T> = anyhow::Result<T>;
