use std::time::SystemTime;
use rand::SeedableRng;

pub struct RandomSeed([u8; 8]);

impl Default for RandomSeed {
    fn default() -> Self {
        Self(unix_timestamp().to_be_bytes())
    }
}

impl AsMut<[u8]> for RandomSeed {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

pub struct Random(RandomSeed);

impl SeedableRng for Random {
    type Seed = RandomSeed;
    fn from_seed(seed: Self::Seed) -> Self {
        Random(seed)
    }
}

/// Return UNIX timestamp in sec.
pub fn unix_timestamp() -> u64 {
    match SystemTime::now().duration_since(SystemTime::UNIX_EPOCH) {
        Ok(n) => n.as_secs(),
        Err(_) => panic!("SystemTime before UNIX EPOCH!"),
    }
}

#[macro_export]
macro_rules! scan {
    ( $string:expr, $sep:expr, $( $x:ty ),+ ) => {{
        let mut iter = $string.split($sep);
        ($(iter.next().and_then(|word| word.parse::<$x>().ok()).unwrap(),)*)
    }}
}

pub use scan;
