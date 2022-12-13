#[derive(Debug, Clone)]
pub struct Error {
    pub message: String,
}

impl std::error::Error for Error {}

impl Error {
    pub fn new(message: &str) -> Self {
        Self {
            message: message.to_string(),
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}
