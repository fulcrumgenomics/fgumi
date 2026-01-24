use std::{fmt, io, result};

use crate::{CompCode, DecompCode};

pub type Result<T> = result::Result<T, Error>;

#[derive(Debug)]
pub enum Error {
    Io(io::Error),

    // isal compression errors
    CompressionError(CompCode),

    // isal decompression errors
    DecompressionError(DecompCode),

    // Anything else not covered, exit code, message
    Other((Option<isize>, String)),
}

impl From<Error> for io::Error {
    fn from(err: Error) -> io::Error {
        io::Error::new(io::ErrorKind::Other, err)
    }
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Error {
        Error::Io(err)
    }
}

impl std::error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // write!(f, msg)
        match self {
            Error::Io(err) => write!(f, "{err}"),
            Error::Other((exit_code, msg)) => write!(f, "{msg} (Exit code: {exit_code:?})"),
            _ => write!(f, "{self:?}"),
        }
    }
}
