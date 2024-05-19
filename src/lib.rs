//! Fcs File manipulation operations.
//!
//! This module contains basic methods to manipulate the contents of fcs files.

use std::collections::HashMap;
use std::io::{self, Read, Seek};
use std::io::{BufReader, SeekFrom};
use std::fs::File;
use std::str;
use byteorder::ReadBytesExt;
use thiserror::Error;

pub use crate::header::read_header;
pub use crate::text::{read_metadata, validate_text};
pub use crate::data::{FlowSample, parse_data, read_events, create_dataframe};

pub mod data;
pub mod header;
pub mod text;

pub const VALID_FCS_VERSIONS: [&str; 2] = ["FCS3.0", "FCS3.1"];

/// Required non-parameter indexed keywords for FCS text segment
pub const REQUIRED_KEYWORDS: [&str; 16] = [
    "$BEGINANALYSIS", // Byte-offset to the beginning of the ANALYSIS segment
    "$BEGINDATA", // Byte-offset to the beginning of the DATA segment
    "$BEGINSTEXT", // Byte-offset to the beginning of a supplemental TEXT segment
    "$BYTEORD", // Byte order for data acquisition computer
    "$DATATYPE", // Type of data in DATA segment (ASCII, integer, floating point)
    "$ENDANALYSIS", // Byte-offset to the last byte of the ANALYSIS segment
    "$ENDDATA", // Byte-offset to the last byte of the DATA segment
    "$ENDSTEXT", // Byte-offset to the last byte of a supplemental TEXT segment
    "$MODE", // Data mode (list mode - preferred, histogram - deprecated)
    "$NEXTDATA", // Byte offset to next data set in the file
    "$PAR", // Number of parameters in an event
    "$TOT", // Total number of events in the dataset
    "$PnB", // Number of bits reserved for parameter number n.
    "$PnE", // Amplification type for parameter n.
    "$PnN", // Short name for parameter n.
    "$PnR" // Range for parameter number n.
];

/// Represents errors that can occur while processing FCS (Flow Cytometry Standard) files.
///
/// This enum covers various error types that might be encountered, including I/O errors,
/// invalid headers, unsupported FCS versions, metadata issues, missing required keywords,
/// and invalid data segments.
///
/// # Variants
///
/// - `IoError`: Represents an I/O error that occurs during file operations.
/// - `InvalidHeader`: Indicates that the FCS header is invalid, possibly due to file corruption
///   or the file not being a valid FCS file.
/// - `InvalidVersion`: Indicates that the FCS version is not supported. The supported versions are
///   FCS3.0 and FCS3.1.
/// - `InvalidMetadata`: Indicates that the FCS metadata is invalid.
/// - `InvalidText`: Indicates that the FCS file is missing a required keyword in its TEXT section.
/// - `InvalidData`: Indicates that the FCS data segment is invalid, with an associated error message.
///
/// # Examples
///
/// ```
/// use thiserror::Error;
/// use std::io;
///
/// #[derive(Debug, Error)]
/// pub enum FcsError {
///     #[error("IO Error: {0}")]
///     IoError(#[from] io::Error),
///     #[error("Invalid FCS Header. File may be corrupted or not a FCS file.")]
///     InvalidHeader,
///     #[error("FCS version `{0}` not supported. Must be either FCS3.0 or FCS3.1")]
///     InvalidVersion(String),
///     #[error("Invalid FCS Metadata")]
///     InvalidMetadata,
///     #[error("FCS file is corrupted. It is missing required keyword {0} in its TEXT section")]
///     InvalidText(String),
///     #[error("Invalid FCS Data: {0}")]
///     InvalidData(String),
/// }
/// ```
#[derive(Debug, Error)]
pub enum FcsError {
    #[error("IO Error: {0}")]
    IoError(#[from] io::Error),
    #[error("Invalid FCS Header. File may be corrupted or not a FCS file.")]
    InvalidHeader,
    #[error("FCS version `{0}` not supported. Must be either FCS3.0 or FCS3.1")]
    InvalidVersion(String),
    #[error("Invalid FCS Metadata")]
    InvalidMetadata,
    #[error("FCS file is corrupted. It is missing required keyword {0} in its TEXT section")]
    InvalidText(String),
    #[error("Invalid FCS Data: {0}")]
    InvalidData(String),
}

/// An object providing access to an FCS file.
///
/// This struct wraps a file handle and provides methods to open the file and read
/// metadata and parameter data from it.
pub struct FcsFile {
    inner: File,
}

impl FcsFile {
    /// Open an FCS file in read-only mode.
    ///
    /// # Arguments
    ///
    /// * `path` - A string slice representing the path to the FCS file.
    ///
    /// # Returns
    ///
    /// A `Result` containing an `FcsFile` object if the file is successfully opened,
    /// or an `FcsError` if there is an issue opening the file.
    ///
    /// # Examples
    ///
    /// ```
    /// use no_bs_fcs::FcsFile;
    /// 
    /// let fcs_file = FcsFile::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
    /// ```
    pub fn open(path: &str) -> Result<FcsFile, FcsError> {
        let file = File::open(path).map_err(FcsError::IoError)?;

        Ok(Self { inner: file })
    }

    /// Create an FcsFile from an existing File object.
    ///
    /// # Arguments
    ///
    /// * `file` - An already opened File object.
    ///
    /// # Returns
    ///
    /// An `FcsFile` that wraps the provided File object.
    ///
    /// # Examples
    ///
    /// ```
    /// use no_bs_fcs::FcsFile;
    /// use std::fs::File;
    /// 
    /// let file = File::open("path/to/file.fcs").unwrap();
    /// let fcs_file = FcsFile::from_file(file);
    /// ```
    pub fn from_file(file: File) -> Self {
        Self { inner: file }
    }

    /// Read the FCS file and return metadata and parameter data in an `FlowSample` struct.
    ///
    /// # Returns
    ///
    /// A `Result` containing an `FlowSample` struct if the file is successfully read,
    /// or an `FcsError` if there is an issue reading the metadata or parameter data.
    ///
    /// # Examples
    ///
    /// ```
    /// use no_bs_fcs::{FcsFile, FcsError, REQUIRED_KEYWORDS};
    /// use no_bs_fcs::header::read_header;
    /// use no_bs_fcs::data::{FlowSample, parse_data, read_events};
    /// use no_bs_fcs::text::{read_metadata, validate_text};
    /// use std::collections::HashMap;
    /// use std::fs::File;
    /// use std::io::BufReader;
    /// use std::io::SeekFrom;
    /// 
    /// let fcs_file = FcsFile::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
    /// let fcs_data = fcs_file.read().unwrap();
    /// println!("{:?}", fcs_data.data);
    /// println!("{:?}", fcs_data.parameters);
    /// ```
    pub fn read(&self) -> Result<FlowSample, FcsError> {
        let mut reader = BufReader::new(&self.inner);
        let metadata = match read_metadata(&mut reader) {
            Ok(metadata) => metadata,
            Err(err) => panic!("There was an error: {}", err)
        };
    
        let flow_sample = match parse_data(&mut reader, &metadata) {
            Ok(sample) => sample,
            Err(err) => panic!("Error parsing Flow sample: {}", err)
        };

        Ok(flow_sample)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fcs_open() {
        let result = FcsFile::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs");
        assert!(result.is_ok(), "FCS file open failed: {:?}", result.err());
    }

    #[test]
    fn test_fcs_open_non_existent_file() {
        let result = FcsFile::open("./examples/non_existent.fcs");
        assert!(result.is_err(), "Opening a non-existent file should fail");
    }

    #[test]
    fn test_fcs_read() {
        let fcs_file = FcsFile::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
        let result = fcs_file.read();
        assert!(result.is_ok(), "FCS file read failed: {:?}", result.err());
    }

    #[test]
    fn test_create_dataframe() {
        let column_titles = vec!["col1".to_string(), "col2".to_string()];
        let data = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let df = data::create_dataframe(&column_titles, &data).unwrap();
        assert_eq!(df.shape(), (3, 2), "DataFrame shape mismatch");
    }
}
