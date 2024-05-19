use std::{
    fmt::{self, Display},
    num::ParseIntError,
    ops::RangeInclusive,
};
use crate::{FcsError, VALID_FCS_VERSIONS, Seek};
use std::io::BufReader;
use std::str;
use std::io::Read;

/// The header segment with its offsets.
///
/// Offsets are zero if the segment doesn't exist (only applies to analysis) or
/// if the offsets don't fit in the header and are instead written to the
/// text segment (applies to all).
#[derive(Debug, PartialEq)]
pub struct Header {
    pub version: String,
    pub text_offsets: RangeInclusive<usize>,
    pub data_offsets: RangeInclusive<usize>,
    pub analysis_offsets: RangeInclusive<usize>,
}

/// Attempts to create a `Header` struct from a string slice representing the FCS file header.
///
/// This implementation parses specific segments of the header string to extract version information
/// and various offsets. The offsets indicate the start and end positions of the text, data, and analysis
/// segments within the FCS file.
///
/// # Arguments
///
/// * `header` - A string slice containing the header segment of the FCS file.
///
/// # Returns
///
/// A `Result` containing a `Header` struct if the parsing is successful, or a `ParseIntError` if
/// any of the numeric parsing operations fail.
///
/// # Errors
///
/// This function will return a `ParseIntError` if any of the required numeric fields in the header
/// string cannot be parsed.
///
/// # Examples
///
/// ```
/// use std::convert::TryFrom;
/// use fcs_rs::header::Header;
///
/// let header_str = "FCS3.0         256    1545       0       0       0       0";
/// let header = Header::try_from(header_str).unwrap();
/// println!("{:?}", header);
/// ```
///
/// # Fields
///
/// The header string is expected to be formatted as follows:
/// - Bytes 0-5: Version (e.g., "FCS3.0")
/// - Bytes 10-17: Text segment start offset
/// - Bytes 18-25: Text segment end offset
/// - Bytes 26-33: Data segment start offset
/// - Bytes 34-41: Data segment end offset
/// - Bytes 42-49: Analysis segment start offset (optional, defaults to 0 if empty)
/// - Bytes 50-57: Analysis segment end offset (optional, defaults to 0 if empty)
impl TryFrom<&str> for Header {
    type Error = ParseIntError;

    fn try_from(header: &str) -> Result<Self, Self::Error> {
        let version = header[0..=5].to_string();
        let text_start = header[10..=17].trim_start().parse::<usize>()?;
        let text_end = header[18..=25].trim_start().parse::<usize>()?;
        let data_start = header[26..=33].trim_start().parse::<usize>()?;
        let data_end = header[34..=41].trim_start().parse::<usize>()?;
        let analysis_start = header[42..=49]
            .trim_start()
            .parse::<usize>()
            .or_else(zero_when_empty)?;
        let analysis_end = header[50..=57]
            .trim_start()
            .parse::<usize>()
            .or_else(zero_when_empty)?;

        Ok(Self {
            version,
            text_offsets: text_start..=text_end,
            data_offsets: data_start..=data_end,
            analysis_offsets: analysis_start..=analysis_end,
        })
    }
}

impl Display for Header {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:<10}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            self.version,
            self.text_offsets.start(),
            self.text_offsets.end(),
            self.data_offsets.start(),
            self.data_offsets.end(),
            self.analysis_offsets.start(),
            self.analysis_offsets.end(),
        )
    }
}

/// A helper function to return zero when the parsing of an empty string fails.
///
/// This function is used to handle optional analysis segment offsets in the FCS header.
///
/// # Arguments
///
/// * `err` - A `ParseIntError` indicating the failure to parse an empty string.
///
/// # Returns
///
/// Always returns `Ok(0)` indicating that the offset should default to zero.
fn zero_when_empty(err: ParseIntError) -> Result<usize, ParseIntError> {
    if err.kind() == &std::num::IntErrorKind::Empty {
        Ok(0)
    } else {
        Err(err)
    }
}

/// Reads the header segment of the FCS file and returns a Header struct.
///
/// The header segment is expected to be 58 bytes long and contains important
/// metadata about the FCS file. This function reads the header, parses it, and
/// validates the FCS version.
///
/// # Arguments
///
/// * `reader` - A mutable reference to a BufReader wrapping the FCS file.
///
/// # Returns
///
/// A Result containing a Header struct or an FcsError.
///
/// # Errors
///
/// This function will return an FcsError if:
/// - There is an I/O error during reading.
/// - The header cannot be parsed as a UTF-8 string.
/// - The header cannot be converted to a Header struct.
/// - The FCS version in the header is not supported.
///
/// # Examples
///
/// ```
/// use std::fs::File;
/// use std::io::BufReader;
/// use std::io::Seek;
/// use fcs_rs::header::read_header;
/// use fcs_rs::{FcsError, VALID_FCS_VERSIONS};
/// 
/// let file = File::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
/// let mut reader = BufReader::new(&file);
/// let header = read_header(&mut reader).unwrap();
/// println!("{:?}", header);
/// ```
pub fn read_header<R: Read + Seek>(reader: &mut BufReader<R>) -> Result<Header, FcsError> {
    let mut buffer = [0u8; 58];
    reader.read_exact(&mut buffer).map_err(FcsError::IoError)?;
    let header_line = str::from_utf8(&buffer).map_err(|_| FcsError::InvalidHeader)?;

    let header = Header::try_from(header_line).map_err(|_| FcsError::InvalidHeader)?;

    let version = &header.version;
    if !VALID_FCS_VERSIONS.contains(&version.as_str()) {
        return Err(FcsError::InvalidVersion(version.to_owned()));
    }

    Ok(header)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header1() {
        let header =
            "FCS3.0         256    1545    1792  202456       0       0";
        let parsed_header = Header {
            version: "FCS3.0".into(),
            text_offsets: 256..=1545,
            data_offsets: 1792..=202456,
            analysis_offsets: 0..=0,
        };

        assert_eq!(parsed_header, Header::try_from(header).unwrap());
    }

    #[test]
    fn header1_spaces() {
        let header =
            "FCS3.0         256    1545    1792  202456                ";
        let parsed_header = Header {
            version: "FCS3.0".into(),
            text_offsets: 256..=1545,
            data_offsets: 1792..=202456,
            analysis_offsets: 0..=0,
        };

        assert_eq!(parsed_header, Header::try_from(header).unwrap());
    }

    #[test]
    fn header2() {
        let header =
            "FCS3.0         256    1545       0       0       0       0";
        let parsed_header = Header {
            version: "FCS3.0".into(),
            text_offsets: 256..=1545,
            data_offsets: 0..=0,
            analysis_offsets: 0..=0,
        };

        assert_eq!(parsed_header, Header::try_from(header).unwrap());
    }

    #[test]
    fn header3() {
        let header =
            "FCS3.0      202451  203140    1792  202450       0       0";
        let parsed_header = Header {
            version: "FCS3.0".into(),
            text_offsets: 202451..=203140,
            data_offsets: 1792..=202450,
            analysis_offsets: 0..=0,
        };

        assert_eq!(parsed_header, Header::try_from(header).unwrap());
    }

    #[test]
    fn write_header1() {
        let header = Header {
            version: "FCS3.0".into(),
            text_offsets: 256..=1545,
            data_offsets: 1792..=202456,
            analysis_offsets: 0..=0,
        };
        let formatted = header.to_string();

        assert_eq!(
            "FCS3.0         256    1545    1792  202456       0       0",
            formatted
        );
    }

    #[test]
    fn write_header2() {
        let header = Header {
            version: "FCS3.0".into(),
            text_offsets: 256..=1545,
            data_offsets: 0..=0,
            analysis_offsets: 0..=0,
        };
        let formatted = header.to_string();

        assert_eq!(
            "FCS3.0         256    1545       0       0       0       0",
            formatted
        );
    }

    #[test]
    fn write_header3() {
        let header = Header {
            version: "FCS3.0".into(),
            text_offsets: 202451..=203140,
            data_offsets: 1792..=202450,
            analysis_offsets: 0..=0,
        };
        let formatted = header.to_string();

        assert_eq!(
            "FCS3.0      202451  203140    1792  202450       0       0",
            formatted
        );
    }
}
