use crate::{
    FcsError, 
    HashMap, 
    BufReader, 
    SeekFrom, 
    read_header, 
    Seek, 
    Read, 
    REQUIRED_KEYWORDS, 
    ReadBytesExt,
    File
};

/// Reads the text segment of the FCS file and returns a HashMap containing metadata.
///
/// The text segment contains key-value pairs of metadata information about the FCS file.
/// This function reads the text segment, parses it, and validates the extracted metadata.
///
/// # Arguments
///
/// * `reader` - A mutable reference to a BufReader wrapping the FCS file.
///
/// # Returns
///
/// A Result containing a HashMap of metadata key-value pairs or an FcsError.
///
/// # Errors
///
/// This function will return an FcsError if:
/// - There is an I/O error during reading.
/// - The text segment cannot be converted to a UTF-8 string.
/// - The metadata validation fails.
///
/// # Examples
///
/// ```
/// use fcs_rs::text::read_metadata;
/// use std::fs::File;
/// use std::io::BufReader;
/// 
/// let file = File::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
/// let mut reader = BufReader::new(&file);
/// let metadata = read_metadata(&mut reader).unwrap();
/// println!("{:?}", metadata);
/// ```
pub fn read_metadata(reader: &mut BufReader<&File>) -> Result<HashMap<String, String>, FcsError> {
    let header = read_header(reader)?;
    let text_offset = header.text_offsets;
    let mut metadata: HashMap<String, String> = HashMap::new();

    let bytes_to_read = text_offset.end() - text_offset.start() - 1;
    let mut buffer = vec![0u8; bytes_to_read];
    reader.seek(SeekFrom::Start(*text_offset.start() as u64)).map_err(FcsError::IoError)?;

    let delimiter = reader.read_u8().map_err(FcsError::IoError)? as char;
    reader.read_exact(&mut buffer).map_err(FcsError::IoError)?;

    let text = String::from_utf8(buffer).map_err(|_| FcsError::InvalidMetadata)?;

    let mut keyword = String::new();
    let mut value = String::new();
    let kv_pairs = text.split(delimiter);

    for kv in kv_pairs {
        if kv.starts_with('$') {
            keyword = kv.to_string();
            value.clear();
        } else {
            value.push_str(kv);
            metadata.insert(keyword.clone(), value.clone());
        }
    }

    validate_text(&metadata)?;
    Ok(metadata)
}

/// Validates that the required keys are present in the text segment of the FCS file.
///
/// This function checks for the presence of required metadata keys in the text segment.
/// If any required key is missing, it returns an `FcsError::InvalidText`.
///
/// # Arguments
///
/// * `text` - A reference to a HashMap containing metadata key-value pairs from the FCS file.
///
/// # Returns
///
/// A Result indicating success or an `FcsError`.
///
/// # Errors
///
/// This function will return an `FcsError::InvalidText` if any required metadata key is missing.
/// It will also return `FcsError::InvalidText` if the `$PAR` key cannot be parsed or is missing.
///
/// # Examples
///
/// ```
/// use std::collections::HashMap;
/// use fcs_rs::text::{read_metadata, validate_text};
/// use fcs_rs::REQUIRED_KEYWORDS;
/// use std::fs::File;
/// use std::io::BufReader;
/// use fcs_rs::FcsError;
/// use fcs_rs::header::read_header;
/// use std::io::SeekFrom;
/// 
/// let file = File::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
/// let mut reader = BufReader::new(&file);
/// let metadata: HashMap<String, String> = read_metadata(&mut reader).unwrap();
/// validate_text(&metadata).unwrap();
/// ```
pub fn validate_text(text: &HashMap<String, String>) -> Result<(), FcsError> {
    let n_params: u32 = text.get("$PAR")
        .ok_or_else(|| FcsError::InvalidText("$PAR".to_string()))?
        .parse()
        .map_err(|_| FcsError::InvalidText("$PAR".to_string()))?;

    for &non_param in &REQUIRED_KEYWORDS[..12] {
        if text.get(non_param).is_none() {
            return Err(FcsError::InvalidText(non_param.to_string()));
        }
    }

    for &param in &REQUIRED_KEYWORDS[12..] {
        for i in 1..=n_params {
            if text.get(&param.replace('n', &i.to_string())).is_none() {
                return Err(FcsError::InvalidText(param.to_string()));
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FcsFile;
    use std::io::BufReader;

    #[test]
    fn test_read_metadata_success() {
        let file = File::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
        let mut reader = BufReader::new(&file);

        let result = read_metadata(&mut reader);
        assert!(result.is_ok(), "Metadata read failed: {:?}", result.err());
        let metadata = result.unwrap();
        assert_eq!(metadata.get("$PAR").unwrap(), "10");
        assert_eq!(metadata.get("$P1N").unwrap(), "Time");
        assert_eq!(metadata.get("$P2N").unwrap(), "FSC-A");
    }

    #[test]
    fn test_validate_text_success() {
        let file = File::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
        let mut reader = BufReader::new(&file);

        let result = read_metadata(&mut reader);
        assert!(result.is_ok(), "Metadata read failed: {:?}", result.err());
        let metadata = result.unwrap();

        let result = validate_text(&metadata);
        assert!(result.is_ok(), "Text validation failed: {:?}", result.err());
    }

    #[test]
    fn test_dataframes_columns() {
        let fcs_file = FcsFile::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();

        let mut flow_sample = fcs_file.read().unwrap();
        // Get the column names using the function
        let keys_vec = flow_sample.get_dataframe_columns();
        
        // Perform arcsinh transformation of data
        match flow_sample.arcsinh_transform(5.0, &keys_vec) {
            Ok(_) => (),
            Err(err) => panic!("{}", err)
        };
    }
}
