use byteorder::{
    BigEndian, LittleEndian,
};
use std::fmt;
use std::io::{Read, Seek};
use crate::{FcsError, HashMap, File, BufReader, SeekFrom};
use polars::prelude::*;

/// Store the names of the parameters in the FCS file.
#[derive(Debug)]
pub struct ColumnNames {
    keys_vec: Vec<String>,
}

/// A structure to represent a flow cytometry sample.
///
/// # Fields
///
/// * `data` - A HashMap containing the data for each channel. The key is the channel name, 
/// and the value is a vector of measurements.
/// * `parameters` - A HashMap containing the parameters of the flow cytometry experiment. 
/// The key is the parameter name, and the value is the parameter value.
#[derive(Debug)]
pub struct FlowSample {
    pub data: DataFrame,
    pub parameters: HashMap<String, String>,
}

impl fmt::Display for FlowSample {
    /// Formats the `FlowSample` for display.
    ///
    /// The display includes general information about the sample such as machine type, 
    /// run times, and volume, as well as details about the measurement axes.
    ///
    /// # Arguments
    ///
    /// * `f` - The formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "
FlowSample:
    Machine: {}
    Begin Time: {}
    End Time: {}
    Date: {}
    File: {}
    Volume run: {}",
            self.parameters.get("$CYT").unwrap_or(&"Unknown".to_string()),
            self.parameters.get("$BTIM").unwrap_or(&"Unknown".to_string()),
            self.parameters.get("$ETIM").unwrap_or(&"Unknown".to_string()),
            self.parameters.get("$DATE").unwrap_or(&"Unknown".to_string()),
            self.parameters.get("$FIL").unwrap_or(&"Unknown".to_string()),
            self.parameters.get("$VOL").unwrap_or(&"Unknown".to_string())
        )?;

        let n_params: u32 = self.parameters.get("$PAR").unwrap_or(&"0".to_string()).parse().unwrap_or(0);
        writeln!(f, "\n    Labels: ")?;
        for i in 1..=n_params {
            if self.parameters.contains_key(&format!("$P{}S", i)) {
                writeln!(f, "        {} ({})", 
                    self.parameters.get(&format!("$P{}N", i)).unwrap(), 
                    self.parameters.get(&format!("$P{}S", i)).unwrap()
                )?
            }
        }
        Ok(())
    }
}

impl FlowSample {
    /// Extracts and returns the column names from a FlowSample as a vector of strings.
    ///
    /// This function retrieves the column names from the `data` field of the provided
    /// `FlowSample` and converts them to a vector of strings.
    ///
    /// # Arguments
    ///
    /// * `flow_sample` - A reference to a FlowSample from which to extract column names.
    ///
    /// # Returns
    ///
    /// A vector of strings representing the column names.
    ///
    /// # Examples
    ///
    /// ```
    /// use fcs_rs::{FcsFile, FcsError};
    /// use fcs_rs::data::ColumnNames;
    /// use polars::prelude::*;
    /// use std::fmt;
    /// 
    /// let fcs_file = FcsFile::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
    /// let flow_sample = fcs_file.read().unwrap();
    /// let column_names = flow_sample.get_dataframe_columns();
    /// println!("{:?}", column_names);
    /// ```
    pub fn get_dataframe_columns(&self) -> Vec<String> {
        let column_names: Vec<String> = self.data.get_column_names().iter().map(|s| s.to_string()).collect();
        let keys_vec: Vec<String> = column_names.clone();

        let names = ColumnNames {
            keys_vec,
        };

        names.keys_vec
    }

    /// Applies the Arcsinh transformation to the data of specified channels.
    ///
    /// The Arcsinh transform is a combination of logarithmic and linear scales.
    ///
    /// # Arguments
    ///
    /// * `cofactor` - A scaling factor applied to the data before transformation.
    /// * `channels` - A vector of channel names to which the transformation will be applied.
    ///
    /// # Returns
    ///
    /// A Result indicating success or an I/O error.
    pub fn arcsinh_transform(
        &mut self, 
        cofactor: f64, 
        channels: &[String]
    ) -> Result<(), Box<dyn std::error::Error>> {
        fn arcsinh(x: f64, cofactor: f64) -> f64 {
            (x / cofactor).ln() + ((x / cofactor).powi(2) + 1.0).sqrt().ln()
        }

        for channel in channels.iter() {
            match self.data.column(channel) {
                Ok(series) => {
                    let transformed_series = series.f64()?
                    .apply(|v| Some(arcsinh(v.expect("REASON"), cofactor)))
                        .into_series();
                    self.data.with_column(transformed_series)?;
                },
                Err(err) => {
                    eprintln!("Error during arcsinh transformation: {}", err)
                }
            }
        }

        Ok(())
    }
}

/// Reads the data segment of the FCS file and returns a FlowSample struct.
///
/// # Arguments
///
/// * `reader` - A mutable reference to a BufReader wrapping the FCS file.
/// * `metadata` - A reference to a HashMap containing metadata from the FCS file.
///
/// # Returns
///
/// A Result containing a FlowSample struct or an FcsError.
///
/// # Errors
///
/// This function will return an FcsError if:
/// - The data mode is not 'L' (list mode).
/// - Required metadata fields are missing or invalid.
/// - The data segment cannot be read or parsed correctly.
/// - Byte order determination fails.
///
/// # Examples
///
/// ```
/// use std::fs::File;
/// use std::io::{self, Read, Seek, BufReader};
/// use std::collections::HashMap;
/// use std::io::SeekFrom;
/// use byteorder::{ReadBytesExt, LittleEndian, BigEndian};
/// use fcs_rs::data::{parse_data, FlowSample};
/// use fcs_rs::{FcsFile, FcsError, REQUIRED_KEYWORDS};
/// use fcs_rs::text::{read_metadata, validate_text};
/// use fcs_rs::header::read_header;
/// 
/// let file = File::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
/// let mut reader = BufReader::new(&file);
/// let metadata = read_metadata(&mut reader).unwrap();
/// let flow_sample = parse_data(&mut reader, &metadata).unwrap();
/// println!("{:?}", flow_sample.data);
/// ```
pub fn parse_data(
    reader: &mut BufReader<&File>, 
    metadata: &HashMap<String, String>,
) -> Result<FlowSample, FcsError> {
    let mode = metadata.get("$MODE")
        .ok_or_else(|| FcsError::InvalidData("Missing $MODE in metadata".to_string()))?;
    if mode != "L" {
        return Err(FcsError::InvalidData("Data must be in list (L) mode".to_string()));
    }

    let data_type = metadata.get("$DATATYPE")
        .ok_or_else(|| FcsError::InvalidData("Missing $DATATYPE in metadata".to_string()))?;
    let n_params = metadata.get("$PAR")
        .ok_or_else(|| FcsError::InvalidData("Missing $PAR in metadata".to_string()))?
        .parse::<usize>()
        .map_err(|_| FcsError::InvalidData("Invalid $PAR value".to_string()))?;
    let n_events = metadata.get("$TOT")
        .ok_or_else(|| FcsError::InvalidData("Missing $TOT in metadata".to_string()))?
        .parse::<usize>()
        .map_err(|_| FcsError::InvalidData("Invalid $TOT value".to_string()))?;
    let data_start = metadata.get("$BEGINDATA")
        .ok_or_else(|| FcsError::InvalidData("Missing $BEGINDATA in metadata".to_string()))?
        .trim()
        .parse::<u64>()
        .map_err(|_| FcsError::InvalidData("Invalid $BEGINDATA value".to_string()))?;
    let byte_order = metadata.get("$BYTEORD")
        .ok_or_else(|| FcsError::InvalidData("Missing $BYTEORD in metadata".to_string()))?;
    let capacity = n_params * n_events;
    if capacity == 0 {
        return Err(FcsError::InvalidData("Fcs file may be corrupted. No data found".to_string()));
    }

    reader.seek(SeekFrom::Start(data_start))?;
    let mut parameters: HashMap<String, Vec<f64>> = HashMap::new();
    let mut events: Vec<f64>;

    for i in 1..=n_params {
        events = if byte_order == "1,2,3,4" {
            read_events::<LittleEndian>(reader, data_type, n_events, i, metadata)?
        } else if byte_order == "4,3,2,1" {
            read_events::<BigEndian>(reader, data_type, n_events, i, metadata)?
        } else {
            return Err(FcsError::InvalidData("Could not determine byte order.".to_string()));
        };

        let id = metadata.get(&format!("$P{}S", i))
            .ok_or_else(|| FcsError::InvalidData(format!("Missing $P{}S in metadata", i)))?;
        parameters.insert(id.to_owned(), events);
    }

    let column_titles = parameters.keys().cloned().collect::<Vec<_>>();
    let data = parameters.values().cloned().collect::<Vec<_>>();

    let fcs_df = create_dataframe(&column_titles, &data)
        .map_err(|_| FcsError::InvalidData("Failed to create DataFrame".to_string()))?;

    let sample = FlowSample {
        data: fcs_df,
        parameters: metadata.to_owned()
    };

    Ok(sample)
}

/// Reads the events data segment of the FCS file and returns a vector of f64 values.
///
/// This function reads the events data based on the specified data type (F, D, or I) and
/// converts the raw data into a vector of f64 values. It handles different byte orders
/// and bit depths as specified in the metadata.
///
/// # Arguments
///
/// * `reader` - A mutable reference to a BufReader wrapping the FCS file.
/// * `data_type` - A string slice indicating the data type ('F' for float, 'D' for double, 'I' for integer).
/// * `n_events` - The number of events to read.
/// * `param_idx` - The parameter index for the current data segment.
/// * `metadata` - A reference to a HashMap containing metadata from the FCS file.
///
/// # Returns
///
/// A Result containing a vector of f64 values representing the events data, or an FcsError.
///
/// # Errors
///
/// This function will return an FcsError if:
/// - There is an I/O error during reading.
/// - The specified data type is not supported.
/// - The bits per parameter for integer data types is not supported or cannot be parsed.
///
/// # Examples
///
/// ```
/// use std::fs::File;
/// use std::io::BufReader;
/// use std::collections::HashMap;
/// use byteorder::LittleEndian;
/// use fcs_rs::text::read_metadata;
/// use fcs_rs::data::read_events;
/// 
/// let file = File::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
/// let mut reader = BufReader::new(&file);
/// let metadata: HashMap<String, String> = read_metadata(&mut reader).unwrap();
/// let events = read_events::<LittleEndian>(&mut reader, "F", 1000, 1, &metadata).unwrap();
/// println!("{:?}", events);
/// ```
pub fn read_events<B: byteorder::ByteOrder>(
    reader: &mut BufReader<&File>, 
    data_type: &str, 
    n_events: usize, 
    param_idx: usize, 
    metadata: &HashMap<String, String>
) -> Result<Vec<f64>, FcsError> {
    let data = match data_type {
        "F" => {
            let mut float_buffer = vec![0; n_events * std::mem::size_of::<f32>()];
            reader.read_exact(&mut float_buffer).map_err(FcsError::IoError)?;
            let mut data = Vec::with_capacity(n_events);
            for i in 0..n_events {
                let float_value = B::read_f32(&float_buffer[i * 4..(i + 1) * 4]) as f64;
                data.push(float_value);
            }
            data
        },
        "D" => {
            let mut float_buffer = vec![0; n_events * std::mem::size_of::<f64>()];
            reader.read_exact(&mut float_buffer).map_err(FcsError::IoError)?;
            let mut data = Vec::with_capacity(n_events);
            for i in 0..n_events {
                let double_value = B::read_f64(&float_buffer[i * 8..(i + 1) * 8]);
                data.push(double_value);
            }
            data
        },
        "I" => {
            let bits_per_param = metadata.get(&format!("$P{}B", param_idx))
                .ok_or_else(|| FcsError::InvalidText(format!("Missing $P{}B in metadata", param_idx)))?
                .parse::<usize>()
                .map_err(|_| FcsError::InvalidData(format!("Invalid bits per param value for $P{}B", param_idx)))?;

            match bits_per_param / 8 {
                2 => {
                    let mut int_buffer = vec![0; n_events * std::mem::size_of::<u16>()];
                    reader.read_exact(&mut int_buffer).map_err(FcsError::IoError)?;
                    let mut data = Vec::with_capacity(n_events);
                    for i in 0..n_events {
                        let int_value = B::read_u16(&int_buffer[i * 2..(i + 1) * 2]) as f64;
                        data.push(int_value);
                    }
                    data
                },
                4 => {
                    let mut int_buffer = vec![0; n_events * std::mem::size_of::<u32>()];
                    reader.read_exact(&mut int_buffer).map_err(FcsError::IoError)?;
                    let mut data = Vec::with_capacity(n_events);
                    for i in 0..n_events {
                        let int_value = B::read_u32(&int_buffer[i * 4..(i + 1) * 4]) as f64;
                        data.push(int_value);
                    }
                    data
                },
                8 => {
                    let mut int_buffer = vec![0; n_events * std::mem::size_of::<u64>()];
                    reader.read_exact(&mut int_buffer).map_err(FcsError::IoError)?;
                    let mut data = Vec::with_capacity(n_events);
                    for i in 0..n_events {
                        let int_value = B::read_u64(&int_buffer[i * 8..(i + 1) * 8]) as f64;
                        data.push(int_value);
                    }
                    data
                },
                16 => {
                    let mut int_buffer = vec![0; n_events * std::mem::size_of::<u128>()];
                    reader.read_exact(&mut int_buffer).map_err(FcsError::IoError)?;
                    let mut data = Vec::with_capacity(n_events);
                    for i in 0..n_events {
                        let int_value = B::read_u128(&int_buffer[i * 16..(i + 1) * 16]) as f64;
                        data.push(int_value);
                    }
                    data
                },
                _ => return Err(FcsError::InvalidData("Bits for param type not supported".to_string())),
            }
        },
        _ => return Err(FcsError::InvalidData("FCS data type not supported. Must be F, D, or I".to_string())),
    };

    Ok(data)
}

/// Creates a DataFrame from column titles and corresponding data vectors.
///
/// This function takes a vector of column titles and a vector of data vectors,
/// converts each data vector to a Series with the corresponding column title,
/// and combines them into a DataFrame.
///
/// # Arguments
///
/// * `column_titles` - A reference to a vector of strings representing the column titles.
/// * `data` - A reference to a vector of vectors containing the data for each column.
///
/// # Returns
///
/// A Result containing a DataFrame or a PolarsError if there is a mismatch in the number of columns and titles.
///
/// # Errors
///
/// This function will return a PolarsError if:
/// - The number of columns does not match the number of column titles.
///
/// # Examples
///
/// ```
/// use fcs_rs::data::create_dataframe;
/// use polars::prelude::*;
/// 
/// let column_titles = vec!["APC-A".to_string(), "FSC-W".to_string()];
/// let data = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
/// let df = create_dataframe(&column_titles, &data).unwrap();
/// println!("{:?}", df);
/// ```
pub fn create_dataframe(column_titles: &[String], data: &[Vec<f64>]) -> Result<DataFrame, PolarsError> {
    // Ensure the number of columns matches the number of column titles
    if column_titles.len() != data.len() {
        return Err(PolarsError::ShapeMismatch(
            "Number of columns does not match number of column titles".into(),
        ));
    }

    let mut series_vec = Vec::with_capacity(column_titles.len());

    for (i, column_data) in data.iter().enumerate() {
        let series = Series::new(&column_titles[i], column_data);
        series_vec.push(series);
    }

    let df = DataFrame::new(series_vec)?;

    Ok(df)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::BufReader;
    use crate::{read_metadata, FcsFile};

    #[test]
    fn test_flow_sample_display() {
        let mut parameters = HashMap::new();
        parameters.insert("$CYT".to_string(), "Test Cytometer".to_string());
        parameters.insert("$BTIM".to_string(), "10:00".to_string());
        parameters.insert("$ETIM".to_string(), "10:30".to_string());
        parameters.insert("$DATE".to_string(), "2022-01-01".to_string());
        parameters.insert("$FIL".to_string(), "test.fcs".to_string());
        parameters.insert("$VOL".to_string(), "500".to_string());
        parameters.insert("$PAR".to_string(), "2".to_string());
        parameters.insert("$P1N".to_string(), "FSC".to_string());
        parameters.insert("$P1S".to_string(), "Forward Scatter".to_string());
        parameters.insert("$P2N".to_string(), "SSC".to_string());
        parameters.insert("$P2S".to_string(), "Side Scatter".to_string());

        let data = DataFrame::new(vec![
            Series::new("FSC", &[1.0, 2.0, 3.0]),
            Series::new("SSC", &[4.0, 5.0, 6.0]),
        ]).unwrap();

        let flow_sample = FlowSample {
            data,
            parameters,
        };

        let expected_display = "
FlowSample:
    Machine: Test Cytometer
    Begin Time: 10:00
    End Time: 10:30
    Date: 2022-01-01
    File: test.fcs
    Volume run: 500
    Labels: 
        FSC (Forward Scatter)
        SSC (Side Scatter)
";

        assert_eq!(format!("{}", flow_sample), expected_display);
    }

    #[test]
    fn test_get_dataframe_columns() {
        let mut parameters = HashMap::new();
        parameters.insert("$PAR".to_string(), "2".to_string());
        parameters.insert("$P1N".to_string(), "FSC".to_string());
        parameters.insert("$P2N".to_string(), "SSC".to_string());

        let data = DataFrame::new(vec![
            Series::new("FSC", &[1.0, 2.0, 3.0]),
            Series::new("SSC", &[4.0, 5.0, 6.0]),
        ]).unwrap();

        let flow_sample = FlowSample {
            data,
            parameters,
        };

        let column_names = flow_sample.get_dataframe_columns();
        assert_eq!(column_names, vec!["FSC", "SSC"]);
    }

    #[test]
    fn test_arcsinh_transform() {
        let fcs_file = FcsFile::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();

        let mut flow_sample = fcs_file.read().unwrap();

        let column_names = flow_sample.get_dataframe_columns();
        
        assert!(flow_sample.arcsinh_transform(5.0, &column_names).is_ok(), "{}", false);   

        let fsc: Vec<f64> = flow_sample.data.column("FSC-H").unwrap().f64().unwrap().into_no_null_iter().collect();
        let ssc: Vec<f64> = flow_sample.data.column("APC-A").unwrap().f64().unwrap().into_no_null_iter().collect();

        assert_eq!(fsc[..=2], vec![22.692050373389968, 21.3191345559628, 5.869126409425455]);
        assert_eq!(ssc[..=2], vec![7.3780713106122935, 23.476587479533844, 23.01206351580942]);
    }

    #[test]
    fn test_parse_data() {
        let file = File::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
        let mut reader = BufReader::new(&file);

        let metadata = read_metadata(&mut reader).expect("Failed to read metadata");
        let flow_sample = parse_data(&mut reader, &metadata).expect("Failed to parse data");

        let fsc: Vec<f64> = flow_sample.data.column("FSC-H").unwrap().f64().unwrap().into_no_null_iter().collect();
        let ssc: Vec<f64> = flow_sample.data.column("APC-A").unwrap().f64().unwrap().into_no_null_iter().collect();

        assert_eq!(fsc[..=2], vec![423142.0, 212991.0, 94.0]);
        assert_eq!(ssc[..=2], vec![200.0, 626392.0, 496565.0]);
    }

    #[test]
    fn test_read_events() {
        let file = File::open("./examples/20200624 LEGENDplex_20200808 CMVMRC5 NY3 pDC.813537.fcs").unwrap();
        let mut reader = BufReader::new(&file);

        let metadata = read_metadata(&mut reader).expect("Failed to read metadata");
        reader.seek(SeekFrom::Start(150)).expect("Failed to seek to data start");

        let events = read_events::<LittleEndian>(&mut reader, "F", 3, 1, &metadata).expect("Failed to read events");
        assert_eq!(events, vec![2.745084202615544e-6, 11018227712.0, 6.37629560262809e-10]);
    }

    #[test]
    fn test_create_dataframe() {
        let column_titles = vec!["col1".to_string(), "col2".to_string()];
        let data = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let df = create_dataframe(&column_titles, &data).unwrap();
        assert_eq!(df.shape(), (3, 2), "DataFrame shape mismatch");
    }

    #[test]
    fn test_create_dataframe_mismatch() {
        let column_titles = vec!["col1".to_string()];
        let data = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let result = create_dataframe(&column_titles, &data);
        assert!(result.is_err(), "DataFrame creation with mismatched columns and data should fail");
    }
}
