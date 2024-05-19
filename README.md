
## Module Overview: `no_bs_fcs`

The `no_bs_fcs` module provides a basic set of tools for manipulating and analyzing Flow Cytometry Standard (FCS) files. It includes methods for reading FCS files, extracting metadata, and processing data segments using arcshinh transformation, all while handling various potential errors gracefully.

Supports FCS versions FCS3.0 and FCS3.1.

### Key Features

- **File Handling**: Open and read FCS files in a structured and efficient manner.
- **Metadata Extraction**: Extract and validate metadata from FCS files, ensuring all required information is available.
- **Data Processing**: Parse data segments from FCS files and convert them into usable formats such as data frames.
- **Error Handling**: Comprehensive error handling to deal with various issues that may arise during file operations.

### Modules

The `no_bs_fcs` module is divided into several sub-modules, each focusing on specific aspects of FCS file handling:

- **data**: Contains structures and functions for handling the data segments of FCS files, including parsing and transformation operations.
- **header**: Includes methods for reading and validating the header segments of FCS files.
- **text**: Provides functions for reading and validating the text segments of FCS files.

### Error Handling

The `FcsError` enum defines various errors that can occur while processing FCS files, including I/O errors, invalid headers, unsupported versions, metadata issues, missing required keywords, and invalid data segments.

### Usage Examples

#### Opening an FCS File

To open an FCS file, use the `FcsFile::open` method:

```rust
use no_bs_fcs::FcsFile;

let fcs_file = FcsFile::open("path/to/file.fcs").unwrap();
```

#### Reading an FCS File

To read an FCS file and extract metadata and data segments:

```rust
use no_bs_fcs::{FcsFile, FcsError};

let fcs_file = FcsFile::open("path/to/file.fcs").unwrap();
let fcs_data = fcs_file.read().unwrap();
println!("{:?}", fcs_data.data); // Prints FCS data in a dataframe
shape: (5_555, 10)
┌───────────┬───────────┬───────────┬───────────┬───┬───────────┬───────────┐
│ FSC-A     ┆ FSC-W     ┆ SSC-SSC-A ┆ Time      ┆ … ┆ FSC-H     ┆ APC-A     │
│ ---       ┆ ---       ┆ ---       ┆ ---       ┆   ┆ ---       ┆ ---       │
│ f64       ┆ f64       ┆ f64       ┆ f64       ┆   ┆ f64       ┆ f64       │
╞═══════════╪═══════════╪═══════════╪═══════════╪═══╪═══════════╪═══════════╡
│ 22.533354 ┆ 6.109114  ┆ 21.705629 ┆ 0.346574  ┆ … ┆ 22.69205  ┆ 7.378071  │
│ 22.114697 ┆ 5.847766  ┆ 6.454475  ┆ 22.523151 ┆ … ┆ 21.319135 ┆ 23.476587 │
│ 5.972638  ┆ 15.808522 ┆ 14.130202 ┆ 22.09955  ┆ … ┆ 5.869126  ┆ 23.012064 │
│ 13.866457 ┆ 23.927268 ┆ 15.367728 ┆ NaN       ┆ … ┆ 5.847766  ┆ 6.438551  │
│ …         ┆ …         ┆ …         ┆ …         ┆ … ┆ …         ┆ …         │
│ 6.762567  ┆ 11.167    ┆ 6.762567  ┆ 22.140327 ┆ … ┆ 2.988391  ┆ 21.279063 │
│ 6.592359  ┆ 22.740847 ┆ 13.824292 ┆ 5.714589  ┆ … ┆ 18.447066 ┆ 3.957715  │
│ 13.03357  ┆ 21.456107 ┆ 23.643495 ┆ 5.759973  ┆ … ┆ 10.021293 ┆ 15.392062 │
│ 22.713507 ┆ 5.847766  ┆ 23.071459 ┆ 11.738598 ┆ … ┆ 23.482336 ┆ 15.443761 │
└───────────┴───────────┴───────────┴───────────┴───┴───────────┴───────────┘
```

View information relating to the sample.
```rust
use no_bs_fcs::{FcsFile, FcsError};

let fcs_file = FcsFile::open("path/to/file.fcs").unwrap();
let fcs_data = fcs_file.read().unwrap();
println!("{}", fcs_data.data); // Prints sample information
FlowSample:
    Machine: 1234567 Attune NxT Acoustic Focusing Cytometer (Lasers: BRVY)
    Begin Time:
    End Time:
    Date:
    File: file.fcs
    Volume run: 250000
    Labels:
        Time (Time)
        FSC-A (FSC-A)
        SSC-A (SSC-A)
        BL1-A (CMV-GFP-FITC-A)
        BL3-A (CD11c-PerCP-Cy5.5-A)
        YL4-A (--PE-Cy7-A)
        YL2-A (TNF alpha ICS-PE-Texas Red-A)
        YL1-A (CD123-R-PE-A)
        RL1-A (IFN alpha ICS-APC-A)
        RL2-A (CD45RA-Alexa Fluor 700-A)
        RL3-A (CD45-APC-Cy7-A)
        VL4-A (--Qdot 700-A)
        VL3-A (--Qdot 605-A)
        VL2-A (Live-Dead-Fixable Aqua-A)
        VL1-A (CD69-Pacific Blue-A)
        FSC-H (FSC-H)
        SSC-H (SSC-H)
        FSC-W (FSC-W)
        SSC-W (SSC-W)
```

#### Perform Arcsinh Transformation on FCS data

To read an FCS file and extract metadata and data segments:

```rust
use no_bs_fcs::{FcsFile, FcsError};

let fcs_file = FcsFile::open("path/to/file.fcs").unwrap();
let fcs_data = fcs_file.read().unwrap();

// Get the column names
let column_names = flow_sample.get_dataframe_columns();
    
// Perform arcsinh transformation of data
flow_sample.arcsinh_transform(5.0, &column_names).unwrap();
```

#### Creating a DataFrame

To create a DataFrame from column titles and corresponding data vectors:

```rust
use no_bs_fcs::data::create_dataframe;
use polars::prelude::*;

let column_titles = vec!["APC-A".to_string(), "FSC-W".to_string()];
let data = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
let df = create_dataframe(&column_titles, &data).unwrap();
println!("{:?}", df);
```

### Conclusion

The `no_bs_fcs` module provides a robust framework for handling FCS files. Users can easily and efficiently work with flow cytometry data.

---

This overview provides a comprehensive look at your module, highlighting its features, usage, and testing. Feel free to adjust the paths, example data, and other details to better match your actual implementation and use cases.