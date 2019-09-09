# GFF compression scheme
This is the GFF compression scheme, a compression platform for GFF3, GTF and expression matrix files.

## Requirement
- GCC compiler

## How to Use the Software:

### Install
This program is used from the command line. After download the GFF compression from https://github.com/qm2/gtf_compressor, you can compile with
```
make
```

### Run
To run the GFF compression scheme, the general command is:
```
./compressor (options)  [inputfile]
```

#### Operating Mode:

-c:           Compress the GFF files without random access

-uc:          Uncompress the GFF files without random access

-r:           Compress the GFF files with support of random access

-q:           Link the compressed GFF files and hashtables and then do searches

-e            Compress the expression matrix files and link relavent information with compressed GFF files 

#### Compression Parameters:

block:        number of blocks used for random access (only available for -r and -e mode)

## Input

 Here, we provide a small test file in the folder data. More sample files can be download from the GENCODE database: https://www.gencodegenes.org/

## Output
The compressed GFF files are called results_*.txt in the results folder

The output decompression GFF files are called decompressed_*.txt in compressed folder

The hastables are compressed and stored in the folder index_tables

The compressed expression matrix files are called matrix_results_*.txt in matrix_results folder

The output decompressed expression matrix files are called matrix_decompressed_*.txt in compressed folder


