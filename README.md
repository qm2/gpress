# GPress
This is the GPress, a framework for querying GFF/GTF and expression files in a compressed form.

## Requirement
- Linux
- GCC and G++ 11 compiler
- At least 2GB of RAM and 20GB of storage

## How to Use the Software:

### Install
This program is used from the command line. After download GPress from https://github.com/qm2/gtf_compressor, you can compile GPress with
```
make
```
You also need to compile the BSC compressor by runing 
```
make
```
in the BSC folder.
### Run
To run the GPress, the general command is:
```
./gpress (options)  [inputfile]
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


