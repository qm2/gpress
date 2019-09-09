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
./compressor (options)  [inputfile]


####Operating Mode:

-c:           Compress the GFF files without random access

-uc:          Uncompress the GFF files without random access

-r:           Compress the GFF files with support of random access

-q:           Link the compressed GFF files and hashtables and then do searches

-e            Compress the expression matrix files and link relavent information with compressed GFF files 

####Compression Parameters:

block:        number of blocks used for random access (only available for -r and -e mode)

## Input

The testing gtf files are download from Genecode. We used the gencode.v30.annotation.gtf and gencode.vM21.annotation.gtf.

## Output
The compressed file for compression is called compressed.txt.* in the toplevel folder

The output decompression file is called result.gtf in toplevel folder


