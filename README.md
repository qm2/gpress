# GFF compression scheme
This is the baseline compressor for gtf files. It simply compresses the file without support for random access.

# Usage:
This program is used from the command line. After building with make, the general format is:

./compressor (options)  [inputfile]


Operating Mode:

-c:           Compress the gtf file

-uc:          Uncompress the files 

Compression Parameters:

-gzip:        Use gzip to compress all columns

-bzip2:       Use bzip2 to compress all columns

-xz:          Use xz to compress all columns

#Input

The testing gtf files are download from Genecode. We used the gencode.v30.annotation.gtf and gencode.vM21.annotation.gtf.

# Output
The compressed file for compression is called compressed.txt.* in the toplevel folder

The output decompression file is called result.gtf in toplevel folder


