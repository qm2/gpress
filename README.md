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
1. To compress the GTF file without random access, run
```
./gpress -cw [inputfile]
```
The compressed file is stored in compressed foler as GTF_compressed_without.tar

2. To decompress the GTF without random access, run 
```
./gpress -dc 
```
The decompressed GTF is stored in the current folder as decompressed_gtf.gtf

3. To compress the GTF with random access, run 
```
./gpress -c [inputfile] <block_size> (e.g. ./gpress -c test_gtf.gtf 500)
```
The compressed GTF file is stored in foler compressed as GTF_compressed.tar.
The associated index tables are also stored in folder compressed.

4. To do queries on compressed GTF file, 
id search:
```
./gpress -q -id <id>
```
The retrieved information is printed in command window.
range search:
```
./gpress -q -range <start> <end> <chromosome>
```
The retrieved information is stored in current folder as range.gtf.

5. To compress and link the expression file, 
```
./gpress -e [inutfile] <block_size>
```
The compressed file is stored in folder compressed as expression_compressed.tar.

6. To do queries on compressed expression file,
```
./gpress -qe <id>
```
The retrieved information is stored in current folder as expression_search.txt. GPress will also print the extra information in command window if it exists in GFF file.

## Input

Here, we provide a small GTF file(test_gtf.gtf) and a small expression file(test_expression.tsv)in the folder data. More sample files can be download from the GENCODE database or other sources.

## Example




