# GPress
This is GPress, a framework for querying GFF/GTF and expression files in a compressed form.

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
in the root folder.

You also need to compile the BSC compressor by runing 
```
make
```
in the **BSC** folder.
### Run
To run the GPress, the general command is:
```
./gpress (options)  [inputfile]
```
in the root folder
#### Operating Mode:
1. To compress the GTF file without random access, run
```
./gpress -cw [inputfile]
```

2. To decompress the GTF without random access, run 
```
./gpress -dc 
```

3. To compress the GTF with random access, run 
```
./gpress -c [inputfile] <block_size> 
```

4. To do queries on compressed GTF file, 

id search:
```
./gpress -q -id <id>
```
range search:
```
./gpress -q -range <start> <end> <chromosome>
```

5. To compress and link the expression file, 
```
./gpress -e [inutfile] <block_size>
```

6. To do queries on compressed expression file,
```
./gpress -qe <id>
```

## Input

Here, we provide a small GTF file (test_gtf.gtf) and a small expression file (test_expression.tsv)in the folder **data**. More sample files can be download from the GENCODE database or other sources.

## Example
Here, we will use the test files to provide an example of how to use GPress.
1. To compress the GTF file without random access, run
```
./gpress -cw data/test_gtf.gtf
```
The compressed file is stored named **GTF_compressed_without.tar** in folder **compressed**.

2. To decompress the GTF without random access, run 
```
./gpress -dc 
```
The decompressed GTF is stored named **decompressed_gtf.gtf** in the root folder.

3. To compress the GTF with random access (500 genes per block), run 
```
./gpress -c data/test_gtf.gtf 500
```
The compressed GTF file is stored named **GTF_compressed.tar** in folder **compressed**.
The associated index tables are also stored in folder **compressed**.

4. To do queries on compressed GTF file, 

id search for exon id "ENSE00003486434.1":
```
./gpress -q -id ENSE00003486434.1
```
The retrieved information is printed in command window.

range search on chromosome 1 with range from 10000 to 100000:
```
./gpress -q -range 10000 100000 1
```
The retrieved information is stored named **range.gtf** in root folder.

5. To compress and link the expression file, 
```
./gpress -e data/test_expression.tsv 500
```
The compressed file is stored named **expression_compressed.tar** in folder **compressed**.

6. To do queries on compressed expression file,
```
./gpress -qe ENST00000009530.11
```
The retrieved information is stored named **expression_search.txt** in root folder.

GPress will also print the extra information in command window if it exists in GFF file. For example, 
```
./gpress -qe ENST00000531822.1
```
will give information from both GTF and expression files.




