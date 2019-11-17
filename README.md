# GPress
This is GPress, a framework for querying GTF, GFF3 and expression files in a compressed form.

## Requirement
- Linux
- GCC and G++ 11 compiler
- At least 2GB of RAM and 20GB of storage

## How to Use the Software:

### Install
This program is used from the command line. After download GPress from https://github.com/qm2/gpress, you can compile GPress with
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
./gpress (options)  [inputfile] [parameters] [folder]
```
in the root folder.

#### Operating Mode:
1. To compress the GTF or GFF3 file without random access, run
```
./gpress -cw [inputfile] [folder]
```
The compressed files are stored in the new specified folder

2. To decompress the GTF or GFF3 file from the specified folder without random access, run 
```
./gpress -dc [filetype] [folder]
```
filetype should be either "gtf" or "gff3".

3. To compress the GTF or GFF3 with random access, run 
```
./gpress -c [inputfile] <block_size> [folder]
```
The compressed files are stored in the new specified folder

block size is 2000 by default

4. To do queries on compressed GTF or GFF3 file from the specified folder, 

id search:
```
./gpress -q -id [id] [folder]
```
range search:
```
./gpress -q -range [start] [end] [chromosome] [folder]
```

By default, the outputs are printed on the command window. 

In order to write the output to a file, 

id search:
```
./gpress -q -id [id] [folder] > [outputfile]
```
range search:
```
./gpress -q -range [start] [end] [chromosome] [folder] > [outputfile]
```

5. To compress and link the expression file, 
```
./gpress -e [inputfile] <block_size> [folder]
```
The compressed expression files are linked to the GFF files in specified folder and are also store in that folder

block size is 2000 by default

6. To do queries on compressed expression file from a specified folder,
```
./gpress -qe [id] [folder]
```
By default, the outputs are printed on the command window. 

In order to write the output to a file, 

```
./gpress -qe [id] [folder] > [outputfile]
```

## Input

Here, we provide a small GTF file (**test_gtf.gtf**), a small GFF3 file (**test_gff3.gff3**) and a small expression file (**test_expression.tsv**)in the folder **data**. More sample files can be download from the GENCODE database or other databases.

## Example
Here, we will use the test files to provide an example of how to use GPress. We will use GTF file for illustration since the GFF3 file is nearly the same except slight format variations. We will store the compressed GTF file and compressed expression file in a new folder **gtf1**
1. To compress the GTF file without random access, run
```
./gpress -cw data/test_gtf.gtf gtf1
```
The compressed file is stored named **GTF_compressed_without.tar** in folder **gtf1**.

2. To decompress the GTF file from folder **gtf1** without random access, run 
```
./gpress -dc gtf gtf1
```
The decompressed GTF is stored named **decompressed_gtf.gtf** in folder **output**.

3. To compress the GTF with random access (500 genes per block), run 
```
./gpress -c data/test_gtf.gtf 500 gtf1
```
The compressed GTF file is stored named **GTF_compressed.tar** in folder **gtf1**.
The associated index tables are also stored in folder **gtf1**.

4. To do queries on compressed GTF file from folder **gtf1**, 

id search for exon id "ENSE00003486434.1":
```
./gpress -q -id ENSE00003486434.1 gtf1
```
The retrieved information is printed in command window.

range search on chromosome 1 with range from 10000 to 100000:
```
./gpress -q -range 10000 100000 1 gtf1 > output/range_search.gtf
```
The retrieved information is stored named **range_search.gtf** in folder **output**.

5. To compress and link the expression file with GFF files in folder **gtf1** (500 genes per block), 
```
./gpress -e data/test_expression.tsv 500 gtf1
```
The compressed file is stored named **expression_compressed.tar** in folder **gtf1**.

6. To do queries on compressed expression file from folder **gtf1**,
```
./gpress -qe ENST00000009530.11 gtf1 > output/expression_search.txt
```
The retrieved information is stored named **expression_search.txt** in folder **output**.

GPress will also print the extra information if it exists in GFF file. For example, 
```
./gpress -qe ENST00000531822.1 gtf1 > output/expression_search2.txt
```
will give information from both GTF and expression files stored as **expression_search2.txt** in folder **output**.




