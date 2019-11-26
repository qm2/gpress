# GPress
This is GPress, a framework for querying GTF, GFF3 and expression files in a compressed form.

## Requirement
- Linux
- GCC and G++ 11 compiler
- At least 2GB of RAM and 20GB of storage

## How to Use the Software:

### Install
This program is used from the command line. First, download GPress using git command:
```
git clone https://github.com/qm2/gpress gpress
```
Then, enter into the root folder:
```
cd gpress
```
you can then compile GPress with:
```
make
```
in the root folder.

You also need to compile the BSC compressor

First, enter into the BSC folder:
```
cd BSC
```
Then, compile the BSC compressor:
```
make
```
in the **BSC** folder.

Finally, go back to the root folder
```
cd ..
```

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

4. To do single id search on compressed GTF or GFF3 file from the specified folder, run 
```
./gpress -q -id [id] [folder]
```
By default, the outputs are printed on the command window. 

In order to write the output to a file, run

```
./gpress -q -id [id] [folder] > [outputfile]
```

5. To do multiple id searches on compressed GTF or GFF3 file from the specified folder, run 
```
./gpress -q -id [folder]
```
Then, the user can enter the ID to search the item and its related parents and children.

6. To do range search on compressed GTF or GFF3 file from the specified folder, run 
```
./gpress -q -range [start] [end] [chromosome] [folder]
```

By default, the outputs are printed on the command window. 

In order to write the output to a file, run
```
./gpress -q -range [start] [end] [chromosome] [folder] > [outputfile]
```

7. To compress and link the expression file, 
```
./gpress -e [inputfile] <block_size> [folder]
```
The compressed expression file is linked to the GFF files in specified folder and are also store in that folder

block size is 2000 by default

8. To do id search on compressed expression file from a specified folder,
```
./gpress -qe [id] [folder]
```
By default, the outputs are printed on the command window. 

In order to write the output to a file, 

```
./gpress -qe [id] [folder] > [outputfile]
```

9. To do range search on compressed expression file based on compressed GFF file from a specified folder ,
```
./gpress -qer [start] [end] [chromosome] [folder]
```
By default, the outputs are printed on the command window. 

In order to write the output to a file, 

```
./gpress -qer [start] [end] [chromosome] [folder] > [outputfile]
```


10. To compress and link the sparse matrix file, 
```
./gpress -sparse [input matrix] [input genes] [input barcodes] <block_size> [folder]
```
The compressed expression files are store in the specified folder

block size is 2000 by default

11. To do id search on compressed sparse expression matrix file from a specified folder,
```
./gpress -qs [id] [folder]
```
By default, the outputs are printed on the command window. 

In order to write the output to a file, 

```
./gpress -qs [id] [folder] > [outputfile]
```

12. To do range search on compressed sparse matrix file based on a compressed GFF file from a specified folder ,
```
./gpress -qsr [start] [end] [chromosome] [folder]
```
By default, the outputs are printed on the command window. 

In order to write the output to a file, 

```
./gpress -qsr [start] [end] [chromosome] [folder] > [outputfile]
```


## Test files

The folder “data” contains a small GTF file (**test_gtf.gtf**), a small GFF3 file (**test_gff3.gff3**), a small expression file resulting from bulk RNA-Seq (**test_expression.tsv**), and the files resulting from a single-cell RNA-Seq experiment (**matrix.mtx, genes.tsv, barcodes.tsv**). More sample files can be downloaded from the GENCODE or other databases.

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

4. To do single id search on compressed GTF file from folder **gtf1**, run
```
./gpress -q -id ENSE00003486434.1 gtf1
```
The retrieved information will be printed in command window by default.

5. To do multiple id searches on compressed GTF file from folder **gtf1**, run
```
./gpress -q -id gtf1
```
Then enter the exon id "ENSE00003486434.1" in the command window, the retrieved information will be printed in command window.
If the user also wants to know all this exon's parents and children contained in the gene, enter "yes" in the command window. Finally, hit "quit" to terminate the program.

6. To do range search on chromosome 1 with range from 10000 to 100000:
```
./gpress -q -range 10000 100000 1 gtf1 
```
The retrieved information is printed in command window.

7. To compress and link the expression file with GFF files in folder **gtf1** (500 genes per block), 
```
./gpress -e data/test_expression.tsv 500 gtf1
```
The compressed file is stored named **expression_compressed.tar** in folder **gtf1**.

8. To do id search on compressed expression file from folder **gtf1**,
```
./gpress -qe ENST00000009530.11 gtf1 
```
The retrieved information is printed in command window.

While parsing the expression file, GPress store all the IDs of the same item from multiple databases. For example, 
```
./gpress -qe OTTHUMT00000374178.1 gtf1 
```
will retrieve same item as above.

GPress will also print the extra information if the searched item exists in the compressed GFF file. For example, 
```
./gpress -qe ENST00000531822.1 gtf1 
```
will give information from both GTF and expression files in command window.

9. To do range search on chromosome 1 with range from 100000 to 100000000 on compressed expression file based on compressed GFF file from folder **gtf1**,
```
./gpress -qer 100000 100000000 1 gtf1
```
The items on chromosome 1 with range from 100000 to 100000000 that exist in both GFF files and expression file are retrieved.

By default, the outputs are printed on the command window.

10. To compress and link the sparse expression matrix file in folder **gtf1** (500 genes per block), 
```
./gpress -sparse data/matrix.mtx data/genes.tsv data/barcodes.tsv 500 gtf1
```
The compressed sparse matrix files are stored named **sparse_compressed.tar** in folder **gtf1**.

11. To do id search on compressed sparse sparse matrix file from folder **gtf1**,
```
./gpress -qs ENSG00000242485 gtf1
```
GPress will also print the extra information if it exists in GFF file.

By default, the outputs are printed in the command window. 

12. To do range search on chromosome 1 with range from 1000 to 1000000 on compressed sparse matrix file based on compressed GFF file from folder **gtf1**,
```
./gpress -qsr 1000 1000000 1 gtf1
```
The items on chromosome 1 with range from 1000 to 1000000 that exist in both GFF files and sparse matrix files are retrieved.

By default, the outputs are printed in the command window.





