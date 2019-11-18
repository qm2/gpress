
#ifndef COMPRESSOR_H
#define COMPRESSOR_H
#include "hash.h"

//this function used to compress the files
int gtf_compressor(FILE* fp, int length, int* chr_table, int* block_min_table,int* block_max_table, int block_size);

//this function used to decompress the files
int gtf_decompressor(FILE* fp, int length, int filetype);

//this function used to compress the expression matrix files
int expression_compressor(FILE* fp, int length, int block_size);

//this function used to compress the sparse matrix files
int sparse_compressor(FILE* fp, int length, int block_size);

//this function compress the GTF without random access
int gtf_compressor2(FILE* fp, int length, int filetype);

//this function used to compress the files
int gff3_compressor(FILE* fp, int length, int* chr_table, int* block_min_table,int* block_max_table, int block_size);

#endif