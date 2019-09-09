
#ifndef COMPRESSOR_H
#define COMPRESSOR_H
#include "hash.h"

//this function used to compress the files
int gtf_compressor(FILE* fp, int length,  int* chr_table, int block_size);
//this function used to decompress the files
int gtf_decompressor(FILE* fp, int length);

//this function used to compress the expression matrix files
int matrix_compressor(FILE* fp, int length, int block_size);

//this function used to compress the sparse expression matrix files
int sparse_matrix_compressor(FILE* fp, int length, int block_size);

#endif