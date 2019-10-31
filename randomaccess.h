
#ifndef RANDOMACCESS_H
#define RANDOMACCESS_H
#include "hash.h"

//helper function to read a specific line
int readTaggedLine(char* filename, int line, char* result);

//this function used access a single item by its id
char* item_search(int block, int block_id);

//this function find all the genes between start_pos and end_pos
int rangeSearch(int start_pos, int end_pos, int chr, int* chr_table);

//this function tries to add more indentifiers in other database to a feature
int add_database_id(char* old_id, char* new_id);

//this function finds the item in expression file
int expressionSearch(int block, int start_id, int end_id);



#endif
