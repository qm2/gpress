
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compressor.h"
#include "hash.h"
#include "randomaccess.h"
#include "linked_list.h"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc , char **argv){
    int *chr_table= (int*)malloc(sizeof(int)*100);
    if ( (argc <= 1 ) || (strcmp("-h", argv[1])==0)){
    	printf("Welcome to GPress: a framework for querying GTF, GFF3 and expression files in a compressed form.\n");
    	// printf("To run the GPress, the general command is:\n");
     //    printf("./gpress (options)  [inputfile] [parameters] [folder] in the root folder.\n");
        int c;
		FILE *file;
		file = fopen("help_message.txt", "r");
		if (file) {
		    while ((c = getc(file)) != EOF){
		        putchar(c);		    	
		    }
		    fclose(file);
		}
		printf("\n");
        return 0;
    }
    if (strcmp("-c", argv[1]) == 0){

        //file pointer for gtf file
        FILE *fp, *fp_chromosome, *fp_min, *fp_max;
        //create the position min and max tables
        int *min_table= (int*)malloc(sizeof(int)*500);
        int *max_table= (int*)malloc(sizeof(int)*500);
        int count_lines = 0;
        char chr;
        fp = fopen(argv[3], "r");
        if(fp == NULL){
            printf("the input file is invalid!\n");
            return 0;
        }
        //count number of lines in the file
        chr = getc(fp);
        while (chr != EOF)
        {
            //Count whenever new line is encountered
            if (chr == '\n')
            {
                count_lines = count_lines + 1;
            }
            //take next character from file.
            chr = getc(fp);
        }
        fclose(fp);
        fp = fopen(argv[3], "r");
        char *dot = strrchr(argv[3], '.');
        int block;
        if(argc == 4){
        	block = 2000;
        }
        else{
        	block = atoi(argv[2]);
        }
        if(!strcmp(dot+1, "gtf")){
            count_lines -=5;
            gtf_compressor(fp, count_lines, chr_table, min_table, max_table, block);
        }
        else if(!strcmp(dot+1, "gff3")){
            count_lines -=7;
            gff3_compressor(fp, count_lines, chr_table, min_table, max_table, block);
        }
        else{
            printf("The input name is invalid!\n");
            return 0;
        }
        count_lines -=5;
        //store the index table for chromosome
        fp_chromosome= fopen("index_tables/data_chr.txt", "w+");
        fwrite(chr_table , sizeof(int) , sizeof(chr_table) , fp_chromosome);
        fclose(fp_chromosome);
        //store the block min position table
        fp_min= fopen("index_tables/data_min.txt", "w+");
        fwrite(min_table , sizeof(int) , sizeof(min_table) , fp_min);
        fclose(fp_min);
        //store the block max position table
        fp_max= fopen("index_tables/data_max.txt", "w+");
        fwrite(max_table , sizeof(int) , sizeof(max_table) , fp_max);
        fclose(fp_max);
        //create the folder specified by the user
        char command7[200];
        snprintf(command7, sizeof(command7), "mkdir %s", argv[argc-1]);
        struct stat st = {0};
        if (stat(argv[argc-1], &st) == -1) {
            system(command7);
        }
        char command1[200];
        char command2[200];
        char command3[200];
        char command4[200];
        char command5[200];
        char command6[200];
        system("rm GTF_parsed/*");
        snprintf(command1, sizeof(command1), "tar -cf %s/GTF_compressed.tar GTF_compressed", argv[argc-1]);
        system(command1);
        system("rm GTF_compressed/*");
        snprintf(command2, sizeof(command2), "BSC/bsc e index_tables/data_key.txt %s/data_key_compressed", argv[argc-1]);
        system(command2);
        snprintf(command3, sizeof(command3), "BSC/bsc e index_tables/data_value.txt %s/data_value_compressed", argv[argc-1]);
        system(command3);
        snprintf(command4, sizeof(command4), "BSC/bsc e index_tables/data_chr.txt %s/data_chr_compressed", argv[argc-1]);
        system(command4);
        snprintf(command5, sizeof(command5), "BSC/bsc e index_tables/data_min.txt %s/data_min_compressed", argv[argc-1]);
        system(command5);
        snprintf(command6, sizeof(command6), "BSC/bsc e index_tables/data_max.txt %s/data_max_compressed", argv[argc-1]);
        system(command6);
        system("rm index_tables/*");
        printf("The compression of GFF file with random access succeeds! The compressed files are in folder %s\n", argv[argc-1]);
        fclose(fp);
    }
    else if (strcmp("-cw", argv[1]) == 0){

        //file pointer for gtf file
        FILE *fp;
        int count_lines = 0;
        char chr;
        fp = fopen(argv[2], "r");
        if(fp == NULL){
            printf("the input file is invalid!\n");
            return 0;
        }
        //count number of lines in the file
        chr = getc(fp);
        while (chr != EOF)
        {
            //Count whenever new line is encountered
            if (chr == '\n')
            {
                count_lines = count_lines + 1;
            }
            //take next character from file.
            chr = getc(fp);
        }
        fclose(fp);
        fp = fopen(argv[2], "r");
        char *dot = strrchr(argv[2], '.');
        if(!strcmp(dot+1, "gtf")){
            count_lines -=5;
            gtf_compressor2(fp, count_lines,0);
        }
        else if(!strcmp(dot+1, "gff3")){
            count_lines -=7;
            gtf_compressor2(fp, count_lines, 1);
        }
        else{
            printf("The input name is invalid!\n");
            return 0;
        }
        //create the folder specified by the user
        char command2[200];
        snprintf(command2, sizeof(command2), "mkdir %s", argv[argc-1]);
        struct stat st = {0};
        if (stat(argv[argc-1], &st) == -1) {
            system(command2);
        }
        char command1[200];
        snprintf(command1, sizeof(command1), "tar -cf %s/GTF_compressed_without.tar GTF_compressed2", argv[argc-1]);
        system(command1);
        system("rm GTF_parsed2/*");
        system("rm GTF_compressed2/*");
        printf("The compression of GTF file without random access succeeds!\n");
        fclose(fp);
    }
    else if (strcmp("-dc", argv[1]) == 0){

        //file pointer for gtf file
        int count_lines = 0;
        char chr;
        FILE *fp, *fp2;

        //run the decompressor
        char command1[200];
        snprintf(command1, sizeof(command1), "tar -xf %s/GTF_compressed_without.tar GTF_compressed2", argv[argc-1]);
        system(command1);
        system("BSC/bsc d GTF_compressed2/gtf_seqname_compressed GTF_parsed2/gtf_seqname.txt");
        system("BSC/bsc d GTF_compressed2/gtf_source_compressed GTF_parsed2/gtf_source.txt");
        system("BSC/bsc d GTF_compressed2/gtf_feature_compressed GTF_parsed2/gtf_feature.txt");
        system("BSC/bsc d GTF_compressed2/gtf_start_compressed GTF_parsed2/gtf_start.txt");
        system("BSC/bsc d GTF_compressed2/gtf_delta_compressed GTF_parsed2/gtf_delta.txt");
        system("BSC/bsc d GTF_compressed2/gtf_score_compressed GTF_parsed2/gtf_score.txt");
        system("BSC/bsc d GTF_compressed2/gtf_strand_compressed GTF_parsed2/gtf_strand.txt");
        system("BSC/bsc d GTF_compressed2/gtf_attribute_compressed GTF_parsed2/gtf_attribute.txt");
        system("BSC/bsc d GTF_compressed2/gtf_gtf_frame_cds_compressed GTF_parsed2/gtf_frame_cds.txt");
        system("BSC/bsc d GTF_compressed2/gtf_frame_start_compressed GTF_parsed2/gtf_frame_start.txt");
        system("BSC/bsc d GTF_compressed2/gtf_frame_stop_compressed GTF_parsed2/gtf_frame_stop.txt");
        fp2 = fopen("GTF_parsed2/gtf_seqname.txt", "r");
        if(fp2 == NULL){
            printf("the compressed files are invalid!\n");
            return 0;
        }
        //count number of lines in the file
        chr = getc(fp2);
        while (chr != EOF)
        {
            //Count whenever new line is encountered
            if (chr == '\n')
            {
                count_lines = count_lines + 1;
            }
            //take next character from file.
            chr = getc(fp2);
        }
        fclose(fp2);


        if(!strcmp(argv[2], "gtf")){
            fp = fopen("output/decompressed_gtf.gtf", "w+");
            gtf_decompressor(fp, count_lines, 0);
            printf("The decompressed GTF file is included in the output/decompressed_gtf.gtf!\n");
        }
        else if(!strcmp(argv[2], "gff3")){
            fp = fopen("output/decompressed_gff3.gff3", "w+");
            gtf_decompressor(fp, count_lines, 1);
            printf("The decompressed GFF3 file is included in the output/decompressed_gff3.gff3!\n");
        }
        else{
            printf("The file type is invalid!\n");
            return 0;
        }
        system("rm GTF_parsed2/*");
        system("rm GTF_compressed2/*");
        fclose(fp);
    }
    else if(strcmp("-q", argv[1]) == 0){

        char hash_key[500];
        char* hash_val;
        char temp[100];

        //recover all the data structures

        char command1[200];
        char command2[200];
        char command3[200];
        snprintf(command1, sizeof(command1), "BSC/bsc d %s/data_key_compressed index_tables/data_key.txt", argv[argc-1]);
        system(command1);
        snprintf(command2, sizeof(command2), "BSC/bsc d %s/data_value_compressed index_tables/data_value.txt", argv[argc-1]);
        system(command2);
        snprintf(command3, sizeof(command3), "BSC/bsc d %s/data_chr_compressed index_tables/data_chr.txt", argv[argc-1]);
        system(command3);

        char command4[200];
        snprintf(command4, sizeof(command4), "tar -xf %s/GTF_compressed.tar GTF_compressed", argv[argc-1]);
        system(command4);

        if(strcmp("-id", argv[2]) == 0){
  	        FILE *fp_hash_key= fopen("index_tables/data_key.txt", "r");
  	        FILE *fp_hash_val= fopen("index_tables/data_value.txt", "r");
            if(argc==4){
      	        hashtable_t *ht = ht_create(3000000);
      	        while(fscanf(fp_hash_key, "%s", hash_key)!=EOF){
      	            hash_val= (char*)malloc(sizeof(char)*50);
      	            fgets(hash_val, 50, fp_hash_val);
      	            hash_val[strlen(hash_val) - 1] = '\0';
      	            ht_put(ht, hash_key, hash_val);
  	            }
	  	        //close the files
	  	        fclose(fp_hash_key);
	  	        fclose(fp_hash_val);
	            char* retval;
	            char* hashval;
	            char hash_split[200];
	            char* s;
	            int block;
	            int block_id;
	            int parsed =0;
	            printf("Welcome to the random access based on ID\n");
	            while(1){
	                printf("Please enter a valid ID for search or enter 'quit' to terminate the program!\n");
	                char input[500];
	                scanf("%[^\n]%*c", input);
	                if(strcmp(input, "quit") == 0){
	                	break;
	                }
	                //get rid of the version after '.'
	                int count;
	                for(count=0; count<strlen(input); count++){
	                    if(input[count] == '.'){
	                        input[count]='\0';
	                        break;
	                    }
	                }
	                hashval= (char*)ht_get(ht, input);
	                if(hashval == NULL){
	                   printf("This ID is not valid!\n");
	                   continue;
	                }
	                strcpy(hash_split, hashval);
	                s= strtok(hash_split, " ");
	                block= atoi(s);
	                s= strtok(NULL, " ");
	                block_id= atoi(s);
	                retval = item_search(block, block_id);
	                printf("The searched item is:\n");
	                printf("%s", retval);
	                printf("Please enter 'yes' if you also want to print out all its parents and children within a gene\n");
	                printf("Otherwise, enter 'no' to skip it\n");
	                while(1){
	                    scanf("%[^\n]%*c", input);
	                    if(strcmp(input, "yes") == 0){
	                	    FILE* fp_gene = fopen("gene_block.gtf", "r");
	                	    char info[1500];
	                	    while(fgets(info, 1024, fp_gene)){
	                            printf("%s", info);
	                	    }
	                	    system("rm info.gtf");
	                	    system("rm gene_block.gtf");
	                        break;
	                    }
	                    else if(strcmp(input, "no") == 0){
	                        break;
	                    }
	                    else{
	                    	printf("the input value is not valid, please enter again!\n");
	                    }
	                }
	                printf("id search succeeds!\n");
	                parsed = 1;
	            }
	            system("rm GTF_compressed/*");
	            if(parsed==1){
	                system("rm GTF_parsed/*");
	            }
          }
          else if(argc == 5 || argc == 6){
              FILE *fp_hash_key= fopen("index_tables/data_key.txt", "r");
    	      FILE *fp_hash_val= fopen("index_tables/data_value.txt", "r");
              char* retval;
              char* no_dot;
              char hash_split[200];
              char* s;
              int block;
              int block_id;
              int exist = 0;
              int family = 0;
              //get rid of the dot in id
              if(!strcmp(argv[3], "-f")){
                  family = 1;
              }
              
              no_dot= strtok(argv[3+family], ".");
              while(fscanf(fp_hash_key, "%s", hash_key)!=EOF){
      	          hash_val= (char*)malloc(sizeof(char)*50);
                  fgets(hash_val, 50, fp_hash_val);
                  if(!strcmp(hash_key, no_dot)){
                      exist = 1;
                      break;
                  }
              }
              if(exist == 0){
                 fprintf(stderr, "This ID is not valid!\n");
                 return 0;
              }
              strcpy(hash_split, hash_val);
              s= strtok(hash_split, " ");
              block= atoi(s);
              s= strtok(NULL, " ");
              block_id= atoi(s);
              retval = item_search(block, block_id);
              printf("%s", retval);
              if(family == 1){
              	    printf("The parents and children of this item are:\n");
                    FILE* fp_gene = fopen("gene_block.gtf", "r");
	                char info[1500];
            	    while(fgets(info, 1024, fp_gene)){
                        printf("%s", info);
            	    }
            	    fclose(fp_gene);            	
              }
              fprintf(stderr, "id search succeeds!\n");
              fclose(fp_hash_key);
              fclose(fp_hash_val);
              system("rm info.gtf");
              system("rm gene_block.gtf"); 
              system("rm GTF_compressed/*");
              if(exist ==1){
                  system("rm GTF_parsed/*");
              }
          }
        }
        else if(strcmp("-range", argv[2]) == 0){
        	//decompress the index tables for min and max positions
        	char command5[200];
	        char command6[200];
	        snprintf(command5, sizeof(command5), "BSC/bsc d %s/data_min_compressed index_tables/data_min.txt", argv[argc-1]);
	        system(command5);
	        snprintf(command6, sizeof(command6), "BSC/bsc d %s/data_max_compressed index_tables/data_max.txt", argv[argc-1]);
	        system(command6);
	        //create the position min and max tables
	        int *min_table= (int*)malloc(sizeof(int)*500);
	        int *max_table= (int*)malloc(sizeof(int)*500);
        	//recover the chromosome table
            FILE *fp_chromosome = fopen("index_tables/data_chr.txt", "rb");
            fread(chr_table, sizeof(int), sizeof(chr_table), fp_chromosome);
        	//recover the block min position table
            FILE *fp_min = fopen("index_tables/data_min.txt", "rb");
            fread(min_table, sizeof(int), sizeof(min_table), fp_min);
        	//recover the block max position table
            FILE *fp_max = fopen("index_tables/data_max.txt", "rb");
            fread(max_table, sizeof(int), sizeof(max_table), fp_max);
            //close the files
            fclose(fp_chromosome);
            fclose(fp_min);
            fclose(fp_max);
            rangeSearch(atoi(argv[3]), atoi(argv[4]), atoi(argv[5])-1, chr_table, min_table, max_table);
            //free the tables
            free(chr_table);
            free(min_table);
            free(max_table);
            fprintf(stderr, "range search succeeds!\n");
            system("rm GTF_compressed/*");
            system("rm GTF_parsed/*");
        }


    }

    else if(strcmp("-e", argv[1]) == 0){
         //file pointer for gtf file
        FILE *fp;
        int count_lines = 0;
        char chr;
        fp = fopen(argv[3], "r");
        if(fp == NULL){
            fprintf(stderr, "the input file is invalid!\n");
            return 0;
        }
        //count number of lines in the file
        chr = getc(fp);
        while (chr != EOF)
        {
            //Count whenever new line is encountered
            if (chr == '\n')
            {
                count_lines = count_lines + 1;
            }
            //take next character from file.
            chr = getc(fp);
        }
        fclose(fp);
        count_lines -=1;
        fp = fopen(argv[3], "r");
        //run the compressor
        int block;
        if(argc == 4){
        	block = 2000;
        }
        else{
        	block = atoi(argv[2]);
        }
        expression_compressor(fp, count_lines, block);
        system("rm expression_parsed/*");
       //create the folder specified by the user
        char command2[200];
        snprintf(command2, sizeof(command2), "mkdir %s", argv[argc-1]);
        struct stat st = {0};
        if (stat(argv[argc-1], &st) == -1) {
            system(command2);
        }
        char command1[200];
        snprintf(command1, sizeof(command1), "tar -cf %s/expression_compressed.tar expression_compressed", argv[argc-1]);
        system(command1);
        system("rm expression_compressed/*");
        printf("compression and linking of expression file succeeds!\n");
        fclose(fp);
    }
    else if(strcmp("-qe", argv[1]) == 0){
        char hash_key[500];
        char* hash_val;
        char temp[100];
        char* retval;
        char command1[200];
        snprintf(command1, sizeof(command1), "BSC/bsc d %s/data_key_compressed index_tables/data_key.txt", argv[argc-1]);
        system(command1);
        char command2[200];
        snprintf(command2, sizeof(command2), "BSC/bsc d %s/data_value_compressed index_tables/data_value.txt", argv[argc-1]);
        system(command2);

        FILE *fp_hash_key= fopen("index_tables/expression_key.txt", "r");
        FILE *fp_hash_val= fopen("index_tables/expression_value.txt", "r");
        hashtable_t *ht = ht_create(3000000);

        while(fscanf(fp_hash_key, "%s", hash_key)!=EOF){
            hash_val= (char*)malloc(sizeof(char)*50);
            fgets(hash_val, 50, fp_hash_val);
            hash_val[strlen(hash_val) - 1] = '\0';
            ht_put(ht, hash_key, hash_val);
        }
        //close the files
        fclose(fp_hash_key);
        fclose(fp_hash_val);
        char command3[200];
        snprintf(command3, sizeof(command3), "tar -xf %s/expression_compressed.tar expression_compressed", argv[argc-1]);
        system(command3);
        char* hashval;
        char* s;
        int block;
        int block_id;
        int block_start_id;
        int block_end_id;
        //get rid of the version after '.'
        int count;
        for(count=0; count<strlen(argv[2]); count++){
            if(argv[2][count] == '.'){
                argv[2][count]='\0';
                break;
            }
        }
        hashval= (char*)malloc(sizeof(char)*100);
        hashval= (char*)ht_get(ht, argv[2]);
        if(hashval == NULL){
            fprintf(stderr, "This ID is not valid!\n");
            return 0;
        }
        s= (char*)malloc(sizeof(char)*50);
        s= strtok(hashval, " ");
        block= atoi(s);
        s= strtok(NULL, " ");
        block_start_id= atoi(s);
        s= strtok(NULL, " ");
        block_end_id= atoi(s);
        expressionSearch(block, block_start_id, block_end_id);
        //check if GTF file contains extra information
        fp_hash_key= fopen("index_tables/data_key.txt", "r");
        fp_hash_val= fopen("index_tables/data_value.txt", "r");
        int exist = 0;
        while(fscanf(fp_hash_key, "%s", hash_key)!=EOF){
            fgets(hash_val, 50, fp_hash_val);
            if(!strcmp(hash_key, argv[2])){
                exist = 1;
                break;
            }
        }
        if(exist == 1){
            char command4[200];
            snprintf(command4, sizeof(command4), "tar -xf %s/GTF_compressed.tar GTF_compressed", argv[argc-1]);
            system(command4);
            s= strtok(hash_val, " ");
            block= atoi(s);
            s= strtok(NULL, " ");
            block_id= atoi(s);
            retval = item_search(block, block_id);
            printf("The item with this id also exists in GFF file:\n");
            printf("%s", retval);
            system("rm GTF_compressed/*");
            system("rm GTF_parsed/*");
        }
        system("rm expression_compressed/*");
        system("rm expression_parsed/*");
        fprintf(stderr, "expression search succeeds!\n");
    }
    else if(strcmp("-sparse", argv[1]) == 0){
        //sort the sparse matrix file
        FILE *fp, *fp_gene;
        FILE *fp2;
        char* s;
        char comments[1000];
        char line[1000];
        fp = fopen(argv[3], "r");
        if(fp == NULL){
            printf("the input .mtx file is not valid\n");
            return 0;
        }
        int i;
        int total_row = 0;
        int row, col, value;
        node** matrix;
        for(i=0; i<3; i++){
            fgets(comments, 1000, fp);
        }
        s= strtok(comments, " ");
        total_row=atoi(s);
        matrix = (node**)malloc(sizeof(node*)*total_row);

        while(fgets(line, 100, fp) != NULL){
            s= strtok(line, " ");
            row = atoi(s);
            row= row - 1;
            s= strtok(NULL, " ");
            col = atoi(s);
            col= col - 1;
            s= strtok(NULL, " ");
            value = atoi(s);
            if(matrix[row] == NULL){
                if (init(&matrix[row], col, value) != 0) {
                    fprintf(stderr, "Failed to init a new linked list");
                    return 0;
                }
            }
            else{
                insert(&matrix[row], col, value);
            }
        }
        fp2 = fopen("sorted.mtx", "w+");
        for(i=0; i<total_row; i++){
            if(matrix[i] != NULL){
                reverse(&matrix[i]);
                node *current = matrix[i];
                while (current) {
                    fprintf(fp2, "%d %d %d\n", i+1, (current->column)+1, current->value);
                    current = current->next;
                }
            }
        }
        fclose(fp);
        fclose(fp2);
        //the sort ends and the compression starts here
        int count_lines = 0;
        char chr;
        fp = fopen("sorted.mtx", "r");
        if(fp == NULL){
            printf("the input file is invalid!\n");
            return 0;
        }
        //count number of lines in the file
        chr = getc(fp);
        while (chr != EOF)
        {
            //count whenever new line is encountered
            if (chr == '\n')
            {
                count_lines = count_lines + 1;
            }
            //take next character from file.
            chr = getc(fp);
        }
        fclose(fp);
        count_lines -=1;
        fp = fopen("sorted.mtx", "r");
        //run the compressor
        int block;
        if(argc == 6){
            block = 2000;
        }
        else{
            block = atoi(argv[2]);
        }
        //open the gene labels file
        fp_gene= fopen(argv[4], "r");
        sparse_compressor(fp, fp_gene, count_lines, block);
        // system("rm sparse_parsed/*");
       //create the folder specified by the user
        char command2[200];
        snprintf(command2, sizeof(command2), "mkdir %s", argv[argc-1]);
        struct stat st = {0};
        if (stat(argv[argc-1], &st) == -1) {
            system(command2);
        }
        char command1[200];
        snprintf(command1, sizeof(command1), "tar -cf %s/sparse_compressed.tar sparse_compressed", argv[argc-1]);
        system(command1);
        //compress the barcode file with BSC compressor
        char command0[200];
        snprintf(command0, sizeof(command0), "BSC/bsc e %s %s/compressed_barcodes", argv[5], argv[argc-1]);
        system(command0);
        system("rm sparse_compressed/*");
        printf("compression and linking of sparse matrix file succeeds!\n");
        fclose(fp);
        system("rm sorted.mtx");
        fclose(fp_gene);
    }
    else if(strcmp("-qs", argv[1]) == 0){
        char hash_key[500];
        char* hash_val;
        char temp[100];
        char* retval;

        FILE *fp_hash_key= fopen("index_tables/sparse_key.txt", "r");
        FILE *fp_hash_val= fopen("index_tables/sparse_value.txt", "r");
        hashtable_t *ht = ht_create(3000000);

        while(fscanf(fp_hash_key, "%s", hash_key)!=EOF){
            hash_val= (char*)malloc(sizeof(char)*50);
            fgets(hash_val, 50, fp_hash_val);
            hash_val[strlen(hash_val) - 1] = '\0';
            ht_put(ht, hash_key, hash_val);
        }
        //close the files
        fclose(fp_hash_key);
        fclose(fp_hash_val);
        char command3[200];
        snprintf(command3, sizeof(command3), "tar -xf %s/sparse_compressed.tar sparse_compressed", argv[argc-1]);
        system(command3);
        char* hashval;
        char* s;
        int block;
        int block_id;
        int block_start_id;
        int block_end_id;
        hashval= (char*)malloc(sizeof(char)*100);
        hashval= (char*)ht_get(ht, argv[2]);
        if(hashval == NULL){
            fprintf(stderr, "This ID is not valid!\n");
            return 0;
        }
        s= (char*)malloc(sizeof(char)*50);
        s= strtok(hashval, " ");
        block= atoi(s);
        s= strtok(NULL, " ");
        block_start_id= atoi(s);
        s= strtok(NULL, " ");
        block_end_id= atoi(s);
        //decompress the barcodes file
        char command0[200];
        snprintf(command0, sizeof(command0), "BSC/bsc d %s/compressed_barcodes search_barcodes.tsv", argv[argc-1]);
        system(command0);
        sparseSearch(block, block_start_id, block_end_id);
        system("rm search_barcodes.tsv");
        //check if GTF file contains extra information
        fp_hash_key= fopen("index_tables/data_key.txt", "r");
        fp_hash_val= fopen("index_tables/data_value.txt", "r");
        int exist = 0;
        while(fscanf(fp_hash_key, "%s", hash_key)!=EOF){
            fgets(hash_val, 50, fp_hash_val);
            if(!strcmp(hash_key, argv[2])){
                exist = 1;
                break;
            }
        }
        if(exist == 1){
            char command4[200];
            snprintf(command4, sizeof(command4), "tar -xf %s/GTF_compressed.tar GTF_compressed", argv[argc-1]);
            system(command4);
            s= strtok(hash_val, " ");
            block= atoi(s);
            s= strtok(NULL, " ");
            block_id= atoi(s);
            retval = item_search(block, block_id);
            printf("The item with this id also exists in GFF file:\n");
            printf("%s", retval);
            system("rm GTF_compressed/*");
            system("rm GTF_parsed/*");
        }
        system("rm sparse_compressed/*");
        system("rm sparse_parsed/*");
        fprintf(stderr, "sparse search succeeds!\n");

    }
    else if(strcmp("-qer", argv[1]) == 0){
    	char hash_key[500];
        char* hash_val;
        char temp[100];

        //recover all the data structures
        int idx = 5;
        char command1[200];
        char command2[200];
        char command3[200];
        snprintf(command1, sizeof(command1), "BSC/bsc d %s/data_key_compressed index_tables/data_key.txt", argv[idx]);
        system(command1);
        snprintf(command2, sizeof(command2), "BSC/bsc d %s/data_value_compressed index_tables/data_value.txt", argv[idx]);
        system(command2);
        snprintf(command3, sizeof(command3), "BSC/bsc d %s/data_chr_compressed index_tables/data_chr.txt", argv[idx]);
        system(command3);

        char command4[200];
        snprintf(command4, sizeof(command4), "tar -xf %s/GTF_compressed.tar GTF_compressed", argv[idx]);
        system(command4);
    	//decompress the index tables for min and max positions
    	char command5[200];
        char command6[200];
        snprintf(command5, sizeof(command5), "BSC/bsc d %s/data_min_compressed index_tables/data_min.txt", argv[idx]);
        system(command5);
        snprintf(command6, sizeof(command6), "BSC/bsc d %s/data_max_compressed index_tables/data_max.txt", argv[idx]);
        system(command6);
		//decompress the sparse matrix files
		char command8[200];
	    snprintf(command8, sizeof(command8), "tar -xf %s/expression_compressed.tar expression_compressed", argv[argc-1]);
	    system(command8);
        //create the position min and max tables
        int *min_table= (int*)malloc(sizeof(int)*500);
        int *max_table= (int*)malloc(sizeof(int)*500);
    	//recover the chromosome table
        FILE *fp_chromosome = fopen("index_tables/data_chr.txt", "rb");
        fread(chr_table, sizeof(int), sizeof(chr_table), fp_chromosome);
    	//recover the block min position table
        FILE *fp_min = fopen("index_tables/data_min.txt", "rb");
        fread(min_table, sizeof(int), sizeof(min_table), fp_min);
    	//recover the block max position table
        FILE *fp_max = fopen("index_tables/data_max.txt", "rb");
        fread(max_table, sizeof(int), sizeof(max_table), fp_max);
        //close the files
        fclose(fp_chromosome);
        fclose(fp_min);
        fclose(fp_max);
        rangeSearch_expression(atoi(argv[2]), atoi(argv[3]), atoi(argv[4])-1, chr_table, min_table, max_table);
        //free the tables
        free(chr_table);
        free(min_table);
        free(max_table);
        fprintf(stderr, "range search for expression file succeeds!\n");
        system("rm GTF_compressed/*");
        system("rm GTF_parsed/*");

    }
    else if(strcmp("-qsr", argv[1]) == 0){
    	char hash_key[500];
        char* hash_val;
        char temp[100];

        //recover all the data structures
        int idx = 5;
        char command1[200];
        char command2[200];
        char command3[200];
        snprintf(command1, sizeof(command1), "BSC/bsc d %s/data_key_compressed index_tables/data_key.txt", argv[idx]);
        system(command1);
        snprintf(command2, sizeof(command2), "BSC/bsc d %s/data_value_compressed index_tables/data_value.txt", argv[idx]);
        system(command2);
        snprintf(command3, sizeof(command3), "BSC/bsc d %s/data_chr_compressed index_tables/data_chr.txt", argv[idx]);
        system(command3);

        char command4[200];
        snprintf(command4, sizeof(command4), "tar -xf %s/GTF_compressed.tar GTF_compressed", argv[idx]);
        system(command4);
    	//decompress the index tables for min and max positions
    	char command5[200];
        char command6[200];
        snprintf(command5, sizeof(command5), "BSC/bsc d %s/data_min_compressed index_tables/data_min.txt", argv[idx]);
        system(command5);
        snprintf(command6, sizeof(command6), "BSC/bsc d %s/data_max_compressed index_tables/data_max.txt", argv[idx]);
        system(command6);
        //decompress the barcodes file
		char command7[200];
		snprintf(command7, sizeof(command7), "BSC/bsc d %s/compressed_barcodes search_barcodes.tsv", argv[argc-1]);
		system(command7);
		//decompress the sparse matrix files
		char command8[200];
	    snprintf(command8, sizeof(command8), "tar -xf %s/sparse_compressed.tar sparse_compressed", argv[argc-1]);
	    system(command8);
        //create the position min and max tables
        int *min_table= (int*)malloc(sizeof(int)*500);
        int *max_table= (int*)malloc(sizeof(int)*500);
    	//recover the chromosome table
        FILE *fp_chromosome = fopen("index_tables/data_chr.txt", "rb");
        fread(chr_table, sizeof(int), sizeof(chr_table), fp_chromosome);
    	//recover the block min position table
        FILE *fp_min = fopen("index_tables/data_min.txt", "rb");
        fread(min_table, sizeof(int), sizeof(min_table), fp_min);
    	//recover the block max position table
        FILE *fp_max = fopen("index_tables/data_max.txt", "rb");
        fread(max_table, sizeof(int), sizeof(max_table), fp_max);
        //close the files
        fclose(fp_chromosome);
        fclose(fp_min);
        fclose(fp_max);
        rangeSearch_sparse(atoi(argv[2]), atoi(argv[3]), atoi(argv[4])-1, chr_table, min_table, max_table);

        //free the tables
        free(chr_table);
        free(min_table);
        free(max_table);
        fprintf(stderr, "range search for sparse matrix succeeds!\n");
        system("rm GTF_compressed/*");
        system("rm GTF_parsed/*");
        system("rm search_barcodes.tsv");

    }

    return 0;
}
