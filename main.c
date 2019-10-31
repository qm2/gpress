
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compressor.h"
#include "hash.h"
#include "randomaccess.h"
#include <time.h>

int main(int argc , char **argv){
    int *chr_table= (int*)malloc(sizeof(int)*100);
    // if(argc < 5){
    //     //no arguments were passed
    //     printf("not enough arguments are passed!\n");
    //     return 0;
    // } 
    if (strcmp("-c", argv[1]) == 0){

        //file pointer for gtf file
        FILE *fp, *fp_chromosome;
        int count_lines = 0;
        char chr;
        fp = fopen(argv[2], "r");
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
        count_lines -=5;
        fp = fopen(argv[2], "r");
        //run the compressor
        gtf_compressor(fp, count_lines, chr_table, atoi(argv[3]));
        fp_chromosome= fopen("index_tables/data_chr.txt", "w+");
        fwrite(chr_table , sizeof(int) , sizeof(chr_table) , fp_chromosome);
        fclose(fp_chromosome);

        system("tar -cvf GTF_compressed.tar GTF_compressed");

        system("tar -cvf index_tables.tar index_tables");   
        system("BSC/bsc e index_tables.tar index_tables_compressed.txt");     
        
        fclose(fp);
    }
    else if(strcmp("-uc", argv[1]) == 0){
        clock_t start, end;
        double cpu_time_used;   
        start = clock();

        char hash_key[500]; 
        char* hash_val;

        char temp[100];
     //    int *chr_table= (int*)malloc(sizeof(int)*100);




     //    //recover all the data structures 
        system("/mnt/e/bsc/bsc d index_tables_compressed.txt index_tables.rar");
        system("tar -xvf index_tables.rar");
        FILE *fp_hash_key= fopen("index_tables/data_key.txt", "r");
        FILE *fp_hash_val= fopen("index_tables/data_value.txt", "r");


        hashtable_t *ht = ht_create(3000000);

        while(fscanf(fp_hash_key, "%s", hash_key)!=EOF){
            hash_val= (char*)malloc(sizeof(char)*50);
            fgets(hash_val, 50, fp_hash_val);
            hash_val[strlen(hash_val) - 1] = '\0';
            ht_put(ht, hash_key, hash_val); 
        }

        //recover the chromosome table
        FILE *fp_chromosome = fopen("index_tables/data_chr.txt", "rb"); 
        fread(chr_table, sizeof(int), sizeof(chr_table), fp_chromosome);

        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf("the link time is: %lf(s)\n", cpu_time_used);

        //close the files
        fclose(fp_hash_key);
        fclose(fp_hash_val);
        fclose(fp_chromosome);
        // system("tar -cvf results.tar results");
        printf("Press enter to continue\n");
        char enter = 0;
        while (1) { 
             enter = getchar();
             if(enter == 'i'){
                    char* retval;
                    char* hashval;
                    char* s;
                    int block;
                    int block_id;
                    hashval= (char*)malloc(sizeof(char)*100);
                    hashval= (char*)ht_get(ht, "ENSMUSE00000808726.1");
                    printf("%s\n", hashval);
                    s= (char*)malloc(sizeof(char)*50);
                    s= strtok(hashval, " ");
                    block= atoi(s);
                    s= strtok(NULL, " ");
                    block_id= atoi(s);
                    retval = item_search(block, block_id);
                    printf("%s\n", retval);

             } 

             else if(enter == 's'){
                 rangeSearch(5298887, 5565649, 1,  chr_table);
             }
             else if(enter == 'm'){
                 rangeSearch(5298887, 72384131, 1,  chr_table);
             }
             else if(enter == 'l'){
                 rangeSearch(5298887, 127780063, 1,  chr_table);
             }
             printf("random access succeed\n");
            }

    }
    else if(strcmp("-a", argv[1]) == 0){
        int total=0;
        double num;

        FILE* fp=fopen("test_c5.txt", "r");
        FILE* fp_new=fopen("test_c5_new.txt", "w");
        FILE* fp_idx=fopen("test_c5_idx.txt", "w");
        while(fscanf(fp, "%lf", &num)!=EOF){
            total++;
            if(num == 0){
                fprintf(fp_idx, "%d\n", total);
            }
            else{
                fprintf(fp_new, "%lf\n", num);
            }
        }
        fclose(fp);
        fclose(fp_new);
        fclose(fp_idx);
        // add_database_id("ENSE00003462276.1", "FLYBASW123434343");
    }
    else if(strcmp("-emc", argv[1]) == 0){
         //file pointer for gtf file
        FILE *fp;
        int count_lines = 0;
        char chr;
        fp = fopen(argv[3], "r");
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
        matrix_compressor(fp, count_lines, atoi(argv[4]));      
        fclose(fp);
    }
    else if(strcmp("-emuc", argv[1]) == 0){
        char hash_key[500]; 
        char* hash_tmp;
        char* hash_val;
        char* s;
        char temp[100];
     //    int *chr_table= (int*)malloc(sizeof(int)*100);

        FILE *fp_hash_key= fopen("index_tables/matrix_key.txt", "r");
        FILE *fp_hash_val= fopen("index_tables/matrix_value.txt", "r");
        hash_val= (char*)malloc(sizeof(char)*100);
        fgets(hash_val, 100, fp_hash_val);
        free(hash_val);
        hashtable_t *ht = ht_create(3000000);
        s= (char*)malloc(sizeof(char)*50);
        s= strcpy(s,"0");
        while(fscanf(fp_hash_key, "%s", hash_key)!=EOF){
            hash_val= (char*)malloc(sizeof(char)*100);
            hash_tmp= (char*)malloc(sizeof(char)*200);            
            fgets(hash_val, 100, fp_hash_val);
            hash_val[strlen(hash_val) - 1] = ' ';
            strcpy(hash_tmp, hash_val);
            strcat(hash_tmp, s);
            // printf("%s  %s\n", hash_key, hash_tmp);
            ht_put(ht, hash_key, hash_tmp);
            s= (char*)malloc(sizeof(char)*50);
            s= strtok(hash_val, " ");
            s= strtok(NULL, " ");
        }
        char* hashval;
        hashval= (char*)malloc(sizeof(char)*100);
        hashval= (char*)ht_get(ht, "ENSTR0000623253.3");
        printf("%s\n", hashval);

    }
    else if(strcmp("-smc", argv[1]) == 0){
        //file pointer for gtf file
        FILE *fp;
        int count_lines = 0;
        char chr;
        fp = fopen(argv[3], "r");
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
        count_lines -=3;
        fp = fopen(argv[3], "r");
        //run the compressor
        sparse_matrix_compressor(fp, count_lines, atoi(argv[4]));       
        fclose(fp);

    }
    return 0;
}