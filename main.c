
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compressor.h"
#include "hash.h"
#include "randomaccess.h"
#include <time.h>

int main(int argc , char **argv){
    int *chr_table= (int*)malloc(sizeof(int)*100);
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
        fp = fopen(argv[2], "r");
        char *dot = strrchr(argv[2], '.');
        if(!strcmp(dot+1, "gtf")){
            count_lines -=5;    
            gtf_compressor(fp, count_lines, chr_table, atoi(argv[3]));
        } 
        else if(!strcmp(dot+1, "gff3")){
            count_lines -=7;    
            gff3_compressor(fp, count_lines, chr_table, atoi(argv[3]));
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

        system("rm GTF_parsed/*");
        system("tar -cvf compressed/GTF_compressed.tar GTF_compressed");
        system("rm GTF_compressed/*"); 
        system("BSC/bsc e index_tables/data_key.txt compressed/data_key_compressed");
        system("BSC/bsc e index_tables/data_value.txt compressed/data_value_compressed");
        system("BSC/bsc e index_tables/data_chr.txt compressed/data_chr_compressed");             
   
        printf("The compression of GFF file with random access succeeds!\n");
        fclose(fp);
    }
    else if (strcmp("-cw", argv[1]) == 0){

        //file pointer for gtf file
        FILE *fp;
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
        system("tar -xvf compressed/GTF_compressed_without.tar GTF_compressed2");
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
            fp = fopen("decompressed_gtf.gtf", "w+");
            gtf_decompressor(fp, count_lines, 0); 
            printf("The decompressed GTF file is included in the decompressed_gtf.gtf!\n");           
        }
        else if(!strcmp(argv[2], "gff3")){
            fp = fopen("decompressed_gff3.gff3", "w+");
            gtf_decompressor(fp, count_lines, 1); 
            printf("The decompressed GFF3 file is included in the decompressed_gff3.gff3!\n");          
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
        system("BSC/bsc d compressed/data_key_compressed index_tables/data_key.txt");
        system("BSC/bsc d compressed/data_value_compressed index_tables/data_value.txt");
        system("BSC/bsc d compressed/data_chr_compressed index_tables/data_chr.txt");  
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

        //close the files
        fclose(fp_hash_key);
        fclose(fp_hash_val);
        fclose(fp_chromosome);
        system("tar -xvf compressed/GTF_compressed.tar GTF_compressed");

        if(strcmp("-id", argv[2]) == 0){
               char* retval;
               char* hashval;
               char* s;
               int block;
               int block_id;
               hashval= (char*)malloc(sizeof(char)*100);
               hashval= (char*)ht_get(ht, argv[3]);
               if(hashval == NULL){
                   printf("This ID is not valid!\n");
                   return 0;
               }
               s= (char*)malloc(sizeof(char)*50);
               s= strtok(hashval, " ");
               block= atoi(s);
               s= strtok(NULL, " ");
               block_id= atoi(s);
               retval = item_search(block, block_id);
               printf("The item with id %s is:\n",argv[3]);
               printf("%s", retval);
               printf("id search succeeds!\n");
               system("rm GTF_compressed/*");
               system("rm GTF_parsed/*");
        }
        else if(strcmp("-range", argv[2]) == 0){
            rangeSearch(atoi(argv[3]), atoi(argv[4]), atoi(argv[5])-1,  chr_table);
            printf("All items on chromosome %s from %s to %s are outputed in the range_search.gtf file\n",argv[5], argv[3], argv[4]);
            printf("range search succeeds!\n");
            system("rm GTF_compressed/*");
            system("rm GTF_parsed/*");
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
    else if(strcmp("-e", argv[1]) == 0){
         //file pointer for gtf file
        FILE *fp;
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
        count_lines -=1;
        fp = fopen(argv[2], "r");
        //run the compressor
        expression_compressor(fp, count_lines, atoi(argv[3]));  
        system("rm expression_parsed/*");
        system("tar -cvf compressed/expression_compressed.tar expression_compressed");
        system("rm expression_compressed/*"); 
        printf("compression and linking of expression file succeeds!\n");
        fclose(fp);
    }
    else if(strcmp("-qe", argv[1]) == 0){
        char hash_key[500]; 
        char* hash_val;
        char temp[100];
        char* retval;

        system("BSC/bsc d compressed/data_key_compressed index_tables/data_key.txt");
        system("BSC/bsc d compressed/data_value_compressed index_tables/data_value.txt");
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
        system("tar -xvf compressed/expression_compressed.tar expression_compressed");
        char* hashval;
        char* s;
        int block;
        int block_id;
        int block_start_id;
        int block_end_id;
        hashval= (char*)malloc(sizeof(char)*100);
        hashval= (char*)ht_get(ht, argv[2]);
        if(hashval == NULL){
            printf("This ID is not valid!\n");
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
            system("tar -xvf compressed/GTF_compressed.tar GTF_compressed");
            s= strtok(hash_val, " ");
            block= atoi(s);
            s= strtok(NULL, " ");
            block_id= atoi(s);
            retval = item_search(block, block_id);
            printf("All the information of item with id %s is outputed in expression_search.txt:\n",argv[2]);
            printf("The item with id %s also exists in GFF file:\n",argv[2]);
            printf("%s", retval);
            system("rm GTF_compressed/*");
            system("rm GTF_parsed/*");
        }
        else{
            printf("All the information of item with id %s is outputed in expression_search.txt:\n",argv[2]);
            printf("The item with this id does not exist in the compressed GFF file\n");
        }
        system("rm expression_compressed/*");
        system("rm expression_parsed/*");
        printf("expression search succeeds!\n");
    }

    return 0;
}