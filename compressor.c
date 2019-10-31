#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compressor.h"
#include "hash.h"
#define BUFFSIZE 1000



int gtf_compressor(FILE* fp, int length, int* chr_table, int block_size){
    int i, j, k, m;
    char empty = ' ';
    char prev_chr[100];
    int new_transcript =0, new_block =1;
    int len;
    int block, gene_numbers=0, item_id;
    //create file pointers for all the output files
    FILE *fp_sq, *fp_src, *fp_fea, *fp_start, *fp_delta, *fp_att, *fp_comments;
    FILE *fp_score, *fp_strand, *fp_frame_cds, *fp_frame_start, *fp_frame_stop, *fp_hash_key, *fp_hash_val; 
    //array to hold each column
    char comments[BUFFSIZE];
    char seqname[BUFFSIZE];
    char source[BUFFSIZE];
    char feature[BUFFSIZE];
    char start[BUFFSIZE];
    char new_start[BUFFSIZE];
    char end[BUFFSIZE];
    char delta[BUFFSIZE];
    char score[BUFFSIZE];
    char strand[BUFFSIZE];
    char frame[BUFFSIZE];
    char attribute[BUFFSIZE];
    char refer[BUFFSIZE];
    char group_id[BUFFSIZE];
    char* group;
    char* id;

    int prev_gene =0, prev_trans =0, prev_exon =0;
    char* gene ="gene";
    char* transcript ="transcript";
    char* exon = "exon";
    char* cds = "CDS";
    char* start_codon = "start_codon";
    char* stop_codon = "stop_codon";
    char* gene_id = "gene_id";
    char* transcript_id = "transcript_id";
    char* exon_id = "exon_id";

    char sq_name[100]; 
    char src_name[100]; 
    char fea_name[100]; 
    char start_name[100];       
    char delta_name[100]; 
    char score_name[100];
    char strand_name[100]; 
    char att_name[100];
    char frame_cds_name[100]; 
    char frame_start_name[100]; 
    char frame_stop_name[100];
    char compress_cmd[100];
    char* sq_name_prefix= "GTF_parsed/gtf_seqname_" ;
    char* src_name_prefix= "GTF_parsed/gtf_source_";
    char* fea_name_prefix= "GTF_parsed/gtf_feature_";
    char* start_name_prefix= "GTF_parsed/gtf_start_";       
    char* delta_name_prefix= "GTF_parsed/gtf_delta_";
    char* score_name_prefix= "GTF_parsed/gtf_score_";
    char* strand_name_prefix= "GTF_parsed/gtf_strand_";
    char* att_name_prefix= "GTF_parsed/gtf_attribute_";
    char* frame_cds_name_prefix= "GTF_parsed/gtf_frame_cds_";
    char* frame_start_name_prefix= "GTF_parsed/gtf_frame_start_";
    char* frame_stop_name_prefix= "GTF_parsed/gtf_frame_stop_";
    char* compress_cmd_prefix= "BSC/bsc e ";
    char* cmd_suffix_sq= " GTF_compressed/compressed_seqname_";
    char* cmd_suffix_src=" GTF_compressed/compressed_source_";    
    char* cmd_suffix_fea=" GTF_compressed/compressed_feature_";
    char* cmd_suffix_start=" GTF_compressed/compressed_start_";
    char* cmd_suffix_delta=" GTF_compressed/compressed_delta_";
    char* cmd_suffix_score=" GTF_compressed/compressed_score_";
    char* cmd_suffix_strand=" GTF_compressed/compressed_strand_";
    char* cmd_suffix_att=" GTF_compressed/compressed_attribure_";
    char* cmd_suffix_frame_cds=" GTF_compressed/compressed_frame_cds_";
    char* cmd_suffix_frame_start=" GTF_compressed/compressed_frame_start_";
    char* cmd_suffix_frame_stop=" GTF_compressed/compressed_frame_stop_";


    char block_number[100];



    //extract each column and write them into the output files
    block= 0;
    sprintf(block_number,"%d", block);
    //create the output files
    strcpy(sq_name, sq_name_prefix);
    strcpy(src_name, src_name_prefix);
    strcpy(fea_name, fea_name_prefix);
    strcpy(start_name, start_name_prefix);
    strcpy(delta_name, delta_name_prefix);
    strcpy(score_name, score_name_prefix);
    strcpy(strand_name, strand_name_prefix);
    strcpy(att_name, att_name_prefix);
    strcpy(frame_cds_name, frame_cds_name_prefix);
    strcpy(frame_start_name, frame_start_name_prefix);
    strcpy(frame_stop_name, frame_stop_name_prefix);

    strcat(strcat(sq_name, block_number),".txt");
    strcat(strcat(src_name, block_number),".txt");
    strcat(strcat(fea_name, block_number),".txt");
    strcat(strcat(start_name, block_number),".txt");
    strcat(strcat(delta_name, block_number),".txt");
    strcat(strcat(score_name, block_number),".txt");
    strcat(strcat(strand_name, block_number),".txt");
    strcat(strcat(att_name, block_number),".txt");
    strcat(strcat(frame_cds_name, block_number),".txt");
    strcat(strcat(frame_start_name, block_number),".txt");
    strcat(strcat(frame_stop_name, block_number),".txt");

    //write the seqname to a file
    fp_sq = fopen(sq_name, "w+");
    //write the source to a file
    fp_src = fopen(src_name, "w+");
    //write the feature to a file
    fp_fea = fopen(fea_name, "w+");
    //write the start to a file
    fp_start = fopen(start_name, "w+");
    //write the difference of start and stop to a file
    fp_delta = fopen(delta_name, "w+");
    //write the score to a file
    fp_score = fopen(score_name, "w+");
    //write the strand to a file
    fp_strand = fopen(strand_name, "w+");
    //write the attribute to a file
    fp_att = fopen(att_name, "w+");
    //write the frames of CDS to a file
    fp_frame_cds = fopen(frame_cds_name, "w+");
    //write the frames of start to a file
    fp_frame_start = fopen(frame_start_name, "w+");
    //write the frames of stop to a file
    fp_frame_stop = fopen(frame_stop_name, "w+");

    fp_hash_key= fopen("index_tables/data_key.txt", "w+");
    fp_hash_val= fopen("index_tables/data_value.txt", "w+");
    //store the first five lines of comments
    for(i=0; i<5; i++){
        fgets(comments, BUFFSIZE, fp);
        //fprintf(fp_att, "%s", comments);           
    }
    item_id = -1;
    strcpy(prev_chr, "chr1");
    int c = 0;
    chr_table[c] = 0;
    for(i=0; i< length; i++){
        item_id++;
        //read in the seqname
        fscanf(fp, "%s", seqname);

        //read in the source
        fscanf(fp, "%s", source);

        //read in the feature
        fscanf(fp, "%s", feature);

        //check if we need to update the block 
        if(!strcmp(feature, gene)){
            gene_numbers++;
            if(gene_numbers == block_size){
                block++;
                new_block= 1;
                gene_numbers=0;
                item_id=0;
                fclose(fp_sq);
                fclose(fp_src);
                fclose(fp_fea);
                fclose(fp_start);
                fclose(fp_delta);
                fclose(fp_score);
                fclose(fp_strand);
                fclose(fp_att);
                fclose(fp_frame_cds);
                fclose(fp_frame_start);
                fclose(fp_frame_stop);
                //compress the files using bsc algorithm
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, sq_name), cmd_suffix_sq);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, src_name), cmd_suffix_src);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, fea_name), cmd_suffix_fea);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, start_name), cmd_suffix_start);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, delta_name), cmd_suffix_delta);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, score_name), cmd_suffix_score);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, strand_name), cmd_suffix_strand);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, att_name), cmd_suffix_att);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, frame_cds_name), cmd_suffix_frame_cds);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, frame_start_name), cmd_suffix_frame_start);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, frame_stop_name), cmd_suffix_frame_stop);
                strcat(compress_cmd, block_number);
                system(compress_cmd);

                sprintf(block_number,"%d", block);
                //create the output files
                strcpy(sq_name, sq_name_prefix);
                strcpy(src_name, src_name_prefix);
                strcpy(fea_name, fea_name_prefix);
                strcpy(start_name, start_name_prefix);
                strcpy(delta_name, delta_name_prefix);
                strcpy(score_name, score_name_prefix);
                strcpy(strand_name, strand_name_prefix);
                strcpy(att_name, att_name_prefix);
                strcpy(frame_cds_name, frame_cds_name_prefix);
                strcpy(frame_start_name, frame_start_name_prefix);
                strcpy(frame_stop_name, frame_stop_name_prefix);
                strcat(strcat(sq_name, block_number),".txt");
                strcat(strcat(src_name, block_number),".txt");
                strcat(strcat(fea_name, block_number),".txt");
                strcat(strcat(start_name, block_number),".txt");
                strcat(strcat(delta_name, block_number),".txt");
                strcat(strcat(score_name, block_number),".txt");
                strcat(strcat(strand_name, block_number),".txt");
                strcat(strcat(att_name, block_number),".txt");
                strcat(strcat(frame_cds_name, block_number),".txt");
                strcat(strcat(frame_start_name, block_number),".txt");
                strcat(strcat(frame_stop_name, block_number),".txt");

                //write the seqname to a file
                fp_sq = fopen(sq_name, "w+");
                //write the source to a file
                fp_src = fopen(src_name, "w+");
                //write the feature to a file
                fp_fea = fopen(fea_name, "w+");
                //write the start to a file
                fp_start = fopen(start_name, "w+");
                //write the difference of start and stop to a file
                fp_delta = fopen(delta_name, "w+");
                //write the score to a file
                fp_score = fopen(score_name, "w+");
                //write the strand to a file
                fp_strand = fopen(strand_name, "w+");
                //write the attribute to a file
                fp_att = fopen(att_name, "w+");
                //write the frames of CDS to a file
                fp_frame_cds = fopen(frame_cds_name, "w+");
                //write the frames of start to a file
                fp_frame_start = fopen(frame_start_name, "w+");
                //write the frames of stop to a file
                fp_frame_stop = fopen(frame_stop_name, "w+");
            }
        }
        //update the chromsome table
        if(strcmp(seqname, prev_chr)){
            strcpy(prev_chr, seqname);
            c++;
            chr_table[c] = block;
        }
        //output the seqname
        fprintf(fp_sq, "%s\n", seqname);
        //output the source
        fprintf(fp_src, "%s\n", source);
        //output the feature
        fprintf(fp_fea, "%s\n", feature);

        //read in the start
        fscanf(fp, "%s", start);
        if(new_block == 1){
            //position[block]= atoi(start);
            new_block= 0;
        }
        //read in the end
        fscanf(fp, "%s", end);

        //calculate the difference between start and stop
        sprintf(delta,"%d", atoi(end) - atoi(start));
        //output the delta
        fprintf(fp_delta, "%s\n", delta);  

        //read in the score
        fscanf(fp, "%s", score);
        fprintf(fp_score, "%s\n", score);

        //read in the strand
        fscanf(fp, "%s", strand);
        //check if current item is gene or not
        if(!strcmp(feature, gene)){
            //output the strands of gene
            fprintf(fp_strand, "%s\n", strand);
        }

        //read in the frame
        fscanf(fp, "%s", frame);
        //check if the current item has frame or not
        //the current item is CDS
        if(!strcmp(feature, cds)){ 
            fprintf(fp_frame_cds,"%s\n", frame);
        }
        //the current item is start codon
        else if(!strcmp(feature, start_codon)){
            fprintf(fp_frame_start,"%s\n", frame);
        }       
        else if(!strcmp(feature, stop_codon)){
            fprintf(fp_frame_stop,"%s\n", frame);
        } 
        //skip the empty spaces
        if((empty=fgetc(fp)) == ' ')
        { 
            empty = fgetc(fp);
        }
        //read in the attribute
        fgets(attribute, BUFFSIZE, fp);
        //ouput the attribute
        fprintf(fp_att, "%s", attribute);
        

        //extract the id from the attribute
        len= strlen(attribute);
        free(id);
        id =(char*) malloc(100*sizeof(char));
        id[0]='\0';
        if(strcmp(feature, gene)  == 0){ 
            //find the starting position of gene id
            memset(refer,0,sizeof(refer));
            m=0;
            for(j= 0; j< len; j++){
                if(attribute[j]==' '){
                    if(!strcmp(refer, gene_id)){
                        break;
                    }
                    m=0;
                    memset(refer,0,sizeof(refer));
                }
                else{
                    refer[m] = attribute[j];  
                    m++;                
                }               
            }
            j+=2;
            for(k=0; attribute[j+k]!='"'; k++){
                id[k] = attribute[j+k];
            }
            id[k]='\0';
            // hash the id
            free(group);
            group= (char*)malloc(sizeof(char)*100);
            group[0]= '\0';
            memset(group_id,0,sizeof(group_id));
            sprintf(group,"%d", block);
            group[strlen(group)]=' ';
            sprintf(group_id, "%d", item_id);
            // if(item_id==0){
            //     printf("group:%s group_id:%s block:%d\n", group, group_id, block);
            // }  
            strcat(group, group_id);    
            fprintf(fp_hash_key, "%s\n", id); 
            if(fprintf(fp_hash_val, "%s\n", group)<0){
                printf("%s\n", group);
            } 
            //  if(item_id==0){
            //     printf("total group:%s\n", group);
            // }  

        }
        else if(strcmp(feature, transcript) == 0){
            //find the starting position of transcript id
            memset(refer,0,sizeof(refer));
            m=0;
            for(j= 0; j< len; j++){
                if(attribute[j]==' '){
                    if(!strcmp(refer, transcript_id)){
                        break;
                    }
                    m=0;
                    memset(refer,0,sizeof(refer));
                }
                else{
                    refer[m] = attribute[j];  
                    m++;                
                }               
            }
            j+=2;
            for(k=0; attribute[j+k]!='"'; k++){
                id[k] = attribute[j+k];
            }
            id[k]='\0';
            // hash the id
            free(group);
            group= (char*)malloc(sizeof(char)*100);
            group[0]= '\0';
            memset(group_id,0,sizeof(group_id));
            sprintf(group,"%d", block);  
            group[strlen(group)]=' ';
            sprintf(group_id, "%d", item_id);
            strcat(group, group_id);    
            fprintf(fp_hash_key, "%s\n", id); 
            if(fprintf(fp_hash_val, "%s\n", group)<0){
                printf("%s\n", group);
            } 
        }
        else if(strcmp(feature, exon) == 0){
            //find the starting position of exon id
            memset(refer,0,sizeof(refer));
            m=0;
            for(j= 0; j< len; j++){
                if(attribute[j]==' '){
                    if(!strcmp(refer, exon_id)){
                        break;
                    }
                    m=0;
                    memset(refer,0,sizeof(refer));
                }
                else{
                    refer[m] = attribute[j];  
                    m++;                
                }               
            }
            j+=2;
            for(k=0; attribute[j+k]!='"'; k++){
                id[k] = attribute[j+k];
            }
            id[k]='\0';
            // hash the id
            free(group);
            group= (char*)malloc(sizeof(char)*100);
            group[0]= '\0';
            memset(group_id,0,sizeof(group_id));
            sprintf(group,"%d", block);    
            group[strlen(group)]=' ';
            sprintf(group_id, "%d", item_id);
            strcat(group, group_id);    
            fprintf(fp_hash_key, "%s\n", id); 
            if(fprintf(fp_hash_val, "%s\n", group)<0){
                printf("%s\n", group);
            } 
        }


        //modify the start for better compression
        if(strcmp(feature, gene)  == 0){ 
            //store the gene start for later uses
            prev_gene = atoi(start);
            prev_trans = prev_gene;
            //write to the new start 
            sprintf(new_start,"%d", atoi(start));
        }
        else if(strcmp(feature, transcript) == 0){
            //store the transcript start for later uses
            //write to the new start 
            sprintf(new_start,"%d", atoi(start) - prev_trans);
            prev_trans = atoi(start);
            prev_exon = prev_trans;
            new_transcript = 1;
        }
        else if(strcmp(feature, exon) == 0){
            //write to the new start 
            //check the strand
            if(strcmp(strand, "+") == 0 || new_transcript == 1){
               sprintf(new_start,"%d", atoi(start) - prev_exon); 
            }
            else{
               sprintf(new_start,"%d", prev_exon - atoi(start)); 
            }   
            prev_exon = atoi(start);    
            new_transcript = 0;    
        }
        else{
            sprintf(new_start,"%d", atoi(start) - prev_trans);
        }
        //output the new start
        fprintf(fp_start, "%s\n", new_start);
    }
    printf("%d\n", block);
    //close all the files
    fclose(fp_sq);
    fclose(fp_src);
    fclose(fp_fea);
    fclose(fp_start);
    fclose(fp_delta);
    fclose(fp_score);
    fclose(fp_strand);
    fclose(fp_att);
    fclose(fp_frame_cds);
    fclose(fp_frame_start);
    fclose(fp_frame_stop);
    return 0;
}

int gtf_decompressor(FILE* fp, int length){
    char seqname[BUFFSIZE];
    char source[BUFFSIZE];
    char feature[BUFFSIZE];
    char start[BUFFSIZE];
    char new_start[BUFFSIZE];
    char end[BUFFSIZE];
    char delta[BUFFSIZE];
    char score[BUFFSIZE];
    char strand[BUFFSIZE];
    char frame[BUFFSIZE];
    char attribute[BUFFSIZE];
    char comments[BUFFSIZE];
    char* gene ="gene";
    char* transcript ="transcript";
    char* exon = "exon";
    char* cds = "CDS";
    char* start_codon = "start_codon";
    char* stop_codon = "stop_codon";

    int prev_gene =0, prev_trans =0, prev_exon =0;
    int new_transcript = 0;

     //create file pointers for all the files
    FILE *fp_sq, *fp_src, *fp_fea, *fp_start, *fp_delta, *fp_att, *fp_comments;
    FILE *fp_score, *fp_strand, *fp_frame_cds, *fp_frame_start, *fp_frame_stop; 
    //open all the files
    fp_sq = fopen("GTF_parsed/gtf_seqname.txt", "r");
    //open the source file
    fp_src = fopen("GTF_parsed/gtf_source.txt", "r");
    //open the feature file
    fp_fea = fopen("GTF_parsed/gtf_feature.txt", "r");
    //open the start file
    fp_start = fopen("GTF_parsed/gtf_start.txt", "r");
    //open the difference of start and stop file
    fp_delta = fopen("GTF_parsed/gtf_delta.txt", "r");
    //open the score file
    fp_score = fopen("GTF_parsed/gtf_score.txt", "r");
    //open the strand file
    fp_strand = fopen("GTF_parsed/gtf_strand.txt", "r");
    //open the attribute file
    fp_att = fopen("GTF_parsed/gtf_attribute.txt", "r");
    //open the frames of CDS file
    fp_frame_cds = fopen("GTF_parsed/gtf_frame_cds.txt", "r");
    //open the frames of start file
    fp_frame_start = fopen("GTF_parsed/gtf_frame_start.txt", "r");
    //open the frames of stop file
    fp_frame_stop = fopen("GTF_parsed/gtf_frame_stop.txt", "r");

    //write the comments to the gtf file
    for(int i=0; i< 5; i++){
        fgets(comments, BUFFSIZE, fp_att);
        //fprintf(fp, "%s", comments);
    }
    //start to combine all columns into the goal file
    for(int i=0; i< length; i++){
        //read in all the rest 
        fscanf(fp_sq, "%s", seqname);
        fscanf(fp_src, "%s", source);
        fscanf(fp_fea, "%s", feature);
        fscanf(fp_start, "%s", start);
        fscanf(fp_delta, "%s", delta);
        fscanf(fp_score, "%s", score);       
        fgets(attribute, BUFFSIZE, fp_att);
        //output the seqname
        fprintf(fp, "%s    ", seqname);
        //output the source
        fprintf(fp, "%s    ", source);
        //output the feature
        fprintf(fp, "%s    ", feature);
        //recover the strand
        if(strcmp(feature, gene)  == 0) {
            fscanf(fp_strand, "%s", strand);
        }   
        //recover the start
        if(strcmp(feature, gene)  == 0){ 
            //store the gene start for later uses
            prev_gene = atoi(start);
            prev_trans = prev_gene;
            //write to the new start 
            sprintf(new_start,"%d", atoi(start));
        }
        else if(strcmp(feature, transcript) == 0){
            //store the transcript start for later uses
            //write to the new start 
            sprintf(new_start,"%d", atoi(start) + prev_trans);
            prev_trans = atoi(new_start);
            prev_exon = prev_trans;
            new_transcript = 1;
        }
        else if(strcmp(feature, exon) == 0){
            //write to the new start 
            //check the strand
            if(strcmp(strand, "+") == 0 || new_transcript == 1){
               sprintf(new_start,"%d", atoi(start) + prev_exon); 
            }
            else{
               sprintf(new_start,"%d", prev_exon - atoi(start)); 
            }   
            prev_exon = atoi(new_start);    
            new_transcript = 0;    
        }
        else{
            sprintf(new_start,"%d", atoi(start) + prev_trans);
        }    
        //output the start
        fprintf(fp, "%s    ", new_start);
        //recover the end
        sprintf(end, "%d", atoi(delta) + atoi(new_start));
        //output the end
        fprintf(fp, "%s    ", end);
        //output the score
        fprintf(fp, "%s    ", score);       
        //output the strand
        fprintf(fp, "%s    ", strand);
        //recover the frame
        if(!strcmp(feature, cds)){ 
            fscanf(fp_frame_cds, "%s", frame);
        }
        //the current item is start codon
        else if(!strcmp(feature, start_codon)){
            fscanf(fp_frame_start, "%s", frame);
        }       
        else if(!strcmp(feature, stop_codon)){
            fscanf(fp_frame_stop, "%s", frame);
        } 
        else{
            strcpy(frame, ".");
        }
        //output the frame
        fprintf(fp, "%s    ", frame);

        //output the attribute
        fprintf(fp, "%s", attribute);        
    }


}

int matrix_compressor(FILE* fp, int length, int block_size){
    int i, j, k, m;
    char empty = ' ';
    char prev_chr[100];
    int new_transcript =0, new_block =1;
    int block, gene_numbers=0, item_id, prev_id;

    //create file pointers for all the output files
    FILE *fp_ets_counts;
    FILE *fp_tpm, *fp_eff_len, *fp_len, *fp_hash_key, *fp_hash_val; 
    //array to hold each column
    char target[BUFFSIZE];
    char sample[BUFFSIZE];
    char ets_counts[BUFFSIZE];
    char tpm[BUFFSIZE];
    char eff_len[BUFFSIZE];
    char len[BUFFSIZE];
    char trash[BUFFSIZE];
    char refer[BUFFSIZE];
    char* group_id=(char*)malloc(sizeof(char)*100);
    char* group=(char*)malloc(sizeof(char)*100);
    char* id;

    int prev_gene =0, prev_trans =0, prev_exon =0;

    char sample_name[100]; 
    char ets_counts_name[100]; 
    char tpm_name[100]; 
    char eff_len_name[100];       
    char len_name[100]; 
    char compress_cmd[100];

    char* sample_name_prefix= "matrix_compressed/matrix_sample_" ;
    char* ets_counts_name_prefix= "matrix_compressed/matrix_ets_counts_";
    char* tpm_name_prefix= "matrix_compressed/matrix_tpm_";
    char* eff_len_name_prefix= "matrix_compressed/matrix_eff_len_";
    char* len_name_prefix= "matrix_compressed/matrix_len_";       

    char* compress_cmd_prefix= "BSC/bsc/bsc e ";
    char* cmd_suffix_sample= " matrix_results/matrix_compressed_sample_";
    char* cmd_suffix_ets_counts=" matrix_results/matrix_compressed_ets_counts_";    
    char* cmd_suffix_tpm=" matrix_results/matrix_compressed_tpm_";
    char* cmd_suffix_eff_len=" matrix_results/matrix_compressed_eff_len_";
    char* cmd_suffix_len=" matrix_results/matrix_compressed_len_";

    char prev_sample[BUFFSIZE];

    char block_number[100];



    //extract each column and write them into the output files
    block= 0;
    sprintf(block_number,"%d", block);
    //create the output files
    strcpy(ets_counts_name, ets_counts_name_prefix);
    strcpy(tpm_name, tpm_name_prefix);
    strcpy(eff_len_name, eff_len_name_prefix);
    strcpy(len_name, len_name_prefix);


    strcat(strcat(ets_counts_name, block_number),".txt");
    strcat(strcat(tpm_name, block_number),".txt");
    strcat(strcat(eff_len_name, block_number),".txt");
    strcat(strcat(len_name, block_number),".txt");

    //write the source to a file
    fp_ets_counts = fopen(ets_counts_name, "w+");
    //write the feature to a file
    fp_tpm = fopen(tpm_name, "w+");
    //write the start to a file
    fp_eff_len = fopen(eff_len_name, "w+");
    //write the difference of start and stop to a file
    fp_len = fopen(len_name, "w+");

    fp_hash_key= fopen("index_tables/matrix_key.txt", "w+");
    fp_hash_val= fopen("index_tables/matrix_value.txt", "w+");
    //store the first line of comments
    for(i=0; i<1; i++){
        fgets(trash, BUFFSIZE, fp);
        //fprintf(fp_att, "%s", comments);           
    }
    item_id = -1;
    prev_id=0;
    for(i=0; i< length; i++){
        item_id++;
        //read in the target
        fscanf(fp, "%s", target);

        //read in the sample
        fscanf(fp, "%s", sample);
        //read in the ets_counts
        fscanf(fp, "%s", ets_counts);
        //read in the tpm   
        fscanf(fp, "%s", tpm);
        //read in the eff_len
        fscanf(fp, "%s", eff_len);
        //read in the len
        fscanf(fp, "%s", len);     
        //extract the transcript ID
        int index;
        for(index=0; index< BUFFSIZE; index++){
            if(sample[index] == '|'){
                break;
            }
        }
        sample[index] = '\0';
        //check if we need to update the block 
        if((i!=0) && (strcmp(sample, prev_sample))){
            gene_numbers++;
            if(gene_numbers == block_size){
                block++;
                new_block= 1;
                gene_numbers=0;
                item_id=0;
                fclose(fp_ets_counts);
                fclose(fp_tpm);
                fclose(fp_eff_len);
                fclose(fp_len);

                //compress the files using bsc algorithm
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, ets_counts_name), cmd_suffix_ets_counts);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, tpm_name), cmd_suffix_tpm);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, eff_len_name), cmd_suffix_eff_len);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, len_name), cmd_suffix_len);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
           
                sprintf(block_number,"%d", block);
                //create the output files
                strcpy(ets_counts_name, ets_counts_name_prefix);
                strcpy(tpm_name, tpm_name_prefix);
                strcpy(eff_len_name, eff_len_name_prefix);
                strcpy(len_name, len_name_prefix);

                strcat(strcat(ets_counts_name, block_number),".txt");
                strcat(strcat(tpm_name, block_number),".txt");
                strcat(strcat(eff_len_name, block_number),".txt");
                strcat(strcat(len_name, block_number),".txt");

                //write the source to a file
                fp_ets_counts = fopen(ets_counts_name, "w+");
                //write the feature to a file
                fp_tpm = fopen(tpm_name, "w+");
                //write the start to a file
                fp_eff_len = fopen(eff_len_name, "w+");
                //write the difference of start and stop to a file
                fp_len = fopen(len_name, "w+");

            }
            // hash the id
            free(group);
            group= (char*)malloc(sizeof(char)*100);
            group[0]= '\0';
            free(group_id);
            group_id= (char*)malloc(sizeof(char)*100);
            group_id[0]= '\0';
            sprintf(group,"%d", block);
            group[strlen(group)]=' ';
            sprintf(group_id, "%d", prev_id);
            strcat(group, group_id); 
            
            fprintf(fp_hash_key, "%s\n", sample); 
            if(fprintf(fp_hash_val, "%s\n", group)<0){
                printf("%s\n", group);
            } 
            prev_id= item_id;
        }

        //output the source
        fprintf(fp_ets_counts, "%s\n", ets_counts);
        //output the feature
        fprintf(fp_tpm, "%s\n",tpm);
        fprintf(fp_eff_len, "%s\n",eff_len);
        fprintf(fp_len, "%s\n",len);


        //read in the last column
        fgets(trash, BUFFSIZE, fp);
        strcpy(prev_sample, sample);

    }
    printf("%d\n", block);
    //close all the files
    fclose(fp_ets_counts);
    fclose(fp_tpm);
    fclose(fp_eff_len);
    fclose(fp_len);
    return 0;
}

int sparse_matrix_compressor(FILE* fp, int length, int block_size){
    int** sparse_matrix_col= (int **)malloc(89796 * sizeof(int *));
    int** sparse_matrix_val= (int **)malloc(89796 * sizeof(int *));   

    //create file pointers for all the output files
    FILE *output;
    int row, col, val;
    char tmp[BUFFSIZE];
    int i, j;

    output= fopen("sparse_matrix.txt", "w+");
    //skip the first line of comments
    for(i=0; i<3; i++){
        fgets(tmp, BUFFSIZE, fp);          
    }
    for(i=0; i< length; i++){
        //read in the row
        fscanf(fp, "%d", &row);

        //read in the column
        fscanf(fp, "%d", &col);
        //read in the value
        fscanf(fp, "%d", &val);   
        if(sparse_matrix_col[row]==NULL){
            sparse_matrix_col[row] = (int *)malloc(1 * sizeof(int));
            sparse_matrix_val[row] = (int *)malloc(1 * sizeof(int));
            sparse_matrix_col[row][0]= col;
            sparse_matrix_val[row][0]= val;
        }
        else{
            int len= sizeof( sparse_matrix_col[row] ) / sizeof( sparse_matrix_col[row][0]) ; 
            sparse_matrix_col[row] = (int *)realloc(sparse_matrix_col[row], sizeof(int)*(len));
            sparse_matrix_val[row] = (int *)realloc(sparse_matrix_val[row], sizeof(int)*(len));
            sparse_matrix_col[row][len]= col;
            sparse_matrix_val[row][len]= val;
        }
 
    }
    for(i=0; i< 89796; i++){
        if(sparse_matrix_col[i]!=NULL){
            for(j=0; j< (sizeof( sparse_matrix_col[i] ) / sizeof( sparse_matrix_col[i][0]) - 1); j++){
                fprintf(output, "%d %d\n",sparse_matrix_col[i][j], sparse_matrix_val[i][j]);                  
            }
 
        }
 
    }
    fclose(output);    

}