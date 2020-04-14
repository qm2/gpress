#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compressor.h"
#include "hash.h"
#define BUFFSIZE 1500



int gtf_compressor(FILE* fp, int length, int* chr_table, int* block_min_table, int* block_max_table, int block_size){
    int i, j, k, m;
    char empty = ' ';
    char prev_chr[100];
    int new_transcript =0, new_block =1;
    int len;
    int block, gene_numbers=0, item_id;
    int prev_gene_end = -1;
    int position_min, position_max;
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
    char* group = NULL;
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
    // strcpy(prev_chr, "chr1");
    int first_chr = 1;
    int c = 0;
    chr_table[c] = 0;
    for(i=0; i< length; i++){
        item_id++;
        //read in the seqname
        fscanf(fp, "%s", seqname);
        if(first_chr == 1){
            first_chr = 0;
            strcpy(prev_chr, seqname);           
        }

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
                prev_gene_end = -1;
                item_id=0;
                //update the min and max table
                block_min_table[block-1]= position_min;
                block_max_table[block-1]= position_max;
                //close the old files and open new files
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
        //check if we need to update the min or max
        if(gene_numbers==0 || atoi(start) < position_min){
            position_min= atoi(start);
        }
        if(gene_numbers==0 || atoi(start) > position_max){
            position_max= atoi(start);
        }
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
            for(k=0; attribute[j+k]!='.' && attribute[j+k]!='"'; k++){
                id[k] = attribute[j+k];
            }
            id[k]='\0';
            // hash the id
            if(group != NULL){
                free(group);
            }
            group= (char*)malloc(sizeof(char)*1000);
  	        snprintf(group, 1000, "%d %d", block, item_id);
            fprintf(fp_hash_key, "%s\n", id);
            fprintf(fp_hash_val, "%s\n", group);
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
            for(k=0; attribute[j+k]!='.' && attribute[j+k]!='"'; k++){
                id[k] = attribute[j+k];
            }
            id[k]='\0';
            // hash the id
            if(group != NULL){
                free(group);
            }
            group= (char*)malloc(sizeof(char)*1000);
  	        snprintf(group, 1000, "%d %d", block, item_id);
            fprintf(fp_hash_key, "%s\n", id);
            fprintf(fp_hash_val, "%s\n", group);
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
            for(k=0; attribute[j+k]!='.' && attribute[j+k]!='"'; k++){
                id[k] = attribute[j+k];
            }
            id[k]='\0';
            // hash the id
            if(group != NULL){
                free(group);
            }
            group= (char*)malloc(sizeof(char)*1000);
  	        snprintf(group, 1000, "%d %d", block, item_id);
            fprintf(fp_hash_key, "%s\n", id);
            fprintf(fp_hash_val, "%s\n", group);
        }


        //modify the start for better compression
        if(strcmp(feature, gene)  == 0){
            //store the gene start for later uses
            prev_gene = atoi(start);
            prev_trans = prev_gene;
            //write to the new start
            if(prev_gene_end == -1){
                sprintf(new_start,"%d", atoi(start));               
            }
            else{
                sprintf(new_start,"%d", atoi(start) - prev_gene_end);                 
            }
            prev_gene_end = atoi(end);
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
                prev_exon = atoi(end);
            }
            else{
                sprintf(new_start,"%d", prev_exon - atoi(end));
                prev_exon = atoi(start);
            }

            new_transcript = 0;
        }
        else{
            sprintf(new_start,"%d", atoi(start) - prev_trans);
        }
        //output the new start
        fprintf(fp_start, "%s\n", new_start);
    }

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
    return 0;
}

int gff3_compressor(FILE* fp, int length, int* chr_table, int* block_min_table,int* block_max_table, int block_size){
    int i, j, k, m;
    char empty = ' ';
    char prev_chr[100];
    int new_transcript =0, new_block =1;
    int len;
    int block, gene_numbers=0, item_id;
    int prev_gene_end = -1;
    int position_min, position_max;
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
    char* group = NULL;
    char* id;
    id =(char*) malloc(100*sizeof(char));

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
    for(i=0; i<7; i++){
        fgets(comments, BUFFSIZE, fp);
    }
    item_id = -1;
    int first_chr = 1;
    int c = 0;
    chr_table[c] = 0;
    for(i=0; i< length; i++){
        item_id++;
        //read in the seqname
        fscanf(fp, "%s", seqname);
        if(first_chr == 1){
            first_chr = 0;
            strcpy(prev_chr, seqname);           
        }

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
                prev_gene_end = -1;
                gene_numbers=0;
                item_id=0;
                //update the min and max table
                block_min_table[block-1]= position_min;
                block_max_table[block-1]= position_max;
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
        //check if we need to update the min or max
        if(gene_numbers==0 || atoi(start) < position_min){
            position_min= atoi(start);
        }
        if(gene_numbers==0 || atoi(start) > position_max){
            position_max= atoi(start);
        }
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
        if((strcmp(feature, gene)  == 0)||(strcmp(feature, transcript) == 0)||(strcmp(feature, exon) == 0)){
            //find the starting position of gene id
            memset(refer,0,sizeof(refer));
            m=0;
            for(j= 0; j< len; j++){
                if(attribute[j]=='E'){
                    break;
                }
            }
            for(k=0; attribute[j+k]!=';'; k++){
                id[k] = attribute[j+k];
            }
            id[k]='\0';
            // hash the id
            if(group != NULL){
                free(group);
            }

            group= (char*)malloc(sizeof(char)*1000);
            snprintf(group, 1000, "%d %d", block, item_id);
            fprintf(fp_hash_key, "%s\n", id);
            fprintf(fp_hash_val, "%s\n", group);
            //  if(item_id==0){
            //     printf("total group:%s\n", group);
            // }

        }



        //modify the start for better compression
        if(strcmp(feature, gene)  == 0){
            //store the gene start for later uses
            prev_gene = atoi(start);
            prev_trans = prev_gene;
            //write to the new start
            if(prev_gene_end == -1){
                sprintf(new_start,"%d", atoi(start));               
            }
            else{
                sprintf(new_start,"%d", atoi(start) - prev_gene_end);                 
            }
            prev_gene_end = atoi(end);
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
                prev_exon = atoi(end);
            }
            else{
                sprintf(new_start,"%d", prev_exon - atoi(end));
                prev_exon = atoi(start);
            }

            new_transcript = 0;
        }
        else{
            sprintf(new_start,"%d", atoi(start) - prev_trans);
        }
        //output the new start
        fprintf(fp_start, "%s\n", new_start);
    }

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
    return 0;
}

int gtf_compressor2(FILE* fp, int length, int filetype){
	int i;
	char empty = ' ';
	int new_transcript =0;
    //create file pointers for all the output files
	FILE *fp_sq, *fp_src, *fp_fea, *fp_start, *fp_delta, *fp_att, *fp_comments;
	FILE *fp_score, *fp_strand, *fp_frame_cds, *fp_frame_start, *fp_frame_stop;
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

	int prev_gene =0, prev_trans =0, prev_exon =0;
    int prev_gene_end = -1;
	char* gene ="gene";
	char* transcript ="transcript";
	char* exon = "exon";
	char* cds = "CDS";
	char* start_codon = "start_codon";
	char* stop_codon = "stop_codon";

    //create the output files
    //write the seqname to a file
    fp_sq = fopen("GTF_parsed2/gtf_seqname.txt", "w+");
    //write the source to a file
    fp_src = fopen("GTF_parsed2/gtf_source.txt", "w+");
    //write the feature to a file
    fp_fea = fopen("GTF_parsed2/gtf_feature.txt", "w+");
    //write the start to a file
    fp_start = fopen("GTF_parsed2/gtf_start.txt", "w+");
    //write the difference of start and stop to a file
    fp_delta = fopen("GTF_parsed2/gtf_delta.txt", "w+");
    //write the score to a file
    fp_score = fopen("GTF_parsed2/gtf_score.txt", "w+");
    //write the strand to a file
    fp_strand = fopen("GTF_parsed2/gtf_strand.txt", "w+");
    //write the attribute to a file
    fp_att = fopen("GTF_parsed2/gtf_attribute.txt", "w+");
    //write the frames of CDS to a file
    fp_frame_cds = fopen("GTF_parsed2/gtf_frame_cds.txt", "w+");
    //write the frames of start to a file
    fp_frame_start = fopen("GTF_parsed2/gtf_frame_start.txt", "w+");
    //write the frames of stop to a file
    fp_frame_stop = fopen("GTF_parsed2/gtf_frame_stop.txt", "w+");

    //store the first five lines of comments

    for(i=0; i<filetype; i++){
        fgets(comments, BUFFSIZE, fp);
        fprintf(fp_att, "%s", comments);
    }
    int j;
    char ch;
    //extract each column and write them into the output files
    for(i=0; i< length; i++){
    	//read in the seqname
        fscanf(fp, "%s", seqname);
        if(seqname[0]=='#'){
            if(seqname[2] == '#'){
                fprintf(fp_sq, "%s\n", seqname);                
            }
            else{
                fgets(comments, BUFFSIZE, fp);
                strcat(seqname, comments);
                fprintf(fp_sq, "%s", seqname);                  
            }

            continue;            
        }
        //output the seqname
        fprintf(fp_sq, "%s\n", seqname);

        //read in the source

        fscanf(fp, "%s", source);
        if(!strcmp(source, "Curated")){
            strcat(source, " Genomic");
            fscanf(fp, "%s", comments);
        }
        //output the source
        fprintf(fp_src, "%s\n", source);

        //read in the feature
        fscanf(fp, "%s", feature);
        //output the source
        fprintf(fp_fea, "%s\n", feature);

        //read in the start
        fscanf(fp, "%s", start);
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
        fprintf(fp_strand, "%s\n", strand);

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


        //modify the start for better compression
        if(strcmp(feature, gene)  == 0){
    		//store the gene start for later uses
    		prev_gene = atoi(start);
    		prev_trans = prev_gene;
    		//write to the new start
            if(prev_gene_end == -1){
                sprintf(new_start,"%d", atoi(start));               
            }
            else{
                sprintf(new_start,"%d", atoi(start) - prev_gene_end);                 
            }
            prev_gene_end = atoi(end);
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
                prev_exon = atoi(end);
            }
            else{
                sprintf(new_start,"%d", prev_exon - atoi(end));
                prev_exon = atoi(start);
            }

            new_transcript = 0;
    	}
    	else{
    	    sprintf(new_start,"%d", atoi(start) - prev_trans);
    	}
    	//output the new start
        fprintf(fp_start, "%s\n", new_start);
    }
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
    system("BSC/bsc e GTF_parsed2/gtf_seqname.txt GTF_compressed2/gtf_seqname_compressed");
    system("BSC/bsc e GTF_parsed2/gtf_source.txt GTF_compressed2/gtf_source_compressed");
    system("BSC/bsc e GTF_parsed2/gtf_feature.txt GTF_compressed2/gtf_feature_compressed");
    system("BSC/bsc e GTF_parsed2/gtf_start.txt GTF_compressed2/gtf_start_compressed");
    system("BSC/bsc e GTF_parsed2/gtf_delta.txt GTF_compressed2/gtf_delta_compressed");
    system("BSC/bsc e GTF_parsed2/gtf_score.txt GTF_compressed2/gtf_score_compressed");
    system("BSC/bsc e GTF_parsed2/gtf_strand.txt GTF_compressed2/gtf_strand_compressed");
    system("BSC/bsc e GTF_parsed2/gtf_attribute.txt GTF_compressed2/gtf_attribute_compressed");
    system("BSC/bsc e GTF_parsed2/gtf_frame_cds.txt GTF_compressed2/gtf_gtf_frame_cds_compressed");
    system("BSC/bsc e GTF_parsed2/gtf_frame_start.txt GTF_compressed2/gtf_frame_start_compressed");
    system("BSC/bsc e GTF_parsed2/gtf_frame_stop.txt GTF_compressed2/gtf_frame_stop_compressed");

    return 0;
}

int gtf_decompressor(FILE* fp, int length, int filetype){
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
    // char comments[BUFFSIZE];
    char* gene ="gene";
    char* transcript ="transcript";
    char* exon = "exon";
    char* cds = "CDS";
    char* start_codon = "start_codon";
    char* stop_codon = "stop_codon";

    int prev_gene =0, prev_trans =0, prev_exon =0;
    int new_transcript = 0;
    int prev_gene_end = -1;

     //create file pointers for all the files
    FILE *fp_sq, *fp_src, *fp_fea, *fp_start, *fp_delta, *fp_att, *fp_comments;
    FILE *fp_score, *fp_strand, *fp_frame_cds, *fp_frame_start, *fp_frame_stop;
    //open all the files
    fp_sq = fopen("GTF_parsed2/gtf_seqname.txt", "r");
    //open the source file
    fp_src = fopen("GTF_parsed2/gtf_source.txt", "r");
    //open the feature file
    fp_fea = fopen("GTF_parsed2/gtf_feature.txt", "r");
    //open the start file
    fp_start = fopen("GTF_parsed2/gtf_start.txt", "r");
    //open the difference of start and stop file
    fp_delta = fopen("GTF_parsed2/gtf_delta.txt", "r");
    //open the score file
    fp_score = fopen("GTF_parsed2/gtf_score.txt", "r");
    //open the strand file
    fp_strand = fopen("GTF_parsed2/gtf_strand.txt", "r");
    //open the attribute file
    fp_att = fopen("GTF_parsed2/gtf_attribute.txt", "r");
    //open the frames of CDS file
    fp_frame_cds = fopen("GTF_parsed2/gtf_frame_cds.txt", "r");
    //open the frames of start file
    fp_frame_start = fopen("GTF_parsed2/gtf_frame_start.txt", "r");
    //open the frames of stop file
    fp_frame_stop = fopen("GTF_parsed2/gtf_frame_stop.txt", "r");

    //write the comments to the gtf file
    int i;

    fgets(attribute, BUFFSIZE, fp_att);
    while(attribute[0] == '#'){
        fprintf(fp, "%s", attribute);
        fgets(attribute, BUFFSIZE, fp_att);
    }
    //start to combine all columns into the goal file
    for(i=0; i< length; i++){
        //read in all the rest
        fscanf(fp_sq, "%s", seqname);
        if(seqname[0]=='#'){
            if(seqname[2] == '#'){
                fprintf(fp, "%s\n", seqname);                
            }
            else{
                fgets(comments, BUFFSIZE, fp_sq);
                strcat(seqname, comments);
                fprintf(fp, "%s", seqname);                  
            }
            continue;  
        }

        fscanf(fp_src, "%s", source);
        if(!strcmp(source, "Curated")){
            strcat(source, " Genomic");
            fscanf(fp_src, "%s", comments);
        }
        fscanf(fp_fea, "%s", feature);
        fscanf(fp_start, "%s", start);
        fscanf(fp_delta, "%s", delta);
        fscanf(fp_score, "%s", score);
        if(i!=0){
            fgets(attribute, BUFFSIZE, fp_att);
        }


        //output the seqname
        fprintf(fp, "%s\t", seqname);              


        //output the source
        fprintf(fp, "%s\t", source);
        //output the feature
        fprintf(fp, "%s\t", feature);
        //recover the strand
        fscanf(fp_strand, "%s", strand);

        //recover the start
        if(strcmp(feature, gene)  == 0){
            //store the gene start for later uses
            if(prev_gene_end == -1){
                 prev_gene = atoi(start);               
            }
            else{
                prev_gene = atoi(start) + prev_gene_end;
            }

            prev_trans = prev_gene;
            //write to the new start
            sprintf(new_start,"%d", prev_gene);
            prev_gene_end =  atoi(delta) + atoi(new_start);
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
               prev_exon = atoi(new_start) +atoi(delta);
            }
            else{
               sprintf(new_start,"%d", prev_exon - atoi(start) - atoi(delta));
               prev_exon = atoi(new_start);
            }

            new_transcript = 0;
        }
        else{
            sprintf(new_start,"%d", atoi(start) + prev_trans);
        }
        //output the start
        fprintf(fp, "%s\t", new_start);
        //recover the end
        sprintf(end, "%d", atoi(delta) + atoi(new_start));
        //output the end
        fprintf(fp, "%s\t", end);
        //output the score
        fprintf(fp, "%s\t", score);
        //output the strand
        fprintf(fp, "%s\t", strand);
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
        fprintf(fp, "%s\t", frame);

        //output the attribute
        fprintf(fp, "%s", attribute);
    }
    return 0;

}

int expression_compressor(FILE* fp, int length, int block_size){
    int i, j, k, m;
    char prev_chr[100];
    int new_transcript =0, new_block =1;
    int block, gene_numbers=0, item_id, prev_id;

    //create file pointers for all the output files
    FILE *fp_sample, *fp_ets_counts;
    FILE *fp_tpm, *fp_eff_len, *fp_len, *fp_hash_key, *fp_hash_val;
    //array to hold each column
    char target[BUFFSIZE];
    char ID_tmp[BUFFSIZE];
    char ID1[BUFFSIZE];
    char ID2[BUFFSIZE];
    char ID3[BUFFSIZE];
    char ID4[BUFFSIZE];
    char ID5[BUFFSIZE];
    char ID6[BUFFSIZE];
    char sample[BUFFSIZE];
    char ets_counts[BUFFSIZE];
    char tpm[BUFFSIZE];
    char eff_len[BUFFSIZE];
    char len[BUFFSIZE];
    char refer[BUFFSIZE];
    char* group=NULL;
    char* id;

    int prev_gene =0, prev_trans =0, prev_exon =0;

    char temp[200];
    char sample_name[200];
    char ets_counts_name[200];
    char tpm_name[200];
    char eff_len_name[200];
    char len_name[200];
    char compress_cmd[200];

    char* sample_name_prefix= "expression_parsed/expression_sample_" ;
    char* ets_counts_name_prefix= "expression_parsed/expression_ets_counts_";
    char* tpm_name_prefix= "expression_parsed/expression_tpm_";
    char* eff_len_name_prefix= "expression_parsed/expression_eff_len_";
    char* len_name_prefix= "expression_parsed/expression_len_";

    char* compress_cmd_prefix= "BSC/bsc e ";
    char* cmd_suffix_sample= " expression_compressed/expression_compressed_sample_";
    char* cmd_suffix_ets_counts=" expression_compressed/expression_compressed_ets_counts_";
    char* cmd_suffix_tpm=" expression_compressed/expression_compressed_tpm_";
    char* cmd_suffix_eff_len=" expression_compressed/expression_compressed_eff_len_";
    char* cmd_suffix_len=" expression_compressed/expression_compressed_len_";

    char prev_ID1[BUFFSIZE];
    char prev_ID2[BUFFSIZE];
    char prev_ID3[BUFFSIZE];
    char prev_ID4[BUFFSIZE];
    char prev_ID5[BUFFSIZE];
    char prev_ID6[BUFFSIZE];

    char block_number[200];



    //extract each column and write them into the output files
    block= 0;
    sprintf(block_number,"%d", block);
    //create the output files
    strcpy(sample_name, sample_name_prefix);
    strcpy(ets_counts_name, ets_counts_name_prefix);
    strcpy(tpm_name, tpm_name_prefix);
    strcpy(eff_len_name, eff_len_name_prefix);
    strcpy(len_name, len_name_prefix);

    strcat(strcat(sample_name, block_number),".txt");
    strcat(strcat(ets_counts_name, block_number),".txt");
    strcat(strcat(tpm_name, block_number),".txt");
    strcat(strcat(eff_len_name, block_number),".txt");
    strcat(strcat(len_name, block_number),".txt");


    fp_sample = fopen(sample_name, "w+");

    fp_ets_counts = fopen(ets_counts_name, "w+");

    fp_tpm = fopen(tpm_name, "w+");

    fp_eff_len = fopen(eff_len_name, "w+");

    fp_len = fopen(len_name, "w+");

    fp_hash_key= fopen("index_tables/expression_key.txt", "w+");
    fp_hash_val= fopen("index_tables/expression_value.txt", "w+");
    //store the first line of comments
    for(i=0; i<1; i++){
        fgets(temp, BUFFSIZE, fp);
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
        //extract the ID or IDs in all specified databases
        //get rid of the '|' first
        //also get rid of the '.'
        char *s, *t;
        s= strtok(target, "|");
        strcpy(ID1, s);
        s= strtok(NULL, "|");
        strcpy(ID2, s);
        s= strtok(NULL, "|");
        strcpy(ID3, s);
        s= strtok(NULL, "|");
        strcpy(ID4, s);
        s= strtok(NULL, "|");
        strcpy(ID5, s);
        s= strtok(NULL, "|");
        strcpy(ID6, s);
        //also get rid of the '.'
        t= strtok(ID1, ".");
        strcpy(ID1, t);
        t= strtok(ID2, ".");
        strcpy(ID2, t);
        t= strtok(ID3, ".");
        strcpy(ID3, t);
        t= strtok(ID4, ".");
        strcpy(ID4, t);

        //check if we need to update the block
        if((i!=0) && (strcmp(ID1, prev_ID1))){
            gene_numbers++;
            if(gene_numbers == block_size){
                block++;
                new_block= 1;
                gene_numbers=0;
                item_id=0;
                fclose(fp_sample);
                fclose(fp_ets_counts);
                fclose(fp_tpm);
                fclose(fp_eff_len);
                fclose(fp_len);

                //compress the files using bsc algorithm
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, sample_name), cmd_suffix_sample);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
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
                strcpy(sample_name, sample_name_prefix);
                strcpy(ets_counts_name, ets_counts_name_prefix);
                strcpy(tpm_name, tpm_name_prefix);
                strcpy(eff_len_name, eff_len_name_prefix);
                strcpy(len_name, len_name_prefix);

                strcat(strcat(sample_name, block_number),".txt");
                strcat(strcat(ets_counts_name, block_number),".txt");
                strcat(strcat(tpm_name, block_number),".txt");
                strcat(strcat(eff_len_name, block_number),".txt");
                strcat(strcat(len_name, block_number),".txt");


                fp_sample = fopen(sample_name, "w+");

                fp_ets_counts = fopen(ets_counts_name, "w+");

                fp_tpm = fopen(tpm_name, "w+");

                fp_eff_len = fopen(eff_len_name, "w+");

                fp_len = fopen(len_name, "w+");

            }
            // hash the id
            if(group!=NULL){
                free(group);
            }

            group= (char*)malloc(sizeof(char)*1000);

  	        snprintf(group, 1000, "%d %d %d", block, prev_id, item_id - 1);


            //write IDs from all the databases in the index table files
            fprintf(fp_hash_key, "%s\n", prev_ID1);
            fprintf(fp_hash_val, "%s\n", group);
            fprintf(fp_hash_key, "%s\n", prev_ID2);
            fprintf(fp_hash_val, "%s\n", group);
            fprintf(fp_hash_key, "%s\n", prev_ID3);
            fprintf(fp_hash_val, "%s\n", group);
            fprintf(fp_hash_key, "%s\n", prev_ID4);
            fprintf(fp_hash_val, "%s\n", group);
            fprintf(fp_hash_key, "%s\n", prev_ID5);
            fprintf(fp_hash_val, "%s\n", group);
            fprintf(fp_hash_key, "%s\n", prev_ID6);
            fprintf(fp_hash_val, "%s\n", group);

            prev_id= item_id;
        }

        fprintf(fp_sample, "%s\n", sample);
        fprintf(fp_ets_counts, "%s\n", ets_counts);
        fprintf(fp_tpm, "%s\n",tpm);
        fprintf(fp_eff_len, "%s\n",eff_len);
        fprintf(fp_len, "%s\n",len);

        //read in the last column
        fgets(temp, BUFFSIZE, fp);
        strcpy(prev_ID1, ID1);
        strcpy(prev_ID2, ID2);
        strcpy(prev_ID3, ID3);
        strcpy(prev_ID4, ID4);
        strcpy(prev_ID5, ID5);
        strcpy(prev_ID6, ID6);

    }
    //compress the files using bsc algorithm
    strcpy(compress_cmd, compress_cmd_prefix);
    strcat(strcat(compress_cmd, sample_name), cmd_suffix_sample);
    strcat(compress_cmd, block_number);
    system(compress_cmd);
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
    //close all the files
    fclose(fp_sample);
    fclose(fp_ets_counts);
    fclose(fp_tpm);
    fclose(fp_eff_len);
    fclose(fp_len);
    return 0;
}

int sparse_compressor(FILE* fp, FILE* fp_gene_id, int length, int block_size){
    int i, j, k, m;
    char prev_chr[100];
    int new_gene =0, new_block =1, prev_row;
    int block, gene_numbers=0, item_id, prev_id;

    //create file pointers for all the output files
    FILE *fp_column, *fp_value;
    FILE *fp_hash_key, *fp_hash_val;
    //array to hold each column
    char row[BUFFSIZE];
    char column[BUFFSIZE];
    char value[BUFFSIZE];
    char refer[BUFFSIZE];
    char gene_identifier[BUFFSIZE];
    char gene_name[BUFFSIZE];
    char* group= NULL;
    char* id;

    int prev_gene =0, prev_trans =0, prev_exon =0;
    int gene_index;
    char temp[200];
    char column_name[200];
    char value_name[200];
    char compress_cmd[200];

    char* column_name_prefix= "sparse_parsed/sparse_column_" ;
    char* value_name_prefix= "sparse_parsed/sparse_value_";

    char* compress_cmd_prefix= "BSC/bsc e ";
    char* cmd_suffix_column= " sparse_compressed/sparse_compressed_column_";
    char* cmd_suffix_value=" sparse_compressed/sparse_compressed_value_";

    char block_number[200];


    //extract each column and write them into the output files
    block= 0;
    sprintf(block_number,"%d", block);
    //create the output files
    strcpy(column_name, column_name_prefix);
    strcpy(value_name, value_name_prefix);

    strcat(strcat(column_name, block_number),".txt");
    strcat(strcat(value_name, block_number),".txt");


    fp_column = fopen(column_name, "w+");

    fp_value = fopen(value_name, "w+");


    fp_hash_key= fopen("index_tables/sparse_key.txt", "w+");
    fp_hash_val= fopen("index_tables/sparse_value.txt", "w+");

    //store the first line of comments
    // for(i=0; i<1; i++){
    //     fgets(temp, BUFFSIZE, fp);
    // }
    item_id = -1;
    int item_id_flag = 0;
    prev_id= 0;
    gene_index= 0;
    for(i=0; i< length; i++){
        item_id++;
        //read in the row
        fscanf(fp, "%s", row);

        //read in the column
        fscanf(fp, "%s", column);
        //read in the value
        fscanf(fp, "%s", value);

        //check if we need to update the block
        if((i!=0) && (atoi(row)!= prev_row)){
            gene_numbers++;
            if(gene_numbers == block_size){
                block++;
                new_block= 1;
                gene_numbers=0;
                item_id_flag=1;
                fclose(fp_column);
                fclose(fp_value);

                //compress the files using bsc algorithm
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, column_name), cmd_suffix_column);
                strcat(compress_cmd, block_number);
                system(compress_cmd);
                strcpy(compress_cmd, compress_cmd_prefix);
                strcat(strcat(compress_cmd, value_name), cmd_suffix_value);
                strcat(compress_cmd, block_number);
                system(compress_cmd);

                sprintf(block_number,"%d", block);
                //create the output files
                strcpy(column_name, column_name_prefix);
                strcpy(value_name, value_name_prefix);

                strcat(strcat(column_name, block_number),".txt");
                strcat(strcat(value_name, block_number),".txt");


                fp_column = fopen(column_name, "w+");
                fp_value = fopen(value_name, "w+");

            }
            //extract the id
            while(gene_index != prev_row){
                fscanf(fp_gene_id, "%s", gene_identifier);
                fscanf(fp_gene_id, "%s", gene_name);
                fgets(temp, BUFFSIZE, fp_gene_id);
                gene_index++;
            }

            // hash the id
            if(group == NULL){
                free(group);
            }
            group= (char*)malloc(sizeof(char)*1000);
  	        snprintf(group, 1000, "%d %d %d", block, prev_id, item_id - 1);
            fprintf(fp_hash_key, "%s\n", gene_identifier);
            fprintf(fp_hash_val, "%s\n", group);
            fprintf(fp_hash_key, "%s\n", gene_name);
            fprintf(fp_hash_val, "%s\n", group);

            if(item_id_flag == 1){
                item_id= 0;
                item_id_flag=0;
            }
            prev_id= item_id;
        }

        fprintf(fp_column, "%s\n", column);
        fprintf(fp_value, "%s\n", value);

        //read in the last column
        fgets(temp, BUFFSIZE, fp);
        prev_row = atoi(row);
    }
    //compress the files using bsc algorithm
    strcpy(compress_cmd, compress_cmd_prefix);
    strcat(strcat(compress_cmd, column_name), cmd_suffix_column);
    strcat(compress_cmd, block_number);
    system(compress_cmd);
    strcpy(compress_cmd, compress_cmd_prefix);
    strcat(strcat(compress_cmd, value_name), cmd_suffix_value);
    strcat(compress_cmd, block_number);
    system(compress_cmd);
    //close all the files
    fclose(fp_column);
    fclose(fp_value);
    return 0;
}
