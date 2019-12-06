
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compressor.h"
#include "hash.h"
#define BUFFSIZE 1000

int readTaggedLine(char* filename, int line, char* result, int gene_start, int gene_end)
{
    char tmp[1500];
    FILE *fp, *fp_gene;
    int i=0;
    fp = fopen(filename, "r");
    fp_gene = fopen("gene_block.gtf", "w+");
    if(fp == NULL || fp_gene==NULL) return -1;
    while(fgets(tmp, 1024, fp))
    {
        if(i>gene_end){
            fclose(fp);
            fclose(fp_gene);
            return 0;
        }
        if(i>=gene_start){
            fprintf(fp_gene, "%s", tmp);
        }
        if(i==line)
        {
            strcpy(result, tmp);
        }
        i++;
    }
    fclose(fp);
    fclose(fp_gene);
    return -1;
}

int add_database_id(char* old_id, char* new_id){
    FILE* fp_key= fopen("index_tables/data_key.txt", "r+");
    FILE* fp_value= fopen("index_tables/data_value.txt", "r+");
    char key[256];
    char* value= (char*)malloc(sizeof(char)*256);
    char tmp[256];
    while(fscanf(fp_key, "%s", key)!= EOF){
        fgets(tmp, 256, fp_value);
        if(!strcmp(old_id, key)){
            strcpy(value, tmp);
        }
    }
    fprintf(fp_key, "%s\n", new_id );
    printf("%s\n", value);
    fprintf(fp_value, "%s", value);
    free(value);
    fclose(fp_key);
    fclose(fp_value);
    return 0;
}

char* item_search(int block, int block_id){
    char* item_info;
    item_info = (char*)malloc(sizeof(char)*1000);
    FILE* fp = fopen("info.gtf", "w+");
    char cmd[500];
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
    char *sq_name_cpy= "BSC/bsc d GTF_compressed/compressed_seqname_";
    char *src_name_cpy="BSC/bsc d GTF_compressed/compressed_source_";
    char *fea_name_cpy="BSC/bsc d GTF_compressed/compressed_feature_";
    char *start_name_cpy="BSC/bsc d GTF_compressed/compressed_start_";
    char *delta_name_cpy="BSC/bsc d GTF_compressed/compressed_delta_";
    char *score_name_cpy="BSC/bsc d GTF_compressed/compressed_score_";
    char *strand_name_cpy="BSC/bsc d GTF_compressed/compressed_strand_";
    char *att_name_cpy="BSC/bsc d GTF_compressed/compressed_attribure_";
    char *frame_cds_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_cds_";
    char *frame_start_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_start_";
    char *frame_stop_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_stop_";
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
    char block_number[50];
    char* cmd_prefix= "BSC/bsc d ";
    char* cmd_suffix_sq= "GTF_parsed/decompressed_seqname.txt";
    char* cmd_suffix_src="GTF_parsed/decompressed_source.txt";
    char* cmd_suffix_fea="GTF_parsed/decompressed_feature.txt";
    char* cmd_suffix_start="GTF_parsed/decompressed_start.txt";
    char* cmd_suffix_delta="GTF_parsed/decompressed_delta.txt";
    char* cmd_suffix_score="GTF_parsed/decompressed_score.txt";
    char* cmd_suffix_strand="GTF_parsed/decompressed_strand.txt";
    char* cmd_suffix_att="GTF_parsed/decompressed_attribure.txt";
    char* cmd_suffix_frame_cds="GTF_parsed/decompressed_frame_cds.txt";
    char* cmd_suffix_frame_start="GTF_parsed/decompressed_frame_start.txt";
    char* cmd_suffix_frame_stop="GTF_parsed/decompressed_frame_stop.txt";
    char* gene ="gene";
    char* transcript ="transcript";
    char* exon = "exon";
    char* cds = "CDS";
    char* start_codon = "start_codon";
    char* stop_codon = "stop_codon";

    FILE *fp_sq, *fp_src, *fp_fea, *fp_start, *fp_delta, *fp_att, *fp_comments;
    FILE *fp_score, *fp_strand, *fp_frame_cds, *fp_frame_start, *fp_frame_stop;

    int prev_gene =0, prev_trans =0, prev_exon =0;
    int new_transcript = 0, prev_gene_end=-1;
    //decompress the block first
    // system("tar -xvf results.tar");
    //recover each file

    sprintf(block_number,"%d", block);
    strcpy(sq_name, sq_name_cpy);
    strcat(strcat(sq_name, block_number), " ");
    strcat(sq_name, cmd_suffix_sq);
    system(sq_name);

    sprintf(block_number,"%d", block);
    strcpy(src_name,src_name_cpy);
    strcat(strcat(src_name, block_number), " ");
    strcat(src_name, cmd_suffix_src);
    system(src_name);

    sprintf(block_number,"%d", block);
    strcpy(fea_name, fea_name_cpy);
    strcat(strcat(fea_name, block_number), " ");
    strcat(fea_name, cmd_suffix_fea);
    system(fea_name);

    sprintf(block_number,"%d", block);
    strcpy(start_name, start_name_cpy);
    strcat(strcat(start_name, block_number), " ");
    strcat(start_name, cmd_suffix_start);
    system(start_name);

    sprintf(block_number,"%d", block);
    strcpy(delta_name, delta_name_cpy);
    strcat(strcat(delta_name, block_number), " ");
    strcat(delta_name, cmd_suffix_delta);
    system(delta_name);

    sprintf(block_number,"%d", block);
    strcpy(score_name, score_name_cpy);
    strcat(strcat(score_name, block_number), " ");
    strcat(score_name, cmd_suffix_score);
    system(score_name);

    sprintf(block_number,"%d", block);
    strcpy(strand_name, strand_name_cpy);
    strcat(strcat(strand_name, block_number), " ");
    strcat(strand_name, cmd_suffix_strand);
    system(strand_name);

    sprintf(block_number,"%d", block);
    strcpy(att_name, att_name_cpy);
    strcat(strcat(att_name, block_number), " ");
    strcat(att_name, cmd_suffix_att);
    system(att_name);

    sprintf(block_number,"%d", block);
    strcpy(frame_cds_name, frame_cds_name_cpy);
    strcat(strcat(frame_cds_name, block_number), " ");
    strcat(frame_cds_name, cmd_suffix_frame_cds);
    system(frame_cds_name);

    sprintf(block_number,"%d", block);
    strcpy(frame_start_name, frame_start_name_cpy);
    strcat(strcat(frame_start_name, block_number), " ");
    strcat(frame_start_name, cmd_suffix_frame_start);
    system(frame_start_name);

    sprintf(block_number,"%d", block);
    strcpy(frame_stop_name, frame_stop_name_cpy);
    strcat(strcat(frame_stop_name, block_number), " ");
    strcat(frame_stop_name, cmd_suffix_frame_stop);
    system(frame_stop_name);


    //open all the files
    fp_sq = fopen(cmd_suffix_sq, "r");
    //open the source file
    fp_src = fopen(cmd_suffix_src, "r");
    //open the feature file
    fp_fea = fopen(cmd_suffix_fea, "r");
    //open the start file
    fp_start = fopen(cmd_suffix_start, "r");
    //open the difference of start and stop file
    fp_delta = fopen(cmd_suffix_delta, "r");
    //open the score file
    fp_score = fopen(cmd_suffix_score, "r");
    //open the strand file
    fp_strand = fopen(cmd_suffix_strand, "r");
    //open the attribute file
    fp_att = fopen(cmd_suffix_att, "r");
    //open the frames of CDS file
    fp_frame_cds = fopen(cmd_suffix_frame_cds, "r");
    //open the frames of start file
    fp_frame_start = fopen(cmd_suffix_frame_start, "r");
    //open the frames of stop file
    fp_frame_stop = fopen(cmd_suffix_frame_stop, "r");
    //start to combine all columns into the goal file
    //the start and end lines for current gene
    int gene_start;
    int gene_end;
    int line_index = 0;
    int reach_goal = 0;
    while(fscanf(fp_sq, "%s", seqname)!=EOF){
        //read in all the rest
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
            if(reach_goal == 1){
                gene_end=line_index - 1;
                break;
            }
            gene_start=line_index;
            fscanf(fp_strand, "%s", strand);
        }
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
        //check if it reaches the desired line
        if(line_index == block_id){
            reach_goal = 1;
        }
        line_index++;
    }
    fclose(fp);
    readTaggedLine("info.gtf", block_id, item_info, gene_start, gene_end);
    return item_info;

}

int rangeSearch(int start_pos, int end_pos, int chr, int* chr_table, int* min_table, int* max_table){
    char cmd[500];
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
    char *sq_name_cpy= "BSC/bsc d GTF_compressed/compressed_seqname_";
    char *src_name_cpy="BSC/bsc d GTF_compressed/compressed_source_";
    char *fea_name_cpy="BSC/bsc d GTF_compressed/compressed_feature_";
    char *start_name_cpy="BSC/bsc d GTF_compressed/compressed_start_";
    char *delta_name_cpy="BSC/bsc d GTF_compressed/compressed_delta_";
    char *score_name_cpy="BSC/bsc d GTF_compressed/compressed_score_";
    char *strand_name_cpy="BSC/bsc d GTF_compressed/compressed_strand_";
    char *att_name_cpy="BSC/bsc d GTF_compressed/compressed_attribure_";
    char *frame_cds_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_cds_";
    char *frame_start_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_start_";
    char *frame_stop_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_stop_";
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
    char block_number[50];
    char* cmd_prefix= "BSC/bsc d ";
    char* cmd_suffix_sq= "GTF_parsed/decompressed_seqname.txt";
    char* cmd_suffix_src="GTF_parsed/decompressed_source.txt";
    char* cmd_suffix_fea="GTF_parsed/decompressed_feature.txt";
    char* cmd_suffix_start="GTF_parsed/decompressed_start.txt";
    char* cmd_suffix_delta="GTF_parsed/decompressed_delta.txt";
    char* cmd_suffix_score="GTF_parsed/decompressed_score.txt";
    char* cmd_suffix_strand="GTF_parsed/decompressed_strand.txt";
    char* cmd_suffix_att="GTF_parsed/decompressed_attribure.txt";
    char* cmd_suffix_frame_cds="GTF_parsed/decompressed_frame_cds.txt";
    char* cmd_suffix_frame_start="GTF_parsed/decompressed_frame_start.txt";
    char* cmd_suffix_frame_stop="GTF_parsed/decompressed_frame_stop.txt";
    char* gene ="gene";
    char* transcript ="transcript";
    char* exon = "exon";
    char* cds = "CDS";
    char* start_codon = "start_codon";
    char* stop_codon = "stop_codon";
    char chr_prefix[20]="chr";
    char chr_str[10];

    sprintf(chr_str,"%d", chr+1);
    strcat(chr_prefix, chr_str);

    // sprintf(chr_prefix,"%d", chr);


    FILE *fp_sq, *fp_src, *fp_fea, *fp_start, *fp_delta, *fp_att, *fp_comments;
    FILE *fp_score, *fp_strand, *fp_frame_cds, *fp_frame_start, *fp_frame_stop;

    int prev_gene =0, prev_trans =0, prev_exon =0, prev_gene_end = -1;
    int new_transcript = 0;
    int i;
    // int start_block = -1;
    //decompress the block first
    // system("tar -xvf results.tar");
    //find the range of the files
    // for(i=chr_table[chr]; i<chr_table[chr+1]; i++){
    //     printf("%d\n", position[i]);
    //     if(start_pos>= position[i] && start_pos<= position[i+1]){
    //      start_block= i;
    //     }

    // }
    // if(start_block < 0){
    //     printf("can't find in this range\n");
    //     return -1;
    // }
    // printf("%d\n", start_block);
    for(i= chr_table[chr]; i<= chr_table[chr+1]; i++){
        //use the index tables to check if the min or max postion
        //of the block fall outside of the searched range
        if(max_table[i] < start_pos || end_pos < min_table[i]){
            //skip the block
            continue;
        }

        //recover each file
        sprintf(block_number,"%d", i);
        strcpy(sq_name, sq_name_cpy);
        strcat(strcat(sq_name, block_number), " ");
        strcat(sq_name, cmd_suffix_sq);
        system(sq_name);

        sprintf(block_number,"%d", i);
        strcpy(src_name,src_name_cpy);
        strcat(strcat(src_name, block_number), " ");
        strcat(src_name, cmd_suffix_src);
        system(src_name);

        sprintf(block_number,"%d", i);
        strcpy(fea_name, fea_name_cpy);
        strcat(strcat(fea_name, block_number), " ");
        strcat(fea_name, cmd_suffix_fea);
        system(fea_name);

        sprintf(block_number,"%d", i);
        strcpy(start_name, start_name_cpy);
        strcat(strcat(start_name, block_number), " ");
        strcat(start_name, cmd_suffix_start);
        system(start_name);

        sprintf(block_number,"%d", i);
        strcpy(delta_name, delta_name_cpy);
        strcat(strcat(delta_name, block_number), " ");
        strcat(delta_name, cmd_suffix_delta);
        system(delta_name);

        sprintf(block_number,"%d", i);
        strcpy(score_name, score_name_cpy);
        strcat(strcat(score_name, block_number), " ");
        strcat(score_name, cmd_suffix_score);
        system(score_name);

        sprintf(block_number,"%d", i);
        strcpy(strand_name, strand_name_cpy);
        strcat(strcat(strand_name, block_number), " ");
        strcat(strand_name, cmd_suffix_strand);
        system(strand_name);

        sprintf(block_number,"%d", i);
        strcpy(att_name, att_name_cpy);
        strcat(strcat(att_name, block_number), " ");
        strcat(att_name, cmd_suffix_att);
        system(att_name);

        sprintf(block_number,"%d", i);
        strcpy(frame_cds_name, frame_cds_name_cpy);
        strcat(strcat(frame_cds_name, block_number), " ");
        strcat(frame_cds_name, cmd_suffix_frame_cds);
        system(frame_cds_name);

        sprintf(block_number,"%d", i);
        strcpy(frame_start_name, frame_start_name_cpy);
        strcat(strcat(frame_start_name, block_number), " ");
        strcat(frame_start_name, cmd_suffix_frame_start);
        system(frame_start_name);

        sprintf(block_number,"%d", i);
        strcpy(frame_stop_name, frame_stop_name_cpy);
        strcat(strcat(frame_stop_name, block_number), " ");
        strcat(frame_stop_name, cmd_suffix_frame_stop);
        system(frame_stop_name);


        //open all the files
        fp_sq = fopen(cmd_suffix_sq, "r");
        //open the source file
        fp_src = fopen(cmd_suffix_src, "r");
        //open the feature file
        fp_fea = fopen(cmd_suffix_fea, "r");
        //open the start file
        fp_start = fopen(cmd_suffix_start, "r");
        //open the difference of start and stop file
        fp_delta = fopen(cmd_suffix_delta, "r");
        //open the score file
        fp_score = fopen(cmd_suffix_score, "r");
        //open the strand file
        fp_strand = fopen(cmd_suffix_strand, "r");
        //open the attribute file
        fp_att = fopen(cmd_suffix_att, "r");
        //open the frames of CDS file
        fp_frame_cds = fopen(cmd_suffix_frame_cds, "r");
        //open the frames of start file
        fp_frame_start = fopen(cmd_suffix_frame_start, "r");
        //open the frames of stop file
        fp_frame_stop = fopen(cmd_suffix_frame_stop, "r");
        prev_gene_end = -1;
       //start to combine all columns into the goal file
        while(fscanf(fp_sq, "%s", seqname)!=EOF){
            //check if the starting position of the block
            //read in all the rest
            fscanf(fp_src, "%s", source);
            fscanf(fp_fea, "%s", feature);
            fscanf(fp_start, "%s", start);
            fscanf(fp_delta, "%s", delta);
            fscanf(fp_score, "%s", score);
            fgets(attribute, BUFFSIZE, fp_att);

            //recover the strand
            if(strcmp(feature, gene)  == 0) {
                fscanf(fp_strand, "%s", strand);
            }
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
            
            //recover the end
            sprintf(end, "%d", atoi(delta) + atoi(new_start));


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
            //check if this should be included or not
            if(atoi(new_start)>= start_pos && atoi(new_start)<= end_pos){
                if(strcmp(chr_prefix, seqname)){
                    continue;
                }
                //output the seqname
                printf("%s    ", seqname);
                //output the source
                printf("%s    ", source);
                //output the feature
                printf("%s    ", feature);
                //output the start
                printf("%s    ", new_start);
                //output the end
                printf("%s    ", end);
                //output the score
                printf("%s    ", score);
                //output the strand
                printf("%s    ", strand);
                //output the frame
                printf("%s    ", frame);
                //output the attribute
                printf("%s", attribute);


            }


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
    }
    return 0;

}

int expressionSearch(int block, int start_id, int end_id){
    char cmd[500];
    char sample[BUFFSIZE];
    char ets[BUFFSIZE];
    char tpm[BUFFSIZE];
    char eff[BUFFSIZE];
    char len[BUFFSIZE];

    char *sample_name_cpy= "BSC/bsc d expression_compressed/expression_compressed_sample_";
    char *ets_name_cpy="BSC/bsc d expression_compressed/expression_compressed_ets_counts_";
    char *tpm_name_cpy="BSC/bsc d expression_compressed/expression_compressed_tpm_";
    char *eff_name_cpy="BSC/bsc d expression_compressed/expression_compressed_eff_len_";
    char *len_name_cpy="BSC/bsc d expression_compressed/expression_compressed_len_";

    char sample_name[100];
    char ets_name[100];
    char tpm_name[100];
    char eff_name[100];
    char len_name[100];

    char block_number[50];
    char* cmd_prefix= "BSC/bsc d ";
    char* cmd_suffix_sample= "expression_parsed/decompressed_sample.txt";
    char* cmd_suffix_ets="expression_parsed/decompressed_ets.txt";
    char* cmd_suffix_tpm="expression_parsed/decompressed_tpm.txt";
    char* cmd_suffix_eff="expression_parsed/decompressed_eff.txt";
    char* cmd_suffix_len="expression_parsed/decompressed_len.txt";


    FILE *fp_sample, *fp_ets, *fp_tpm, *fp_eff, *fp_len;


    sprintf(block_number,"%d", block);
    strcpy(sample_name, sample_name_cpy);
    strcat(strcat(sample_name, block_number), " ");
    strcat(sample_name, cmd_suffix_sample);
    system(sample_name);

    sprintf(block_number,"%d", block);
    strcpy(ets_name,ets_name_cpy);
    strcat(strcat(ets_name, block_number), " ");
    strcat(ets_name, cmd_suffix_ets);
    system(ets_name);

    sprintf(block_number,"%d", block);
    strcpy(tpm_name, tpm_name_cpy);
    strcat(strcat(tpm_name, block_number), " ");
    strcat(tpm_name, cmd_suffix_tpm);
    system(tpm_name);
    sprintf(block_number,"%d", block);
    strcpy(eff_name, eff_name_cpy);
    strcat(strcat(eff_name, block_number), " ");
    strcat(eff_name, cmd_suffix_eff);
    system(eff_name);

    sprintf(block_number,"%d", block);
    strcpy(len_name, len_name_cpy);
    strcat(strcat(len_name, block_number), " ");
    strcat(len_name, cmd_suffix_len);
    system(len_name);



    //open all the files
    fp_sample = fopen(cmd_suffix_sample, "r");
    //open the source file
    fp_ets = fopen(cmd_suffix_ets, "r");
    //open the feature file
    fp_tpm = fopen(cmd_suffix_tpm, "r");
    //open the start file
    fp_eff = fopen(cmd_suffix_eff, "r");
    //open the difference of start and stop file
    fp_len = fopen(cmd_suffix_len, "r");
    int i=0;
    printf("sample    est_count    tpm    eff_len    len\n");
    while(fscanf(fp_sample, "%s", sample)!=EOF){
        //read in all the rest
        fscanf(fp_ets, "%s", ets);
        fscanf(fp_tpm, "%s", tpm);
        fscanf(fp_eff, "%s", eff);
        fscanf(fp_len, "%s", len);
        if(i>= start_id && i<=end_id){
            //output the seqname
            printf("%s    ", sample);
            //output the source
            printf("%s    ", ets);
            //output the feature
            printf("%s    ", tpm);
            //output the start
            printf("%s    ", eff);
            //recover the end
            printf("%s\n", len);
        }
        i++;
    }

    return 0;
}

int sparseSearch(int block, int start_id, int end_id){
    FILE* fp = fopen("output/sparse_search.txt", "w+");
    char cmd[500];
    char column[BUFFSIZE];
    char value[BUFFSIZE];
    char cell_id[BUFFSIZE];
    char temp[200];

    char *column_name_cpy= "BSC/bsc d sparse_compressed/sparse_compressed_column_";
    char *value_name_cpy= "BSC/bsc d sparse_compressed/sparse_compressed_value_";

    char column_name[100];
    char value_name[100];

    char block_number[50];
    char* cmd_prefix= "BSC/bsc d ";
    char* cmd_suffix_column= "sparse_parsed/decompressed_sample.txt";
    char* cmd_suffix_value="sparse_parsed/decompressed_ets.txt";


    FILE *fp_column, *fp_value, *fp_cell;


    sprintf(block_number,"%d", block);
    strcpy(column_name, column_name_cpy);
    strcat(strcat(column_name, block_number), " ");
    strcat(column_name, cmd_suffix_column);
    system(column_name);

    sprintf(block_number,"%d", block);
    strcpy(value_name,value_name_cpy);
    strcat(strcat(value_name, block_number), " ");
    strcat(value_name, cmd_suffix_value);
    system(value_name);

    //open all the files
    fp_column = fopen(cmd_suffix_column, "r");
    //open the source file
    fp_value = fopen(cmd_suffix_value, "r");
    //open the cell id file
    fp_cell = fopen("search_barcodes.tsv", "r");

    int i= 0;
    int col_num= 0;
    printf("column index    cell id    value\n");
    while(fscanf(fp_column, "%s", column)!=EOF){
        //read in all the rest
        fscanf(fp_value, "%s", value);
        if(i>= start_id && i<=end_id){
            //extract the cell id
            while(col_num != atoi(column)){
                fscanf(fp_cell, "%s", cell_id);
                fgets(temp, BUFFSIZE, fp_cell);
                col_num++;
            }
            //output the column
            printf("%s    ", column);
            //output the cell id
            printf("%s    ", cell_id);
            //output the matrix value
            printf("%s\n", value);
        }
        i++;
    }
    fclose(fp);
    fclose(fp_column);
    fclose(fp_value);
    fclose(fp_cell);
    return 0;
}


int rangeSearch_sparse(int start_pos, int end_pos, int chr, int* chr_table, int* min_table, int* max_table){
    char cmd[500];
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
    char comments[BUFFSIZE];
    char *sq_name_cpy= "BSC/bsc d GTF_compressed/compressed_seqname_";
    char *src_name_cpy="BSC/bsc d GTF_compressed/compressed_source_";
    char *fea_name_cpy="BSC/bsc d GTF_compressed/compressed_feature_";
    char *start_name_cpy="BSC/bsc d GTF_compressed/compressed_start_";
    char *delta_name_cpy="BSC/bsc d GTF_compressed/compressed_delta_";
    char *score_name_cpy="BSC/bsc d GTF_compressed/compressed_score_";
    char *strand_name_cpy="BSC/bsc d GTF_compressed/compressed_strand_";
    char *att_name_cpy="BSC/bsc d GTF_compressed/compressed_attribure_";
    char *frame_cds_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_cds_";
    char *frame_start_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_start_";
    char *frame_stop_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_stop_";
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
    char block_number[50];
    char* cmd_prefix= "BSC/bsc d ";
    char* cmd_suffix_sq= "GTF_parsed/decompressed_seqname.txt";
    char* cmd_suffix_src="GTF_parsed/decompressed_source.txt";
    char* cmd_suffix_fea="GTF_parsed/decompressed_feature.txt";
    char* cmd_suffix_start="GTF_parsed/decompressed_start.txt";
    char* cmd_suffix_delta="GTF_parsed/decompressed_delta.txt";
    char* cmd_suffix_score="GTF_parsed/decompressed_score.txt";
    char* cmd_suffix_strand="GTF_parsed/decompressed_strand.txt";
    char* cmd_suffix_att="GTF_parsed/decompressed_attribure.txt";
    char* cmd_suffix_frame_cds="GTF_parsed/decompressed_frame_cds.txt";
    char* cmd_suffix_frame_start="GTF_parsed/decompressed_frame_start.txt";
    char* cmd_suffix_frame_stop="GTF_parsed/decompressed_frame_stop.txt";
    char* gene_id = "gene_id";
    char* transcript_id = "transcript_id";
    char* exon_id = "exon_id";
    char* gene ="gene";
    char* transcript ="transcript";
    char* exon = "exon";
    char* cds = "CDS";
    char* start_codon = "start_codon";
    char* stop_codon = "stop_codon";
    char chr_prefix[20]="chr";
    char chr_str[10];
    char* id;

    sprintf(chr_str,"%d", chr+1);
    strcat(chr_prefix, chr_str);

    // sprintf(chr_prefix,"%d", chr);


    FILE *fp_sq, *fp_src, *fp_fea, *fp_start, *fp_delta, *fp_att, *fp_comments;
    FILE *fp_score, *fp_strand, *fp_frame_cds, *fp_frame_start, *fp_frame_stop;

    int prev_gene =0, prev_trans =0, prev_exon =0;
    int new_transcript = 0, prev_gene_end = -1;
    int i, j, k, m;
    int len;

    //create the hashtable for sparse matrix
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

    for(i= chr_table[chr]; i<= chr_table[chr+1]; i++){
        //use the index tables to check if the min or max postion
        //of the block fall outside of the searched range
        if(max_table[i] < start_pos || end_pos < min_table[i]){
            //skip the block
            continue;
        }

        //recover each file
        sprintf(block_number,"%d", i);
        strcpy(sq_name, sq_name_cpy);
        strcat(strcat(sq_name, block_number), " ");
        strcat(sq_name, cmd_suffix_sq);
        system(sq_name);

        sprintf(block_number,"%d", i);
        strcpy(src_name,src_name_cpy);
        strcat(strcat(src_name, block_number), " ");
        strcat(src_name, cmd_suffix_src);
        system(src_name);

        sprintf(block_number,"%d", i);
        strcpy(fea_name, fea_name_cpy);
        strcat(strcat(fea_name, block_number), " ");
        strcat(fea_name, cmd_suffix_fea);
        system(fea_name);

        sprintf(block_number,"%d", i);
        strcpy(start_name, start_name_cpy);
        strcat(strcat(start_name, block_number), " ");
        strcat(start_name, cmd_suffix_start);
        system(start_name);

        sprintf(block_number,"%d", i);
        strcpy(delta_name, delta_name_cpy);
        strcat(strcat(delta_name, block_number), " ");
        strcat(delta_name, cmd_suffix_delta);
        system(delta_name);

        sprintf(block_number,"%d", i);
        strcpy(score_name, score_name_cpy);
        strcat(strcat(score_name, block_number), " ");
        strcat(score_name, cmd_suffix_score);
        system(score_name);

        sprintf(block_number,"%d", i);
        strcpy(strand_name, strand_name_cpy);
        strcat(strcat(strand_name, block_number), " ");
        strcat(strand_name, cmd_suffix_strand);
        system(strand_name);

        sprintf(block_number,"%d", i);
        strcpy(att_name, att_name_cpy);
        strcat(strcat(att_name, block_number), " ");
        strcat(att_name, cmd_suffix_att);
        system(att_name);

        sprintf(block_number,"%d", i);
        strcpy(frame_cds_name, frame_cds_name_cpy);
        strcat(strcat(frame_cds_name, block_number), " ");
        strcat(frame_cds_name, cmd_suffix_frame_cds);
        system(frame_cds_name);

        sprintf(block_number,"%d", i);
        strcpy(frame_start_name, frame_start_name_cpy);
        strcat(strcat(frame_start_name, block_number), " ");
        strcat(frame_start_name, cmd_suffix_frame_start);
        system(frame_start_name);

        sprintf(block_number,"%d", i);
        strcpy(frame_stop_name, frame_stop_name_cpy);
        strcat(strcat(frame_stop_name, block_number), " ");
        strcat(frame_stop_name, cmd_suffix_frame_stop);
        system(frame_stop_name);


        //open all the files
        fp_sq = fopen(cmd_suffix_sq, "r");
        //open the source file
        fp_src = fopen(cmd_suffix_src, "r");
        //open the feature file
        fp_fea = fopen(cmd_suffix_fea, "r");
        //open the start file
        fp_start = fopen(cmd_suffix_start, "r");
        //open the difference of start and stop file
        fp_delta = fopen(cmd_suffix_delta, "r");
        //open the score file
        fp_score = fopen(cmd_suffix_score, "r");
        //open the strand file
        fp_strand = fopen(cmd_suffix_strand, "r");
        //open the attribute file
        fp_att = fopen(cmd_suffix_att, "r");
        //open the frames of CDS file
        fp_frame_cds = fopen(cmd_suffix_frame_cds, "r");
        //open the frames of start file
        fp_frame_start = fopen(cmd_suffix_frame_start, "r");
        //open the frames of stop file
        fp_frame_stop = fopen(cmd_suffix_frame_stop, "r");
       //start to combine all columns into the goal file
        while(fscanf(fp_sq, "%s", seqname)!=EOF){
            //read in all the rest
            fscanf(fp_src, "%s", source);
            fscanf(fp_fea, "%s", feature);
            fscanf(fp_start, "%s", start);
            fscanf(fp_delta, "%s", delta);
            fscanf(fp_score, "%s", score);
            fgets(attribute, BUFFSIZE, fp_att);

            //recover the strand
            if(strcmp(feature, gene)  == 0) {
                fscanf(fp_strand, "%s", strand);
            }
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
            //recover the end
            sprintf(end, "%d", atoi(delta) + atoi(new_start));


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
            //check if this should be included or not
            if(atoi(new_start)>= start_pos && atoi(new_start)<= end_pos){
                if(strcmp(chr_prefix, seqname)){
                    continue;
                }
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
                    for(k=0; attribute[j+k]!='.'; k++){
                        id[k] = attribute[j+k];
                    }
                    id[k]='\0';

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
                    for(k=0; attribute[j+k]!='.'; k++){
                        id[k] = attribute[j+k];
                    }
                    id[k]='\0';

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
                    for(k=0; attribute[j+k]!='.'; k++){
                        id[k] = attribute[j+k];
                    }
                    id[k]='\0';

                }
                else{
                    continue;
                }

                //print out the information in the sparse matrix

                char* hashval;
                char* s;
                int block;
                int block_id;
                int block_start_id;
                int block_end_id;
                hashval= (char*)malloc(sizeof(char)*100);
                hashval= (char*)ht_get(ht, id);
                if(hashval == NULL){
                    continue;
                }

                s= (char*)malloc(sizeof(char)*200);

                s= strtok(hashval, " ");
                block= atoi(s);

                s= strtok(NULL, " ");
                block_start_id= atoi(s);

                s= strtok(NULL, " ");
                block_end_id= atoi(s);
                //output the information in sparse matrix
                printf("The infomation for gene %s in sparse matrix：\n", id);
                sparseSearch(block, block_start_id, block_end_id);
                //this item's information in GFF file
                //output the seqname
                printf("The infomation for gene %s in GFF file：\n", id);
                printf("%s    ", seqname);
                //output the source
                printf("%s    ", source);
                //output the feature
                printf("%s    ", feature);
                //output the start
                printf("%s    ", new_start);
                //output the end
                printf("%s    ", end);
                //output the score
                printf("%s    ", score);
                //output the strand
                printf("%s    ", strand);
                //output the frame
                printf("%s    ", frame);
                //output the attribute
                printf("%s\n", attribute);


            }

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
    }
    return 0;

}

int rangeSearch_expression(int start_pos, int end_pos, int chr, int* chr_table, int* min_table, int* max_table){
    char cmd[500];
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
    char comments[BUFFSIZE];
    char *sq_name_cpy= "BSC/bsc d GTF_compressed/compressed_seqname_";
    char *src_name_cpy="BSC/bsc d GTF_compressed/compressed_source_";
    char *fea_name_cpy="BSC/bsc d GTF_compressed/compressed_feature_";
    char *start_name_cpy="BSC/bsc d GTF_compressed/compressed_start_";
    char *delta_name_cpy="BSC/bsc d GTF_compressed/compressed_delta_";
    char *score_name_cpy="BSC/bsc d GTF_compressed/compressed_score_";
    char *strand_name_cpy="BSC/bsc d GTF_compressed/compressed_strand_";
    char *att_name_cpy="BSC/bsc d GTF_compressed/compressed_attribure_";
    char *frame_cds_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_cds_";
    char *frame_start_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_start_";
    char *frame_stop_name_cpy="BSC/bsc d GTF_compressed/compressed_frame_stop_";
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
    char block_number[50];
    char* cmd_prefix= "BSC/bsc d ";
    char* cmd_suffix_sq= "GTF_parsed/decompressed_seqname.txt";
    char* cmd_suffix_src="GTF_parsed/decompressed_source.txt";
    char* cmd_suffix_fea="GTF_parsed/decompressed_feature.txt";
    char* cmd_suffix_start="GTF_parsed/decompressed_start.txt";
    char* cmd_suffix_delta="GTF_parsed/decompressed_delta.txt";
    char* cmd_suffix_score="GTF_parsed/decompressed_score.txt";
    char* cmd_suffix_strand="GTF_parsed/decompressed_strand.txt";
    char* cmd_suffix_att="GTF_parsed/decompressed_attribure.txt";
    char* cmd_suffix_frame_cds="GTF_parsed/decompressed_frame_cds.txt";
    char* cmd_suffix_frame_start="GTF_parsed/decompressed_frame_start.txt";
    char* cmd_suffix_frame_stop="GTF_parsed/decompressed_frame_stop.txt";
    char* gene_id = "gene_id";
    char* transcript_id = "transcript_id";
    char* exon_id = "exon_id";
    char* gene ="gene";
    char* transcript ="transcript";
    char* exon = "exon";
    char* cds = "CDS";
    char* start_codon = "start_codon";
    char* stop_codon = "stop_codon";
    char chr_prefix[20]="chr";
    char chr_str[10];
    char* id;

    sprintf(chr_str,"%d", chr+1);
    strcat(chr_prefix, chr_str);

    // sprintf(chr_prefix,"%d", chr);


    FILE *fp_sq, *fp_src, *fp_fea, *fp_start, *fp_delta, *fp_att, *fp_comments;
    FILE *fp_score, *fp_strand, *fp_frame_cds, *fp_frame_start, *fp_frame_stop;

    int prev_gene =0, prev_trans =0, prev_exon =0;
    int new_transcript = 0;
    int i, j, k, m;
    int len;

    //create the hashtable for sparse matrix
    char hash_key[500];
    char* hash_val;
    char temp[100];
    char* retval;

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

    for(i= chr_table[chr]; i<= chr_table[chr+1]; i++){
        //use the index tables to check if the min or max postion
        //of the block fall outside of the searched range
        if(max_table[i] < start_pos || end_pos < min_table[i]){
            //skip the block
            continue;
        }

        //recover each file
        sprintf(block_number,"%d", i);
        strcpy(sq_name, sq_name_cpy);
        strcat(strcat(sq_name, block_number), " ");
        strcat(sq_name, cmd_suffix_sq);
        system(sq_name);

        sprintf(block_number,"%d", i);
        strcpy(src_name,src_name_cpy);
        strcat(strcat(src_name, block_number), " ");
        strcat(src_name, cmd_suffix_src);
        system(src_name);

        sprintf(block_number,"%d", i);
        strcpy(fea_name, fea_name_cpy);
        strcat(strcat(fea_name, block_number), " ");
        strcat(fea_name, cmd_suffix_fea);
        system(fea_name);

        sprintf(block_number,"%d", i);
        strcpy(start_name, start_name_cpy);
        strcat(strcat(start_name, block_number), " ");
        strcat(start_name, cmd_suffix_start);
        system(start_name);

        sprintf(block_number,"%d", i);
        strcpy(delta_name, delta_name_cpy);
        strcat(strcat(delta_name, block_number), " ");
        strcat(delta_name, cmd_suffix_delta);
        system(delta_name);

        sprintf(block_number,"%d", i);
        strcpy(score_name, score_name_cpy);
        strcat(strcat(score_name, block_number), " ");
        strcat(score_name, cmd_suffix_score);
        system(score_name);

        sprintf(block_number,"%d", i);
        strcpy(strand_name, strand_name_cpy);
        strcat(strcat(strand_name, block_number), " ");
        strcat(strand_name, cmd_suffix_strand);
        system(strand_name);

        sprintf(block_number,"%d", i);
        strcpy(att_name, att_name_cpy);
        strcat(strcat(att_name, block_number), " ");
        strcat(att_name, cmd_suffix_att);
        system(att_name);

        sprintf(block_number,"%d", i);
        strcpy(frame_cds_name, frame_cds_name_cpy);
        strcat(strcat(frame_cds_name, block_number), " ");
        strcat(frame_cds_name, cmd_suffix_frame_cds);
        system(frame_cds_name);

        sprintf(block_number,"%d", i);
        strcpy(frame_start_name, frame_start_name_cpy);
        strcat(strcat(frame_start_name, block_number), " ");
        strcat(frame_start_name, cmd_suffix_frame_start);
        system(frame_start_name);

        sprintf(block_number,"%d", i);
        strcpy(frame_stop_name, frame_stop_name_cpy);
        strcat(strcat(frame_stop_name, block_number), " ");
        strcat(frame_stop_name, cmd_suffix_frame_stop);
        system(frame_stop_name);


        //open all the files
        fp_sq = fopen(cmd_suffix_sq, "r");
        //open the source file
        fp_src = fopen(cmd_suffix_src, "r");
        //open the feature file
        fp_fea = fopen(cmd_suffix_fea, "r");
        //open the start file
        fp_start = fopen(cmd_suffix_start, "r");
        //open the difference of start and stop file
        fp_delta = fopen(cmd_suffix_delta, "r");
        //open the score file
        fp_score = fopen(cmd_suffix_score, "r");
        //open the strand file
        fp_strand = fopen(cmd_suffix_strand, "r");
        //open the attribute file
        fp_att = fopen(cmd_suffix_att, "r");
        //open the frames of CDS file
        fp_frame_cds = fopen(cmd_suffix_frame_cds, "r");
        //open the frames of start file
        fp_frame_start = fopen(cmd_suffix_frame_start, "r");
        //open the frames of stop file
        fp_frame_stop = fopen(cmd_suffix_frame_stop, "r");
       //start to combine all columns into the goal file
        while(fscanf(fp_sq, "%s", seqname)!=EOF){
            //read in all the rest
            fscanf(fp_src, "%s", source);
            fscanf(fp_fea, "%s", feature);
            fscanf(fp_start, "%s", start);
            fscanf(fp_delta, "%s", delta);
            fscanf(fp_score, "%s", score);
            fgets(attribute, BUFFSIZE, fp_att);

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
            //recover the end
            sprintf(end, "%d", atoi(delta) + atoi(new_start));


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
            //check if this should be included or not
            if(atoi(new_start)>= start_pos && atoi(new_start)<= end_pos){
                if(strcmp(chr_prefix, seqname)){
                    continue;
                }
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
                    for(k=0; attribute[j+k]!='.'; k++){
                        id[k] = attribute[j+k];
                    }
                    id[k]='\0';

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
                    for(k=0; attribute[j+k]!='.'; k++){
                        id[k] = attribute[j+k];
                    }
                    id[k]='\0';

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
                    for(k=0; attribute[j+k]!='.'; k++){
                        id[k] = attribute[j+k];
                    }
                    id[k]='\0';

                }
                else{
                    continue;
                }

                //print out the information in the sparse matrix
                char* hashval;
                char* s;
                int block;
                int block_id;
                int block_start_id;
                int block_end_id;
                hashval= (char*)malloc(sizeof(char)*100);
                hashval= (char*)ht_get(ht, id);
                if(hashval == NULL){
                    continue;
                }
                s= (char*)malloc(sizeof(char)*50);
                s= strtok(hashval, " ");
                block= atoi(s);
                s= strtok(NULL, " ");
                block_start_id= atoi(s);
                s= strtok(NULL, " ");
                block_end_id= atoi(s);
                //output the information in sparse matrix
                printf("The infomation for item with ID %s in expression file：\n", id);
                expressionSearch(block, block_start_id, block_end_id);
                //this item's information in GFF file
                //output the seqname
                printf("The infomation for item with ID %s in GFF file：\n", id);
                printf("%s    ", seqname);
                //output the source
                printf("%s    ", source);
                //output the feature
                printf("%s    ", feature);
                //output the start
                printf("%s    ", new_start);
                //output the end
                printf("%s    ", end);
                //output the score
                printf("%s    ", score);
                //output the strand
                printf("%s    ", strand);
                //output the frame
                printf("%s    ", frame);
                //output the attribute
                printf("%s\n", attribute);


            }

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
    }
    return 0;

}
