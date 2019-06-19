#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compressor.h"
#define BUFFSIZE 1000


int gtf_compressor(FILE* fp, int length){
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
	char* gene ="gene";
	char* transcript ="transcript";
	char* exon = "exon";
	char* cds = "CDS";
	char* start_codon = "start_codon";
	char* stop_codon = "stop_codon";

    //create the output files
    //write the seqname to a file
    fp_sq = fopen("compressed/gtf_seqname.txt", "w+");
    //write the source to a file
    fp_src = fopen("compressed/gtf_source.txt", "w+");
    //write the feature to a file
    fp_fea = fopen("compressed/gtf_feature.txt", "w+");
    //write the start to a file
    fp_start = fopen("compressed/gtf_start.txt", "w+");
    //write the difference of start and stop to a file
    fp_delta = fopen("compressed/gtf_delta.txt", "w+");
    //write the score to a file
    fp_score = fopen("compressed/gtf_score.txt", "w+");
    //write the strand to a file
    fp_strand = fopen("compressed/gtf_strand.txt", "w+");
    //write the attribute to a file
    fp_att = fopen("compressed/gtf_attribute.txt", "w+");
    //write the frames of CDS to a file
    fp_frame_cds = fopen("compressed/gtf_frame_cds.txt", "w+");
    //write the frames of start to a file
    fp_frame_start = fopen("compressed/gtf_frame_start.txt", "w+");
    //write the frames of stop to a file
    fp_frame_stop = fopen("compressed/gtf_frame_stop.txt", "w+");

    //store the first five lines of comments
    for(i=0; i<5; i++){
        fgets(comments, BUFFSIZE, fp);
        fprintf(fp_att, "%s", comments);           
    }
    //extract each column and write them into the output files
    for(i=0; i< length; i++){
    	//read in the seqname
        fscanf(fp, "%s", seqname);
        //output the seqname
        fprintf(fp_sq, "%s\n", seqname);

        //read in the source
        fscanf(fp, "%s", source);
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
    fp_sq = fopen("compressed/gtf_seqname.txt", "r");
    //open the source file
    fp_src = fopen("compressed/gtf_source.txt", "r");
    //open the feature file
    fp_fea = fopen("compressed/gtf_feature.txt", "r");
    //open the start file
    fp_start = fopen("compressed/gtf_start.txt", "r");
    //open the difference of start and stop file
    fp_delta = fopen("compressed/gtf_delta.txt", "r");
    //open the score file
    fp_score = fopen("compressed/gtf_score.txt", "r");
    //open the strand file
    fp_strand = fopen("compressed/gtf_strand.txt", "r");
    //open the attribute file
    fp_att = fopen("compressed/gtf_attribute.txt", "r");
    //open the frames of CDS file
    fp_frame_cds = fopen("compressed/gtf_frame_cds.txt", "r");
    //open the frames of start file
    fp_frame_start = fopen("compressed/gtf_frame_start.txt", "r");
    //open the frames of stop file
    fp_frame_stop = fopen("compressed/gtf_frame_stop.txt", "r");

    //write the comments to the gtf file
    for(int i=0; i< 5; i++){
        fgets(comments, BUFFSIZE, fp_att);
        fprintf(fp, "%s", comments);
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