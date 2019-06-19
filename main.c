
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compressor.h"


int main(int argc , char **argv){
    if(argc < 4){
        //no arguments were passed
        printf("not enough arguments are passed!\n");
        return 0;
    } 
    if (strcmp("-c", argv[1]) == 0){
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
        count_lines -=5;
        fp = fopen(argv[3], "r");
        //run the compressor
        gtf_compressor(fp, count_lines);
        if(strcmp("-gzip", argv[2]) == 0){
            //start to compress all the files using gzip
            system("tar -cvzf compressed.tar.gz compressed");           
        }
        else if(strcmp("-bzip2", argv[2]) == 0){
            //start to compress all the files using bzip2
            system("tar -cjvf compressed.tar.bz2 compressed");          
        }
        else if(strcmp("-xz", argv[2]) == 0){
            //start to compress all the files using xz
            system("tar -cJvf compressed.tar.xz compressed");         
        }

    	
        fclose(fp);
    }
    else if(strcmp("-uc", argv[1]) == 0){
        FILE *fp;
        char chr;
        int count_lines = 0;
        if(strcmp("-gzip", argv[2]) == 0){
            //start to decompress all the files using gzip
            system("tar -xvzf compressed.tar.gz");           
        }
        else if(strcmp("-bzip2", argv[2]) == 0){
            //start to decompress all the files using bzip2
            system("tar -xjvf compressed.tar.bz2");          
        }
        else if(strcmp("-xz", argv[2]) == 0){
            //start to decompress all the files using xz
            system("tar -xJvf compressed.tar.xz");         
        }
        fp = fopen("compressed/gtf_seqname.txt", "r");
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

        fp = fopen("result.gtf", "w+");
        gtf_decompressor(fp, count_lines);

    }
	return 0;
}