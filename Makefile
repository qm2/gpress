
output: main.o compressor.o hash.o randomaccess.o
	gcc main.o compressor.o hash.o randomaccess.o -o compressor

main.o: main.c
	gcc -c main.c

compressor.o: compressor.c
	gcc -c compressor.c

hash.o: hash.c
	gcc -c hash.c

randomaccess.o: randomaccess.c
	gcc -c randomaccess.c

clean:
	rm *.o compressor

