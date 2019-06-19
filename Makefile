
output: main.o compressor.o
	gcc main.o compressor.o -o compressor

main.o: main.c
	gcc -c main.c

compressor.o: compressor.c
	gcc -c compressor.c

clean:
	rm *.o compressor
