CC = gcc

erfinv.o: erfinv.c
	$(CC) -O2 -c erfinv.c -o erfinv.o

test: erfinv_test.c erfinv.o
	$(CC) -O0 -ggdb3 erfinv_test.c erfinv.o -o test
	./test

clean:
	rm -rf test erfinv.o
