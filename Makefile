CC = g++ --std=c++11

erfinv.o: erfinv.cpp
	$(CC) -O2 -c erfinv.cpp -o erfinv.o

test: erfinv_test.cpp erfinv.o
	$(CC) -O0 -ggdb3 erfinv_test.cpp erfinv.o -o test
	./test

clean:
	rm -rf test erfinv.o
