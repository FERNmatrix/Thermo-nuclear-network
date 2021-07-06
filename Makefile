CC=g++
CFLAGS=-I.

#Compiler flags, inclde lgsl and lgslcblas for gsl libs
LIBS=-lgsl -lgslcblas
all: EMATS.cpp NATS.cpp
	$(CC) -o EMATS EMATS.cpp $(LIBS)
	$(CC) -o NATS NATS.cpp $(LIBS)

.PHONY: clean

clean:
	rm -f EMATS *.o
	rm -f NATS *.o
