CC=g++
CFLAGS=-I. -g3

#compiler flags
LIBS=-lgsl -lgslcblas -lgmp -g
all: explicitMatrix.cpp  #EMATS.cpp
#	$(CC) -o EMATS EMATS.cpp $(CFLAGS) $(LIBS)
	$(CC) -o explicitMatrix explicitMatrix.cpp  $(CFLAGS) $(LIBS)
	$(CC) -o EM.out explicitMatrix.cpp $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f EMATS *.o
	rm -f explicitMatrix *.o
