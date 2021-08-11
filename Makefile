CC=g++
CFLAGS=-I. -g3

#compiler flags
LIBS=-lgsl -lgslcblas
all: EMATS.cpp NATS.cpp dtrange.cpp #explicitMatrix.cpp
	$(CC) -o EMATS EMATS.cpp $(CFLAGS) $(LIBS)
	$(CC) -o NATS NATS.cpp $(CFLAGS) $(LIBS)
######	$(CC) -o explicitMatrix explicitMatrix.cpp  $(CFLAGS) $(LIBS)
	$(CC) -o dtrange dtrange.cpp $(CFLAGS) $(LIBS)
#NATS: NATS.cpp
#	$(CC) -o NATS NATS.cpp $(LIBS)

.PHONY: clean

clean:
	rm -f EMATS *.o
	rm -f NATS *.o
	rm -f explicitMatrix *.o
