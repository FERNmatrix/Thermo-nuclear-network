CC=g++
CFLAGS=-I. -g3

#compiler flags
LIBS=-lgsl -lgslcblas -lgmp
all: EMATS.cpp NEM.cpp #NATS.cpp  EMATS.cpp CEM.cpp #dtrange.cpp CATS.cpp  explicitMatrix.cpp
	$(CC) -o EMATS EMATS.cpp $(CFLAGS) $(LIBS)
	$(CC) -o NEM NEM.cpp $(CFLAGS) $(LIBS)
#	$(CC) -o NATS NATS.cpp $(CFLAGS) $(LIBS)
#	$(CC) -o CEM CEM.cpp $(CFLAGS) $(LIBS)
#	$(CC) -o explicitMatrix explicitMatrix.cpp  $(CFLAGS) $(LIBS)
#	$(CC) -o dtrange dtrange.cpp $(CFLAGS) $(LIBS)
#	$(CC) -o CATS CATS.cpp $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f EMATS *.o
	rm -f NEM *.o
	rm -f NATS *.o
	rm -f explicitMatrix *.o
	rm -f dtrange *.o
	rm -f CATS *.o
	rm -f CEM *.o
