CC=g++
CFLAGS=-I.

#LIBS=-lgsl -lgslcbla
LIBS=-lgsl
EMATS: EMATS.cpp
	$(CC) -o EMATS EMATS.cpp $(LIBS)

.PHONY: clean

clean:
	rm -f EMATS *.o
