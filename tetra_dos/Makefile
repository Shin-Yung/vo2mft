CC=gcc
CFLAGS=-Wall -O3 -DHAVE_INLINE
LDFLAGS=-lgsl -lgslcblas -lm

all: elHamiltonian.o

clean:
	rm *.o *.out

elHamiltonian.o: elHamiltonian.c elHamiltonian.h environment.h
	$(CC) $(CFLAGS) -c elHamiltonian.c -o elHamiltonian.o
