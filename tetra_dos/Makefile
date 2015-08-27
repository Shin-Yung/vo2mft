CC=gcc
CFLAGS=-Wall -O3 -DHAVE_INLINE
LDFLAGS=-lgsl -lgslcblas -lm
OBJFILES=bstrlib/bstrlib.o bstrlib/bstraux.o elHamiltonian.o DosValues.o ctetra/submesh.o ctetra/dos.o ctetra/numstates.o ctetra/fermi.o ctetra/weights.o ctetra/sum.o ctetra/ecache.o ctetra/tetra.o

all: bstrlib.o bstraux.o elHamiltonian.o DosValues.o RunDosValues.out

clean:
	rm *.o *.out

bstrlib.o:
	$(CC) $(CFLAGS) -c bstrlib/bstrlib.c -o bstrlib/bstrlib.o

bstraux.o:
	$(CC) $(CFLAGS) -c bstrlib/bstraux.c -o bstrlib/bstraux.o

elHamiltonian.o: elHamiltonian.c elHamiltonian.h environment.h
	$(CC) $(CFLAGS) -c elHamiltonian.c -o elHamiltonian.o

DosValues.o: DosValues.c DosValues.h elHamiltonian.o 
	$(CC) $(CFLAGS) -c DosValues.c -o DosValues.o

RunDosValues.out: RunDosValues.c DosValues.o bstrlib/bstrlib.o bstrlib/bstraux.o
	$(CC) $(CFLAGS) -c RunDosValues.c -o RunDosValues.o
	$(CC) RunDosValues.o -o RunDosValues.out $(OBJFILES) $(LDFLAGS)
