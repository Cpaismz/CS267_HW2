CC = CC
MPCC = CC
OPENMP = -openmp
CFLAGS = -O3
LIBS =


TARGETS = mpi serial autograder

all:	$(TARGETS)

autograder: autograder.o common.o matrixCells.o
	$(CC) -o $@ $(LIBS) autograder.o common.o matrixCells.o
serial: serial.o common.o matrixCells.o
	$(CC) -o $@ $(LIBS) serial.o common.o matrixCells.o
mpi: mpi.o common.o matrixCells.o
	$(CC) -o $@ $(LIBS) mpi.o common.o matrixCells.o


autograder.o: autograder.cpp common.h matrixCells.cpp
	$(CC) -c $(CFLAGS) autograder.cpp matrixCells.cpp 
serial.o: serial.cpp common.h matrixCells.cpp 
	$(CC) -c $(CFLAGS) serial.cpp
mpi.o: mpi.cpp common.h matrixCells.cpp 
	$(CC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt

