CC = mpicxx
MPCC = CC
OPENMP = -openmp
CFLAGS = -g -std=c++11
LIBS =


TARGETS = mpi serial autograder

all:	$(TARGETS)

autograder: autograder.o common.o matrixCells.o
	$(CC) -o $@ $(LIBS) autograder.o common.o matrixCells.o
serial: serial.o common.o matrixCells.o
	$(CC) -o $@ $(LIBS) serial.o common.o matrixCells.o
mpi: mpi.o common.o matrixCells.o
	$(CC) -o $@ $(LIBS) mpi.o common.o matrixCells.o


autograder.o: autograder.cpp common.h matrixCells.h
	$(CC) -c $(CFLAGS) autograder.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
mpi.o: mpi.cpp common.h
	$(CC) -c $(CFLAGS) mpi.cpp
matrixCells.o: common.h matrixCells.h matrixCells.cpp
	$(CC) -c $(CFLAGS) matrixCells.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt

