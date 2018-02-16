CC = CC
MPCC = CC
OPENMP = -openmp
CFLAGS = -O3
LIBS =


TARGETS = serial openmp autograder

all:	$(TARGETS)


serial: serial.o common.o matrixCells.o
	$(CC) -o $@ $(LIBS) $(OPENMP) serial.o common.o matrixCells.o
autograder: autograder.o common.o matrixCells.o
	$(CC) -o $@ $(LIBS) autograder.o common.o matrixCells.o
openmp: openmp.o common.o matrixCells.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o matrixCells.o



autograder.o: autograder.cpp common.h matrixCells.cpp
	$(CC) -c $(CFLAGS) autograder.cpp matrixCells.cpp 
openmp.o: openmp.cpp common.h matrixCells.cpp 
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h matrixCells.cpp
	$(CC) -c $(OPENMP) $(CFLAGS) serial.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt

