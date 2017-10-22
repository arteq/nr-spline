CXX=g++
FLAGS=-g
LFLAGS=-lm

all: spline.o nrutil.o
	$(CXX) ${FLAGS} ${LFLAGS} -o spline spline.o nrutil.o

spline.o:
	$(CXX) ${FLAGS} -c spline.cpp

nrutil.o:
	$(CXX) ${FLAGS} -c nrutil.c

clean:
	rm -f spline spline.o nrutil.o
	rm -f spline.dat max.dat min.dat mod.dat
