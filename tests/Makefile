# DEBUG
CXXFLAGS=-pedantic -Wall -Wsign-compare -std=c++98 -g
LDFLAGS=-g -v

# NODEBUG
#CXXFLAGS=-pedantic -Wall -Wsign-compare -std=c++98 -O3 -DNDEBUG
#LDFLAGS=-v 

IFLAGS=-I/usr/local/include/eigen3 -I../include/ -I../tools/ -I/usr/include/ -I/usr/local/include/
#CXX=g++47


all: test_01


test_01: 
	${CXX} ${CXXFLAGS} ${IFLAGS} -c test_01.cpp
	${CXX} test_01.o -o test_01.a ${LDFLAGS}


clean:
	rm -Rf *.a
	rm -Rf *.o

.PHONY: lexls
