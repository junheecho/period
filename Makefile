# Makefile for user_programs

prefix=/usr/local
exec_prefix=/usr/local

CC = gcc -std=c11
CFLAGS = -g -O2
CPP = gcc -E
CPPFLAGS =     -I/usr/local/include
CXX = g++ -std=c++11
CXXCPP = g++ -E -std=c++11
CXXFLAGS = -g -O2
LDFLAGS = -Xlinker -rpath -Xlinker /usr/local/lib
LDLIBS =  -L/usr/local/lib -liRRAM -lmpfr -lgmp -lm -lpthread

BIN = period

all: $(BIN)

clean:
	rm -f $(BIN)
