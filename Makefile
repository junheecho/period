# Makefile for user_programs

prefix=/usr/local
exec_prefix=/usr/local

CC = gcc -std=c11
CFLAGS = -g -O2
CPP = gcc -E
CPPFLAGS = -I/usr/local/include
CXX = g++ -std=c++11
CXXCPP = g++ -E -std=c++11
CXXFLAGS = -g -O2
LDFLAGS = -Xlinker -rpath -Xlinker /usr/local/lib
LDLIBS =  -L/usr/local/lib -liRRAM -lmpfr -lgmp -lm -lpthread


bin/period-rand:
	g++ --std=c++11 -g -O2 -Wall -I/usr/local/include -Xlinker -rpath -Xlinker /usr/local/lib  src/period-rand.cpp src/random.cpp src/semi-algebraic.cpp  -L/usr/local/lib -liRRAM -lmpfr -lm -lgmp -lpthread -o bin/period-rand

bin/period-rand2:
	g++ --std=c++11 -g -O2 -Wall -I/usr/local/include -Xlinker -rpath -Xlinker /usr/local/lib  src/period-rand2.cpp src/random.cpp src/semi-algebraic.cpp  -L/usr/local/lib -liRRAM -lmpfr -lm -lgmp -lpthread -o bin/period-rand2

bin/period-det:
	g++ --std=c++11 -g -O2 -Wall -I/usr/local/include -Xlinker -rpath -Xlinker /usr/local/lib  src/period-det.cpp src/semi-algebraic.cpp  -L/usr/local/lib -liRRAM -lmpfr -lm -lgmp -lpthread -o bin/period-det

bin/ref-pi:
	g++ --std=c++11 -g -O2 -Wall -I/usr/local/include -Xlinker -rpath -Xlinker /usr/local/lib  src/ref-pi.cpp  -L/usr/local/lib -liRRAM -lmpfr -lm -lgmp -lpthread -o bin/ref-pi

bin/ref-ln2:
	g++ --std=c++11 -g -O2 -Wall -I/usr/local/include -Xlinker -rpath -Xlinker /usr/local/lib  src/ref-ln2.cpp  -L/usr/local/lib -liRRAM -lmpfr -lm -lgmp -lpthread -o bin/ref-ln2

all: bin/period-rand bin/period-rand2 bin/period-det bin/ref-pi bin/ref-ln2

clean:
	rm -rf bin/period-rand bin/period-rand2 bin/period-det bin/ref-pi bin/ref-ln2
