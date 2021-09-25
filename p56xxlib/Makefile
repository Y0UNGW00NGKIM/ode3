ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
CXXFLAGS  += $(ROOTCFLAGS) -I$(P56xxLIBS)/inc -MMD  -Wall -O3
LIBS       = $(ROOTLIBS)
GLIBS      = $(ROOTGLIBS)
GXX	   = g++ $(CXXFLAGS)



all: RKn.o libP56xx.a libP56xx.so

RKn.o: RKn.cpp RKn.hpp
	$(GXX) -c RKn.cpp

libP56xx.a: RKn.o
	ar crs libP56xx.a RKn.o

libP56xx.so: RKn.o
	g++ RKn.o -shared -o libP56xx.so


clean:
	rm -f RKn.o libP56xx.a libP56xx.so
