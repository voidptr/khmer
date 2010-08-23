CXXFLAGS=-g -fPIC -O2 -Wall -pg -I$(KMERDIR)

KMERDIR=../../lib

all: cerulean

clean:
	cd ../../ && make clean
	rm -f cerulean

cerulean: 
	cd ../../ && make lib_files 
	make ceruleanexe

ceruleanexe: cerulean.cpp $(KMERDIR)/hashtable.o $(KMERDIR)/ktable.o $(KMERDIR)/parsers.o
	$(CXX) -o cerulean cerulean.cpp $(KMERDIR)/hashtable.o $(KMERDIR)/ktable.o $(KMERDIR)/parsers.o

#ktable.o: $(KMERDIR)/ktable.cc $(KMERDIR)/ktable.hh

#hashtable.o: $(KMERDIR)/hashtable.cc $(KMERDIR)/hashtable.hh
