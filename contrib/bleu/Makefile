CXXFLAGS=-g -fPIC -O2 -Wall -pg -I$(KMERDIR)

KMERDIR=../../lib

all: bleu

clean:
	cd ../../ && make clean
	rm -f bleu

bleu: 
	cd ../../ && make lib_files 
	make bleuexe

bleuexe: bleu.cpp $(KMERDIR)/hashtable.o $(KMERDIR)/ktable.o $(KMERDIR)/parsers.o
	$(CXX) -o bleu bleu.cpp $(KMERDIR)/hashtable.o $(KMERDIR)/ktable.o $(KMERDIR)/parsers.o

#ktable.o: $(KMERDIR)/ktable.cc $(KMERDIR)/ktable.hh

#hashtable.o: $(KMERDIR)/hashtable.cc $(KMERDIR)/hashtable.hh
