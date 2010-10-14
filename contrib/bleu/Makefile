CXXFLAGS=-fPIC -O3 -Wall -pg -I$(KMERDIR)

KMERDIR=../../lib
EXTLIBDIR=../external_lib

all: bleu

clean:
	cd ../../ && make clean
	cd ../external_lib && make clean
	rm -f bleu

bleu: 
	cd ../../ && make lib_files
	cd ../external_lib && make 
	make bleuexe

bleuexe: bleu.cpp $(KMERDIR)/hashtable.o $(KMERDIR)/ktable.o $(KMERDIR)/parsers.o $(EXTLIBDIR)/cBitArray.o
	$(CXX) -o bleu bleu.cpp $(KMERDIR)/hashtable.o $(KMERDIR)/ktable.o $(KMERDIR)/parsers.o $(EXTLIBDIR)/cBitArray.o

#ktable.o: $(KMERDIR)/ktable.cc $(KMERDIR)/ktable.hh

#hashtable.o: $(KMERDIR)/hashtable.cc $(KMERDIR)/hashtable.hh
