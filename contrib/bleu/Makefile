###### Bleu Makefile #######
CXX=g++
CXXFLAGS=-O3 -m64 -Wall

KMERDIR=../../lib
EXTLIBDIR=../external_lib

all: bleu

clean:
	cd ../../ && make clean
	cd ../external_lib && make clean
	rm -f bleu

bleu: 
	cd ../../lib/ && make parsers.o
	cd ../external_lib && make 
	make bleuexe
	
tests: 
	cd ../../ && make lib_files
	cd ../external_lib && make 
	make bleutestsexe
	cd ../tests/ && ../bleu/bleutests

bleuexe: bleu.cpp $(KMERDIR)/parsers.o $(EXTLIBDIR)/cBitArray.o
	$(CXX) $(CXXFLAGS) -o bleu bleu.cpp $(KMERDIR)/parsers.o $(EXTLIBDIR)/cBitArray.o

bleutestsexe: bleutests.cpp $(KMERDIR)/parsers.o $(EXTLIBDIR)/cBitArray.o
	$(CXX) $(CXXFLAGS) -o bleu bleu.cpp $(KMERDIR)/parsers.o $(EXTLIBDIR)/cBitArray.o
