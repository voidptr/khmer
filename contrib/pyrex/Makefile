all: bleu-build bleu.so

all_all: clean all tests

tests: all
	nosetests

clean:
	# remove bleu.pyx intermediate files.
	rm -fr *.so bleu.cpp	
	# reset C++ source build
	cd ../bleu && make clean

bleu.so: bleu.cpp
	python setup.py build_ext -i

bleu.cpp: bleu.pyx

bleu-build:
	cd ../bleu && make all
