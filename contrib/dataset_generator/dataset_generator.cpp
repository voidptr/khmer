// dataset_generator.cpp : main project file.

#include <iostream>
#include <time.h>
#include <math.h>

#include "generate_genomes.hpp"
#include "sample_genomes.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 6)
  {
    cout << argv[0] << "<number of genomes> <genome max size> <genomes min size> <output read length> <output size (#reads)> [starting genome #]" << endl;
    return 1;
  }
//  cout << "Executable: " << argv[0] << endl;
//  cout << "Parameters: # genomes=" << argv[1] 
//       << ", genome max size=" << argv[2] 
//       << ", genome min size=" << argv[3] 
//       << ", output reads length=" << argv[4]
//       << ", output reads #=" << argv[5];
//  cout << endl;

  

  int lStart = 0;
  if ( argc > 6 )
    lStart = atoi(argv[6]);

  unsigned long lSeed = time(NULL) * clock() + pow(lStart, 3);
  cerr << ":" << lSeed << endl;
  srand(lSeed);
  datagen::GenerateGenomes gg( atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), lStart );
  gg.generate();

  datagen::SampleGenomes sg( gg, atoi(argv[4]), atoi(argv[5]) );
  sg.sample_output();

//  std::cout << "DONE" << std::endl;

  return 0;
}



