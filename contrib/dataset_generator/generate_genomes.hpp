/*
 *  generate_genomes.hpp
 *  dataset_generator
 *
 *  Created by Rosangela Canino-Koning on 2/4/11.
 *  Copyright 2011 VoidPtr. All rights reserved.
 *
 */

#include <iostream>
#include <vector>

namespace datagen {
  
  using namespace std;
  
  class GenerateGenomes
  {
  public:
    vector< string > mGenomes;
    
    string lNucleotides[4];
    
    int mNumGenomes;
    int mMaxSize;
    int mMinSize;
    int mStartNumber;

  public:
    
    GenerateGenomes(int aNumGenomes, int aMaxSize, int aMinSize, int aStartNumber=0)
    { 
      lNucleotides[0] = 'a';
      lNucleotides[1] = 'c';
      lNucleotides[2] = 'g';
      lNucleotides[3] = 't';
    
      mNumGenomes = aNumGenomes;
      mMaxSize = aMaxSize;
      mMinSize = aMinSize;
      
      mStartNumber = aStartNumber;
    }
    
    void generate()
    {
      for (int i = 0; i < mNumGenomes; ++i)
      {
        int lGenomeSize = mMinSize + (rand() % (mMaxSize - mMinSize));
        mGenomes.push_back("");
        
        
        for (int j = 0; j < lGenomeSize; ++j)
        {
          mGenomes[i].append( lNucleotides[ rand() % 4 ] );
        } 
      }
    }
  };
}
