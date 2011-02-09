/*
 *  sample_genomes.hpp
 *  dataset_generator
 *
 *  Created by Rosangela Canino-Koning on 2/4/11.
 *  Copyright 2011 VoidPtr. All rights reserved.
 *
 */

#include <iostream>

namespace datagen {
  
  using namespace std;
  
  class SampleGenomes
  {
  private:
    int mReadLength;
    int mNumReads;
    GenerateGenomes * mGenomes;
    
  public:
    
    SampleGenomes( GenerateGenomes & aGenomes, int aReadLength, int aNumReads)
    { 
      mGenomes = &aGenomes;
      mReadLength = aReadLength;
      mNumReads = aNumReads;
    }
    
    void sample_output()
    {
      for (int i = 0; i < mNumReads; ++i)
      {
        int lGenome = rand() % mGenomes->mNumGenomes;
        int lGenomeSize = mGenomes->mGenomes[lGenome].length();
        int lStartPos = rand() % (lGenomeSize - mReadLength);
                
        cout << ">Genome" << lGenome + mGenomes->mStartNumber << ":pos" << lStartPos << '\t' << endl;
        cout << mGenomes->mGenomes[lGenome].substr( lStartPos, mReadLength ) << endl;
      }
    }
  };
}