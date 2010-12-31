/*
 *  SequenceHashArbitrary_Lookup.hpp
 *  bleu
 *
 *  Created by Rosangela Canino-Koning on 12/26/10.
 *
 */

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <assert.h>

using namespace std;

#define MAX_K_SIZE 200

#define CHUNK_SIZE 1
#define CHUNK_PERMUTATIONS 4 //pow(4, CHUNK_SIZE)

//#define CHUNK_SIZE 2
//#define CHUNK_PERMUTATIONS 16 //pow(4, CHUNK_SIZE)

//#define CHUNK_SIZE 4
//#define CHUNK_PERMUTATIONS 256 //pow(4, CHUNK_SIZE)

#define PRIME_GROUPS MAX_K_SIZE/CHUNK_SIZE

namespace bleu {
  class SequenceHashArbitrary_Lookup {
  
  private:  
    unsigned long long CurrentPrimes[CHUNK_SIZE];
    static const unsigned long long Primes[200];
    static const unsigned long long Nucleotides[4];
  
  public:
    unsigned long long PrimesLookup[PRIME_GROUPS][CHUNK_PERMUTATIONS];

    SequenceHashArbitrary_Lookup()
    {
      assert (CHUNK_PERMUTATIONS == pow(4, CHUNK_SIZE) );

      cout << "Generating Lookup Table" << endl;
      memset(CurrentPrimes, 0, sizeof(unsigned long long) * CHUNK_SIZE);
      memset(PrimesLookup, 0, sizeof(unsigned long long) * PRIME_GROUPS * CHUNK_PERMUTATIONS);

      for ( int i = 0; i < PRIME_GROUPS; ++i )
      {
        for ( int k = 0; k < CHUNK_SIZE; ++k )
          CurrentPrimes[k] = Primes[(i * CHUNK_SIZE) + k];

        for ( int j = 0; j < CHUNK_PERMUTATIONS; ++j )
        {
          unsigned long long lNucleotide = j; // AA, AC, AG, AT, CA, CC, etc.
        
          PrimesLookup[i][j] = 1;
          for ( int l = 0; l < CHUNK_SIZE; ++l )
          {      
            PrimesLookup[i][j] *= (unsigned long long ) pow( CurrentPrimes[l], Nucleotides[ lNucleotide & 3 ] );
            
            lNucleotide >>= 2; // move the next proper one to the left most position.
          }
          
        }
      }
      cout << " done." << endl;
    }
  };

  const unsigned long long SequenceHashArbitrary_Lookup::Primes[200] =
  {
    3, 5, 7, 11, 13, 17, 19, 23, 29, // 9
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71, // 19
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, // 29
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173, // 39 
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229, // 49
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281, // 59
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349, // 69
    353, 359, 367, 373, 379, 383, 389, 397, 401, 409, // 79
    419, 421, 431, 433, 439, 443, 449, 457, 461, 463, // 89
    467, 479, 487, 491, 499, 503, 509, 521, 523, 541, // 99
    547, 557, 563, 569, 571, 577, 587, 593, 599, 601, // 109
    607, 613, 617, 619, 631, 641, 643, 647, 653, 659, // 119
    661, 673, 677, 683, 691, 701, 709, 719, 727, 733, // 129
    739, 743, 751, 757, 761, 769, 773, 787, 797, 809, // 139
    811, 821, 823, 827, 829, 839, 853, 857, 859, 863, // 149
    877, 881, 883, 887, 907, 911, 919, 929, 937, 941, // 159
    947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, // 169
    1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, // 179 
    1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, // 189
    1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229 // 200
  };

  const unsigned long long SequenceHashArbitrary_Lookup::Nucleotides[4] = 
  {
    1, 2, 3, 4
  };

}

