/*
 *  SequenceHashArbitrary.hpp
 *  bleu
 *
 *  Created by Rosangela Canino-Koning on 11/15/10.
 *
 */

// test validity
#define is_valid_dna_character(ch) ( ch == 'A' || ch == 'a' || \
ch == 'C' || ch == 'c' || \
ch == 'G' || ch == 'g' || \
ch == 'T' || ch == 't' )

// bit representation of A/T/C/G.
#define twobit_representation(ch) ((toupper(ch)) == 'A' ? 0LL : \
(toupper(ch)) == 'C' ? 1LL : \
(toupper(ch)) == 'G' ? 2LL : 3LL)

#define reversetwobit_representation(n) ((n) == 0 ? 'A' : \
(n) == 1 ? 'C' : \
(n) == 2 ? 'G' : 'T')

#define complement_twobit(n) ((~n & 3))

#define twobit_complement(ch) ((toupper(ch)) == 'A' ? 3LL : \
(toupper(ch)) == 'T' ? 0LL : \
(toupper(ch)) == 'C' ? 2LL : 1LL)

// choose wisely between forward and rev comp.
#define uniqify_forward_reverse(f, r) ((f) < (r) ? (f) : (r))

namespace bleu {
 
  class SequenceHashArbitrary {
  private:
//    unsigned long long forward_hash;
//    unsigned long long reverse_hash;
    
    string sequence;
    
    static const unsigned long long Primes[199];
    static const unsigned long long Nucleotides[4];
    
  public:
    unsigned long long canonical_hash;
    
  public:
    SequenceHashArbitrary()
    {
//      forward_hash = 0;
//      reverse_hash = 0;
      sequence = "";
      canonical_hash = 0;
    }
  
    SequenceHashArbitrary(string aSequence)//, unsigned char aWordLength)
    {
      sequence = aSequence;      
      canonical_hash = hash(aSequence);//, //word_length, 
    }
    
    bool operator==( SequenceHashArbitrary & aOther)
    {
      return ( aOther.canonical_hash == canonical_hash );
    }
    
    bool operator<( SequenceHashArbitrary & aOther )
    {
      return ( canonical_hash < aOther.canonical_hash );
    }
    
    bool operator>( SequenceHashArbitrary & aOther )
    {
      return ( canonical_hash > aOther.canonical_hash );
    }
    
    bool SeqEquals( SequenceHashArbitrary & aOther )
    {
      return ( sequence == aOther.sequence );
    }
    
  private:

    unsigned long long rotatebitsclockwise( unsigned long long aHash )
    {
      unsigned long long b = aHash >> 1;
      b ^= aHash << 63;
      
      return b;
    } 

    unsigned long long hash(string & aSeq)
    {
      unsigned long long lHash = pow(Primes[0], 
        Nucleotides[twobit_representation(aSeq[0])]);
      unsigned long long lHashRev = pow(Primes[0], 
        Nucleotides[twobit_complement(aSeq[ aSeq.length() - 1 ])]);
      
      lHash = rotatebitsclockwise( lHash );
      lHashRev = rotatebitsclockwise( lHashRev );
      
      for (int i = 1, j = aSeq.length() - 2; i < aSeq.length(); ++i, --j)
      {
        unsigned long long lTimes =  pow(Primes[i], 
          Nucleotides[twobit_representation(aSeq[i])]);
        lHash = mul( lHash, lTimes );

        
        unsigned long long lRevTimes = pow(Primes[i], 
          Nucleotides[twobit_complement(aSeq[j])]); 
        lHashRev = mul( lHashRev, lRevTimes );
        
        lHash = rotatebitsclockwise( lHash );
        lHashRev = rotatebitsclockwise( lHashRev );
      }
      
      unsigned long long lFinalHash = uniqify_forward_reverse(lHash, lHashRev);

      return lFinalHash;
    }
    
    unsigned long long mul(unsigned long long x, unsigned long long y)
    {
      unsigned long long r = 0;
      
      for ( ; y != 0 ; y >>= 1, x <<= 1 )
      {    
        if ( y & 1 )
          r += x;
      } 
      return r;
      
    }
  };
  
  const unsigned long long SequenceHashArbitrary::Primes[199] =
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
    1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 // 199
  };
  
  const unsigned long long SequenceHashArbitrary::Nucleotides[4] = 
  {
    1, 2, 3, 4
  };

}

  
