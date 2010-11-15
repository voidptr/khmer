/*
 *  HashIntoType.hpp
 *  bleu
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
 
  class SequenceHash {
  private:
    unsigned long long forward_hash;
    unsigned long long reverse_hash;
    
    unsigned char word_length;
    
    string sequence;
    
    unsigned long long bitmask;
  public:
    unsigned long long canonical_hash;
    
  public:
    SequenceHash()
    {
      forward_hash = 0;
      reverse_hash = 0;
      sequence = "";
      word_length = 0;
      canonical_hash = 0;
      
      bitmask = 0;
      for (unsigned int i = 0; i < word_length; i++) {
        bitmask = (bitmask << 2) | 3;
      }
    }
  
    SequenceHash(string & aSequence, unsigned char aWordLength)
    {
      forward_hash = 0;
      reverse_hash = 0;

      assert (aWordLength <= 32);
      
      word_length = aWordLength;
      sequence = aSequence;
      
      canonical_hash = hash(aSequence.c_str(), word_length, forward_hash, reverse_hash);  
      
      bitmask = 0;
      for (unsigned int i = 0; i < word_length; i++) {
        bitmask = (bitmask << 2) | 3;
      }
    }
    
    void Hash(string & aSequence, unsigned char aWordLength)
    {
      forward_hash = 0;
      reverse_hash = 0;
      
      assert (aWordLength <= 32);
      
      word_length = aWordLength;
      sequence = aSequence;
      
      canonical_hash = hash(aSequence.c_str(), word_length, forward_hash, reverse_hash);
      
      bitmask = 0;
      for (unsigned int i = 0; i < word_length; i++) {
        bitmask = (bitmask << 2) | 3;
      }
      
      cout << "hai" << endl;
    }
    
    void Hash(const char * aSequence, unsigned char aWordLength)
    {
      forward_hash = 0;
      reverse_hash = 0;
      
      assert (aWordLength <= 32);
      
      word_length = aWordLength;
      sequence = string(aSequence);
      
      canonical_hash = hash(aSequence, word_length, forward_hash, reverse_hash);
      
      bitmask = 0;
      for (unsigned int i = 0; i < word_length; i++) {
        bitmask = (bitmask << 2) | 3;
      }
      
            cout << "hai" << endl;
    }
    
    void MoveForward(const char & aNucleotide)
    {
      canonical_hash = move_hash_forward(forward_hash, reverse_hash, aNucleotide);
      sequence = sequence.substr(1).append(&aNucleotide);
    }
  
  private:
    
//    string _revhash(unsigned long long hash, WordLength k)
//    {
//      std::string s = "";
//      
//      unsigned int val = hash & 3;
//      s += reversetwobit_representation(val);
//      
//      for (WordLength i = 1; i < k; i++) {
//        hash = hash >> 2;
//        val = hash & 3;
//        s += reversetwobit_representation(val);
//      }
//      
//      reverse(s.begin(), s.end());
//      
//      return s;
//    }
    
    unsigned long long hash(const char * kmer, const WordLength k, unsigned long long& _h, unsigned long long& _r)
    {
      // sizeof(HashIntoType) * 8 bits / 2 bits/base  
      assert(k <= sizeof(HashIntoType)*4);
      
      HashIntoType h = 0, r = 0;
      
      h |= twobit_representation(kmer[0]);
      for (WordLength i = 1; i < k; i++) {
        h = h << 2;
        h |= twobit_representation(kmer[i]);
      }
      
      r = ~h; // bit complement
      
      // swap consecutive pairs
      r = ((r >> 2)  & 0x3333333333333333ULL) | ((r & 0x3333333333333333ULL) << 2);
      // swap nibbles 
      r = ((r >> 4)  & 0x0F0F0F0F0F0F0F0FULL) | ((r & 0x0F0F0F0F0F0F0F0FULL) << 4);
      // swap bytes
      r = ((r >> 8)  & 0x00FF00FF00FF00FFULL) | ((r & 0x00FF00FF00FF00FFULL) << 8);
      // swap 2-byte words
      r = ((r >> 16) & 0x0000FFFF0000FFFFULL) | ((r & 0x0000FFFF0000FFFFULL) << 16);
      // swap 2-word pairs
      r = ( r >> 32                      ) | ( r                       << 32);
      
      _h = h;
      _r = r;
      
      return uniqify_forward_reverse(h, r);
    }    
    
    unsigned long long move_hash_forward( unsigned long long & aOldForwardHash, unsigned long long & aOldReverseHash, const char & aNextNucleotide )
    {
      unsigned long long lTwoBit = twobit_representation( aNextNucleotide );
      
      aOldForwardHash = aOldForwardHash << 2; // left-shift the previous hash over
      aOldForwardHash |= lTwoBit; // 'or' in the current nucleotide
      aOldForwardHash &= bitmask; // mask off the 2 bits we shifted over.
      
      // now handle reverse complement
      aOldReverseHash = aOldReverseHash >> 2;
      aOldReverseHash |= (complement_twobit(lTwoBit) << (word_length*2 - 2));
      
      // pick the better bin of the forward or reverse hashes
      return uniqify_forward_reverse(aOldForwardHash, aOldReverseHash);
    }
  };
}