/*
 *  SequenceHashArbitrary.hpp
 *  bleu
 *
 *  Created by Rosangela Canino-Koning on 11/15/10.
 *
 */

namespace bleu {
 
  // class SequenceHashArbitrary
  //
  // Respresents a DNA sequence of arbitrary length (currently limited to 200)
  // as a 64 bit hash value.
  class SequenceHashArbitrary {
  private:
    
    string sequence;    
    
  public:
    unsigned long long canonical_hash;
    
  public:
    // default constructor, to satisfy the Pair<T> class requirements
    // don't use it.
    SequenceHashArbitrary()
    {
      sequence = "";
      canonical_hash = 0;
    }
  
    // constructor
    // takes a sequence, and hashes it.
    // note, the sequence is used as is, using the length as K.
    SequenceHashArbitrary(string & aSequence, unsigned long long aHash)
    {
      sequence = aSequence;      
      canonical_hash = aHash;
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
    
    // compare the actual sequence
    bool SeqEquals( SequenceHashArbitrary & aOther )
    {
      return ( sequence == aOther.sequence );
    }
  };
}

  
