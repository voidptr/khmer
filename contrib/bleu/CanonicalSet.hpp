/*
 *  CanonicalSet.hpp
 *  bleu
 *
 *  Created by Rosangela Canino-Koning on 10/8/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  Set.hpp
 *  bleu
 *
 *  Created by Rosangela Canino-Koning on 10/4/10.
 *
 */

#include "../../lib/khmer.hh"
#include <list>



#define HASHES 8
#define BIN_SIZE 50


namespace bleu {
  
  using namespace khmer;
  using namespace std;
  
  class CanonicalSet
  {
    typedef CanonicalSet * SetHandle;
    typedef SetHandle * SetPointer;
    typedef unsigned short SetOffset;
    
    friend class CanonicalSetManager;    

  public:
    vector<SetPointer> BackReferences;
    bool JoinOfConvenience;
    SetOffset PrimarySetOffset;
    unsigned long long KmerCount;

    
  public:
    SetPointer Self;
    
    CanonicalSet( SetOffset aStartingSetOffset )      
    {        
      PrimarySetOffset = aStartingSetOffset;
      
      Self = new SetHandle();
      *Self = this;
      
      //BackReferences.push_back( Self );
      
      KmerCount = 0;
      
      JoinOfConvenience = false;
    }
    
    void OutputInfo()
    {
      cout << "Set Info: Set=" << this
      << " KmerCount=" << KmerCount
      << " PrimaryOffset=" << PrimarySetOffset
      << " BackRefCt=" << BackReferences.size()
      
      << endl;

    }
    
//    SetOffset GetPrimaryOffset() const // the main one to go to
//    {
//      return PrimarySetOffset;
//    }
    
    void ApplyJoinOfConvenience()
    {
      JoinOfConvenience = true;
    }

  public:
    
//    void Increment( unsigned long long aCount = 1 )
//    {
//      KmerCount+=aCount;
//    }
    
//    unsigned long long GetKmerCount()
//    {
//      return KmerCount;
//    }
    
    struct CompSet {
      bool operator()(SetPointer aS1, SetPointer aS2)
      {
        if ( (aS1 == NULL || *aS1 == NULL ) && ( aS2 == NULL || *aS2 == NULL ) )          
          return false;          
        else if ( (aS1 == NULL || *aS1 == NULL ) )
          return true;
        else if ( ( aS2 == NULL || *aS2 == NULL ) )
          return false;
        else            
          return (*aS1)->KmerCount > (*aS2)->KmerCount; 
      }
    };      
      
    };
  }