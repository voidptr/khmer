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

#define BIT_COUNT_PARTITION 1000

#define HASHES 8
#define BIN_SIZE 50

namespace bleu {
  
  typedef unsigned short SetID;    
  
  using namespace khmer;
  using namespace std;
  
//  class CanonicalSetManager;
  
  class CanonicalSet
  {
    typedef CanonicalSet * SetHandle;
    typedef SetHandle * SetPointer;
    typedef unsigned short SetOffset;
    
    friend class CanonicalSetManager;    
  private:
    
//    unsigned long long KmerCount;
//    unsigned long long TotalKmerCount;
    
    SetOffset PrimarySetOffset;
    
  public:
    vector<SetPointer> BackReferences;
    
  public:
    SetPointer Self;
    
    CanonicalSet( SetOffset aStartingSetOffset = 0 )      
    {        
      PrimarySetOffset = aStartingSetOffset;

      Self = new SetHandle();
      *Self = NULL;
      
      *Self = this;
      
      BackReferences.push_back( Self );
    }

    
    SetOffset GetPrimaryOffset() const // the main one to go to
    {
      return PrimarySetOffset;
    }
    
//    ~CanonicalSet()
//    {
//      *Self = NULL;
//      Self = NULL;
//    }
    

//    struct CompSet {
//      bool operator()(SetHandle * aS1, SetHandle * aS2)
//      {
//        if ( (aS1 == NULL || *aS1 == NULL ) && ( aS2 == NULL || *aS2 == NULL ) )          
//          return false;          
//        else if ( (aS1 == NULL || *aS1 == NULL ) )
//          return true;
//        else if ( ( aS2 == NULL || *aS2 == NULL ) )
//          return false;
//        else            
//          return (*aS1)->GetTotalKmerCount() > (*aS2)->GetTotalKmerCount(); 
//      }
//    };
    
  };
}