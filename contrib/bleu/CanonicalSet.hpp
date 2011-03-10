/*
 *  CanonicalSet.hpp
 *  bleu
 *
 *  Created by Rosangela Canino-Koning on 10/8/10.
 *
 */

#include <list>
#include <vector>

#define BIN_SIZE 50


namespace bleu {
  
  using namespace std;
  
  class CanonicalSet
  {
    typedef CanonicalSet * SetHandle;
    typedef SetHandle * SetPointer;
    typedef unsigned long long SetOffset;
    
    friend class CanonicalSetManager;    

  public:
    vector<SetPointer> BackReferences;
    bool JoinOfConvenience;
    unsigned long long KmerCount;
  private:
    SetOffset PrimarySetOffset;
    
  public:
    SetPointer Self;
    
    CanonicalSet( SetOffset aStartingSetOffset )      
    {        
      PrimarySetOffset = aStartingSetOffset;
      assert( PrimarySetOffset > 0 );
      
      Self = new SetHandle();
      *Self = this;
      
      //BackReferences.push_back( Self );
      
      KmerCount = 0;
      
      JoinOfConvenience = false;
    }
    
    SetOffset GetPrimarySetOffset()
    {
      assert( PrimarySetOffset > 0 );
      return PrimarySetOffset;
    }
    
    void OutputInfo()
    {
      cout << "Set Info: Set=" << this
      << " KmerCount=" << KmerCount
      << " PrimaryOffset=" << PrimarySetOffset
      << " BackRefCt=" << BackReferences.size()
      
      << endl;

    }
    
    void ApplyJoinOfConvenience()
    {
      JoinOfConvenience = true;
    }

  public:
    
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