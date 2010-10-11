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
  
  using namespace khmer;
  using namespace std;
  
  class CanonicalSet
  {
    typedef CanonicalSet * SetHandle;
    typedef SetHandle * SetPointer;
    typedef unsigned short SetOffset;
    
    friend class CanonicalSetManager;    
  private:
    
    SetOffset PrimarySetOffset;
    unsigned long long KmerCount;
    
    bool AmFoster;
    bool StoringHashes;
    
  public:
    SetHandle Parent;
    vector<SetHandle> Fosters;
    set<HashIntoType> Hashes;
    vector<SetPointer> BackReferences;
    
  public:
    SetPointer Self;
    
    CanonicalSet( SetOffset aStartingSetOffset = 0, bool aAmFoster = false, SetHandle aParent = NULL )      
    {        
      PrimarySetOffset = aStartingSetOffset;
      
      Self = new SetHandle();
      *Self = NULL;
      
      *Self = this;
      
      BackReferences.push_back( Self );
      
      KmerCount = 0;
      
      AmFoster = aAmFoster;
      Parent = aParent;
      
      StoringHashes = true; // everyone starts out storing hashes.
    }
    
    void OutputInfo()
    {
      cout << "Set Info: Set=" << this
      << " AmFosterChild=" << AmFosterChild() 
      << " Parent=" << Parent
      << " AmStoringHashes=" << AmStoringHashes()
      << " HashCt=" << Hashes.size()
      << " KmerCount=" << KmerCount
      << " PrimaryOffset=" << GetPrimaryOffset()
      << " BackRefCt=" << BackReferences.size() << endl;
      
      if ( !AmFoster )
      {
        for( int i = 0; i < Fosters.size(); ++i )
        {
          cout << " ";
          Fosters[i]->OutputInfo();
        }
      }
    }
    
    SetOffset GetPrimaryOffset() const // the main one to go to
    {
      if ( AmFoster )
        return Parent->GetPrimaryOffset();
      
      return PrimarySetOffset;
    }
    
    bool AmStoringHashes()
    {
      return StoringHashes;
    }
    
    bool AmFosterChild()
    {
      return AmFoster;
    }
    
    void AddToSet( HashIntoType aHash ) // this can quietly fail, no big.
    {
      if ( StoringHashes )
        Hashes.insert( aHash );
    }
    
    SetHandle FindResponsibleSet( HashIntoType aHash )
    {
      if ( Fosters.size() == 0 )
        return this; // nobody here but us chickens, (so it must be me)
      else
      {
        for ( int i = 0; i < Fosters.size(); ++i )
        {
          if ( Fosters[i]->HaveThisHashSaved( aHash ) )
            return Fosters[i];
        }
      }
      return this; // wasn't one of them, so it's me.
    }
    
  private:
    bool HaveThisHashSaved( HashIntoType aHash )
    {
      if ( Hashes.find( aHash ) != Hashes.end() )
        return true; // I've got it
      
      return false;
    }
    
  public:
    
    bool AcceptFosterChild( SetHandle aChild )
    { 
      if ( AmFoster ) // can't very well accept a foster child if I"m a foster myself.
      {
        cout << "AcceptFosterChild, AmFoster" << endl; 
        cout << "This=";
        OutputInfo();
        cout << "Child=";
        aChild->OutputInfo();
        return false;
      }
      
      if ( aChild->BecomeFoster( this ) ) // if it succeeded
      {
        Fosters.push_back( aChild );
        return true;
      }
      
      cout << "AcceptFosterChild, Child Couldn't become Foster" << endl;
      cout << "This=";
      OutputInfo();
      cout << "Child=";
      aChild->OutputInfo();
      return false;
    }
    
    bool TakeFosterChild( SetHandle aChild ) // take it from someone else
    {
      if ( AmFoster ) // can't accept a foster if I'm a foster myself
      {
        cout << "TakeFosterChild, AmFoster" << endl;     
        cout << "This=";
        OutputInfo();
        cout << "Child=";
        aChild->OutputInfo();
        assert(0); //return false;
      }
      
      if ( !aChild->AmFoster ) // it's gotta be someone else's kid already
      {
        cout << "TakeFosterChild, Child !AmFoster" << endl; 
        cout << "This=";
        OutputInfo();
        cout << "Child=";
        aChild->OutputInfo();
        assert(0); //return false;
      }
      
      bool lWorked = false;
      for ( int i = 0; i < aChild->Parent->Fosters.size(); ++i )
      {
        if ( aChild->Parent->Fosters[i] == aChild )
        {
          aChild->Parent->Fosters.erase( aChild->Parent->Fosters.begin() + i );
          lWorked = true;
          break;
        }
      }
      
      if ( !lWorked )
      {
        cout << "TakeFosterChild, Foster Release Failed" << endl; 
        cout << "This=";
        OutputInfo();
        cout << "Child=";
        aChild->OutputInfo();
        assert(0); // return false;
      }
      
      Fosters.push_back( aChild );
      aChild->Parent = this; // reassign the parent.
      
      return true;
    }
    
    bool EmancipateFosterChild( SetHandle aChild, SetOffset aNewSetOffset )
    {
      if ( AmFoster || Fosters.size() == 0 ) // can't very well remove a foster child I don't have
      {
        cout << "EmancipateFosterChild - Foster.size() == 0 || AmFoster" << endl;
        cout << " This=";
        OutputInfo();
        cout << " Child=";
        aChild->OutputInfo();
        cout << " NewOffset=" << aNewSetOffset << endl;
        assert(0);  //  return false;
      }
      
      for ( int i = 0; i < Fosters.size(); ++i )
      {
        if ( Fosters[i] == aChild ) // found it
        {
          if ( aChild->BecomeIndependent( aNewSetOffset ) ) // if it succeeded
          {
            Fosters.erase( Fosters.begin() + i );
            return true;
          }
          else // removal failed. wut?
          {
            cout << "EmancipateFosterChild, Foster Release Failed" << endl; 
            cout << " This=";
            OutputInfo();
            cout << " Child=";
            aChild->OutputInfo();
            cout << " NewOffset=" << aNewSetOffset << endl;
            assert(0);   //  return false;
          }
        }
      }      
      
      cout << "EmancipateFosterChild - Couldn't find the child to emancipate" << endl;
      cout << " This=";
      OutputInfo();
      cout << " Child=";
      aChild->OutputInfo();
      cout << " NewOffset=" << aNewSetOffset << endl;
      assert(0);   //  return false; // couldn't find it.
    }
    
  private:
    
    bool BecomeFoster( SetHandle aParent )
    {
      if ( Fosters.size() > 0 || AmFoster || !StoringHashes ) // can't become someone's foster if I'm already fostering, or I'm already fostered to someone else.
      {
        cout << "BecomeFoster - Foster.size() > 0 || AmFoster || !StoringHashes" << endl;
        cout << " This=";
        OutputInfo();
        cout << " Parent=";
        aParent->OutputInfo();
        assert(0); // return false;
      }
      
      PrimarySetOffset = 0;
      AmFoster = true;
      Parent = aParent;
      
      return true;
    }
    
    bool BecomeIndependent( SetOffset aNewSetOffset )
    {
      if ( !AmFoster ) // can't become free if you're already free
      {
        cout << "BecomeIndependent, !AmFoster" << endl;
        cout << " This=";
        OutputInfo();
        cout << " NewOffset=" << aNewSetOffset << endl;
        assert(0); //return false;
      }
      
      PrimarySetOffset = aNewSetOffset;
      AmFoster = false;
      Parent = NULL;
      
      return true;
    }
    
  public:
    
    bool StopStoringHashes()
    {
      if ( AmFoster )
      {
        cout << "StopStoringHashes, !AmFoster" << endl;
        cout << " This=";
        OutputInfo();
        assert(0); //return false; // nowai
      }
      
      StoringHashes = false;
      Hashes.clear();
      
      return true;
    }
    
    bool ShouldStopStoringHashes()
    {
      if ( Hashes.size() > 100 ) // I'm really too big.
        return true;
      
      return false; // I'm all set, thanks.
    }
    
  public:
    
    void Increment( unsigned long long aCount = 1 )
    {
      KmerCount+=aCount;
    }
    
    unsigned long long GetKmerCount()
    {
      return KmerCount;
    }
    
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
    
    bool AmValid()
    {
      if ( AmFoster )
      {
        if ( Fosters.size() > 0 )
        {
          cout << "AmValid, Foster - Have Foster Children" << endl;
          cout << " This=";
          OutputInfo();
          assert(0); 
        }
        //  return false;
        
        if ( Parent == NULL )
        {
          cout << "AmValid, Foster - No Parent" << endl;
          cout << " This=";
          OutputInfo();
          
          assert(0); //
          //return false;
          
          if ( !StoringHashes )
          {
            cout << "AmValid, Foster - Not Storing Hashes" << endl;
            cout << " This=";
            OutputInfo();
            
            assert(0); //
            //return false;
          }
          
          if ( !FoundInParent() )
          {
            cout << "AmValid, Foster - Not Found in Parent" << endl;
            cout << " This=";
            OutputInfo();
            assert(0); //
            //return false;
          }
        }
        else
        {
          if ( Parent != NULL )
          {
            cout << "AmValid, !Foster - Have Parent" << endl;
            cout << " This=";
            OutputInfo();
            assert(0); //
            //return false;
          }
        }   
      }
      if ( StoringHashes )
        if ( KmerCount != Hashes.size() )
        {
            cout << "AmValid, Storing Hashes, but Unmatched KMer to Hash Count" << endl;
            cout << " This=";
            OutputInfo();
            assert(0); //
            //return false; 
          }
        
      return true;          
      
    }
      
    private:
      
      bool FoundInParent()
      {
        if ( !AmFoster )
        {
          cout << "FoundInParent, Called by Not AmFoster" << endl;
          cout << " This=";
          OutputInfo();
          assert( 0 );
        }
        
        
        for ( int i = 0; i < Parent->Fosters.size(); ++i )
        {
          if ( Parent->Fosters[i] == this )
            return true;
        }
        
        cout << "FoundInParent, Not Found" << endl;
        cout << " This=";
        OutputInfo();      
        return false;
      }  
      
      
    };
  }