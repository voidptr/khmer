///*
// *  ChainedSet.hpp
// *  bleu
// *
// *  Created by Rosangela Canino-Koning on 10/4/10.
// *
// */
//
//#include "../../lib/khmer.hh"
//#include <list>
//
//#define BIT_COUNT_PARTITION 1000
//
//#define HASHES 8
//#define BIN_SIZE 50
//
//namespace bleu {
//
//  typedef unsigned short SetID;    
//  
//  using namespace khmer;
//  using namespace std;
//    
//  class ChainedSetsManager;
//  
//  class ChainedSet
//  {
//    typedef Set * SetHandle;
//    typedef unsigned short SetOffset;
//    
//    friend class ChainedSetsManager;    
//  private:
//    
//    unsigned long long KmerCount;
//    unsigned long long TotalKmerCount;
//    
//    unsigned long long CurrentBin;
//    list<SetHandle*>::iterator BinLocation; // this will be invalid until the SetsManager assigns it.
//    
//    unsigned int NodeCount;
//    
//    vector<HashIntoType> Hashes; // in this set
//    bool StoringHashes; // this vector contains stuff
//    
//    vector<SetHandle> ChildSets;
//    vector<SetHandle> FosterSets;
//    Set* ParentSet;
//    
//    vector<SetOffset> SetOffsetsIOccupy;
//    
//    SetHandle Root; // the top node in this tree, the main gateway node to the tree    
//    SetHandle BiologicalPatriarch; // the top of this dingle-dangle of this tree. If we're fostered, this isn't the root node.
//    
//    SetHandle * Self;
//    
//    bool AmInBin;
//        
//  public:
//    ChainedSet( SetOffset aStartingSetOffset = 0 )      
//    {        
//      SetOffsetsIOccupy.push_back( aStartingSetOffset );
//      
//      KmerCount = 0;
//      TotalKmerCount = 0; // this stays zero unless you are a root or a patriarch
//      NodeCount = 1; // this stays 1 unless you are a root or a patriarch.
//
//      CurrentBin = 0;
//      
//      Root = this; // default that you are your own root set. You aren't anyone's child or foster child.
//      BiologicalPatriarch = this; // I'm my own grandpa!
//      
//      ParentSet = NULL; // I'm root      
//      
//      StoringHashes = true;
//      
//      Self = new SetHandle();
//      *Self = this;
//      
//      AmInBin = false;
//    }
//    
//    void output_info()
//    {
//      cout << " ID=" << GetPrimaryGatewaySetOffset() << " Addr=" << this << " Bin=" << CurrentBin << " TotKmCt=" << GetTotalKmerCount() << " TotNdCt=" << GetTreeSize() << " TotFsCt=" << GetTotalFosterCount();
//      
//      if ( AmRoot() )
//        cout << " ROOT";
//      
//      if ( AmFostered() )
//        cout << " FOSTERED";
//      
//      if ( AmFosterRoot() )
//        cout << " FOSTEREDROOT";
//      
//      if ( AmBiologicalPatriarch() )
//        cout << " BIOLOGICALPATRIARCH";
//      
//      cout << endl;
//      
//      if ( !AmFostered() ) // only do this for the root        
//        for ( int i = 0; i < BiologicalPatriarch->FosterSets.size(); ++i )
//        {
//          cout << " Foster";
//          BiologicalPatriarch->FosterSets[i]->output_info();
//        }
//    }
//    
//    SetOffset GetPrimaryGatewaySetOffset() const // the main one to go to
//    {
//      return Root->SetOffsetsIOccupy[0]; // pull the root set's primary one.
//    }
//    
//    ~ChainedSet()
//    {
//      *Self = NULL;
//      Self = NULL;
//    }
//    
//    bool AmRoot()
//    {
//      if ( ParentSet == NULL )
//      {
//        assert( Root == this && BiologicalPatriarch == this ); // we didn't fuck this up.
//        return true;
//      }
//
//      return false;
//    }
//    
//    bool AmElgibleForFostering()
//    {
//      return AmRoot() && !AmFostered();
//    }
//    
//    bool AmFostered() // this is either if you're part of a foster tree, or the root of that branch.
//    {
//      if ( Root != BiologicalPatriarch )
//      {
//        assert( StoringHashes ); // if you're fostered, you should always be storing hashes.
//        return true;
//      }
//        
//      
//      return false;
//    }
//    
//    bool AmFosterRoot()
//    {
//      if ( AmFostered() && AmBiologicalPatriarch() )
//      {
//        assert ( ParentSet != NULL ); // I'd better have a parent.
//        return true;
//      }
//        
//      
//      return false;
//    }
//    
//    bool AmBiologicalPatriarch()
//    {
//      if ( BiologicalPatriarch == this )
//        return true;
//      
//      return false;
//    }
//    
//    bool HaveSpaceForFoster()
//    {
//      if ( !AmFostered() && GetTotalFosterCount() < 50 )
//        return true;
//      
//      return false;
//    }
//    
//    int GetTreeSize()
//    {
//      return BiologicalPatriarch->NodeCount; // this always works.
//    }
//    
//    int GetTotalKmerCount()
//    {
//      return BiologicalPatriarch->TotalKmerCount;
//    }
//    
//    int GetTotalFosterCount()
//    {
//      return BiologicalPatriarch->FosterSets.size();
//    }
//    
//    unsigned long long GetIndividualKmerCount()
//    {
//      return KmerCount; // um, shrug.
//    }
//    
//    void IncrementCount( int lAmount = 1)
//    {
//      KmerCount += lAmount;
//      
//      BiologicalPatriarch->TotalKmerCount += lAmount; // keep the dude up to date.
//    }
//    
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
//    
//  };
//}