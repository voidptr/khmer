/*
 *  CanonicalSetManager.h
 *  bleu
 *
 *  Created by Rosangela Canino-Koning on 10/8/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "../../lib/hashtable.hh"
#include "../../lib/parsers.hh"
#include "../external_lib/cBitArray.h"
#include "CanonicalSet.hpp"
#include <algorithm>

#define BIT_COUNT_PARTITION 1000

#define HASHES 8
#define BIN_SIZE 50
#define FOSTER_COUNT_MULTIPLIER 5

#define CACHESIZE 10

#define SETS_SIZE 65535

namespace bleu {
  
  using namespace khmer;
  using namespace std;
  
  typedef CanonicalSet * SetHandle;
  typedef SetHandle * SetPointer;
  
  typedef unsigned long long HashBin;
  typedef unsigned long long SetOffsetBin;
  typedef unsigned short SetOffset;
  
  class CanonicalSetsManager
  {
  private:
    // sizes set during constructor
    cBitArray * _hash_table_preliminary[HASHES]; // two-dimensional hash-table
    cBitArray * _hash_table[HASHES]; // two-dimensional hash-table
    
    unsigned long long _tablesizes[HASHES]; // the sizes of the hash tables
    unsigned long long * _hash_table_bit_counts_lookup[HASHES]; // the number of bits set in the hash table, by every 10k entries
    unsigned int _last_set_offset;
    
    // sizes set during prep 1 (based on processing of _hash_table)
    unsigned long long _hash_table_total_bit_counts[HASHES];
    unsigned short * _set_offsets[HASHES]; 
    
    // contents set during (based on the processing of the reads)  
    vector<SetPointer> _sets; // array of sets
    
    
    // caching system so I don't have to re-modulo. Dunno if this is faster than modulo or not. will test.
    pair<HashIntoType, HashBin> _hash_bin_cache[HASHES][CACHESIZE];
    int _hash_bin_cache_last_used_index[HASHES];
    
    // caching system so I don't have to keep re-counting the bits. Dunno if this is faster than counting or not.
    pair<HashBin, SetOffsetBin> _set_offset_bin_cache[HASHES][CACHESIZE];
    int _set_offset_bin_cache_last_used_index[HASHES];
    
    vector<SetOffset> _released_set_offsets;
    
    //map<unsigned long long, list<SetHandle*> > _binned_sets; // our sets, sorted by size, in heap form
    vector<SetPointer> _sorted_sets;
    
  public:
    CanonicalSetsManager( unsigned long long aMaxMemory )
    {
      _last_set_offset = 0; // init to zero
      
      _tablesizes[0] = get_first_prime_below( aMaxMemory / HASHES );
      for ( int i = 1; i < HASHES; ++i )      
      {
        _tablesizes[i] = get_first_prime_below(_tablesizes[i-1]); 
      }
      
      for ( int j = 0; j < HASHES; ++j )      
      {
        _hash_table_preliminary[j] = new cBitArray( _tablesizes[j] );
        _hash_table_preliminary[j]->Clear();      
        
        _hash_table[j] = new cBitArray( _tablesizes[j] );
        _hash_table[j]->Clear(); 
        
        _hash_table_bit_counts_lookup[j] = new unsigned long long[(_tablesizes[j] / BIT_COUNT_PARTITION)+1];
        memset(_hash_table_bit_counts_lookup[j], 0, ((_tablesizes[j] / BIT_COUNT_PARTITION)+1) * sizeof(unsigned long long));
        
        // THESE WILL GET SET LATER.            
        _hash_table_total_bit_counts[j] = 0;
        _set_offsets[j] = NULL; // a table the size of the number of bits in the hash table.
        
        _hash_bin_cache_last_used_index[j] = 0;
        memset(_hash_bin_cache[j], 0, sizeof(pair<HashIntoType,HashBin>) * CACHESIZE);
        
        _set_offset_bin_cache_last_used_index[j] = 0;
        memset(_set_offset_bin_cache[j], 0, sizeof(pair<HashBin,SetOffsetBin>) * CACHESIZE);
      }
      
      _sets.resize(SETS_SIZE, NULL);
    }
    //
    // first pass through the reads -- determine if reads are interesting 
    //   enough to warrant getting space to call out a set
    //    
    
    // mark a kmer/hash as "interesting" if it appears more than once.
    void seen_hash( HashIntoType aHash )
    {
      for (int i = 0; i < HASHES; ++i )
      {
        unsigned long long lHashBin = HashToHashBinCached(aHash, i); 
        
        if ( _hash_table_preliminary[i]->Get(lHashBin) == true )
          _hash_table[i]->Set(lHashBin, true);
        else
          _hash_table_preliminary[i]->Set(lHashBin, true);
      }
    }
    
    //
    // second pass through the reads -- actually perform the set assignments
    //
    
    // check whether a kmer can have a set. (I keep wanting to type can_has_set) Damned lolcats.
    bool can_have_set( HashIntoType aHash )
    {
      bool lCanHaveASet = false;
      for ( int i = 0; i < HASHES; ++i )
      {
        HashBin lHashBin = HashToHashBinCached(aHash, i);
        
        if ( i == 0 )
          lCanHaveASet = _hash_table[i]->Get(lHashBin) == true;
        else if ( _hash_table[i]->Get(lHashBin) != lCanHaveASet ) // we have a disagreement.
          return false;
      }
      
      return lCanHaveASet;
    }    
    // find a set for this hash. if there isn't one, create it.
    SetHandle get_set( HashIntoType aHash )
    { 
      SetHandle lHandle = NULL;
      
      if ( has_existing_set( aHash ) )
        lHandle = get_existing_set( aHash );
      else
      {
        lHandle = create_set();
        add_to_set(lHandle, aHash);
      }
      
      
      return lHandle;
    }
    // check whether a kmer has already been assigned a set.
    bool has_existing_set( HashIntoType aHash )
    {
      for ( int i = 0; i < HASHES; ++i )
      {
        SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBinCached(aHash, i), i );
        if ( _set_offsets[i][ lBin ] == 0 )
          return false;
      }
      
      return true;
    }
    // find a set for this hash that already exists.
    SetHandle get_existing_set( HashIntoType aHash )
    {      
      map<SetHandle, int> lRepresented;
      
      for ( int i = 0; i < HASHES; ++i )
      {
        SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBinCached(aHash, i), i );
        
                
        if ( !(_set_offsets[i][ lBin ] > 0) )
        {
          cout << "get_existing_set - set offset is 0" << endl;
          cout << "i=" << i;
          cout << "lBin=" << lBin;
          assert( 0 ); // this should work
        }
        
        lRepresented[ SetOffsetToSet( _set_offsets[i][ lBin ] )->FindResponsibleSet( aHash ) ]++;
      }
      
      // figure out what the consensus set was
      SetHandle lMostRepresented = NULL; 
      for ( map<SetHandle, int>::iterator lSet = lRepresented.begin(); lSet != lRepresented.end(); ++lSet ) // count the votes
      {
        if ( lMostRepresented == NULL || lSet->second > lRepresented[ lMostRepresented ] ) // because maps are sorted, this will end up being the set with the lowest pointer, and the highest count.
          lMostRepresented = lSet->first;
      }     
      
      return lMostRepresented;      
    }
    
    // find a bucket for this hash that already exists. Identical to get_existing_set, but without the drill-down
    SetHandle get_existing_bucket( HashIntoType aHash )
    {      
      map<SetHandle, int> lRepresented;
      
      for ( int i = 0; i < HASHES; ++i )
      {
        SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBinCached(aHash, i), i );
        
        if ( !(_set_offsets[i][ lBin ] > 0) )// somehow we fucked this up.
        {
          cout << "get_existing_buck - set offset is 0" << endl;
          cout << "i=" << i;
          cout << "lBin=" << lBin;
          assert( 0 ); // this should work
        }
        
        lRepresented[ SetOffsetToSet( _set_offsets[i][ lBin ] ) ]++;
      }
      
      // figure out what the consensus set was
      SetHandle lMostRepresented = NULL; 
      for ( map<SetHandle, int>::iterator lSet = lRepresented.begin(); lSet != lRepresented.end(); ++lSet ) // count the votes
      {
        if ( lMostRepresented == NULL || lSet->second > lRepresented[ lMostRepresented ] ) // because maps are sorted, this will end up being the set with the lowest pointer, and the highest count.
          lMostRepresented = lSet->first;
      }     
      
      return lMostRepresented;      
    }
    bool sets_are_disconnected( SetHandle aSet1, SetHandle aSet2 )
    {
      if ( aSet1 == aSet2 )
        return false; // they're already in the same set.
      else
        return true; // one of them is 
    }    
    // put together two sets.
    SetHandle bridge_sets( SetHandle aEncounteredSet, SetHandle aOriginatingSet )
    {
      if ( !( aOriginatingSet != aEncounteredSet) )// we really really shouldn't be the same set.
      {
        cout << "bridge set - originating set == encoutnered set" << endl;
        cout << "originating set=";
        aOriginatingSet->OutputInfo();
        cout << "encountered set=";
        aEncounteredSet->OutputInfo();
        assert( 0 ); // this should work
      }
      
      
      SetHandle lDominatingSet = NULL;

      // four scenarios, 
      // O->!f, E->!f -- lower offset dominates
      // O->f, E->f -- lower offset parent dominates
      // O->!f, E->f -- O dominates
      // O->f, E->!f -- E dominates

      if ( (!aOriginatingSet->AmFosterChild() && !aEncounteredSet->AmFosterChild()) || 
           (aOriginatingSet->AmFosterChild() && aEncounteredSet->AmFosterChild() ) ) // both or neither are fostered
      {
        if ( aOriginatingSet->GetPrimaryOffset() < aEncounteredSet->GetPrimaryOffset() )
        {
          join( aOriginatingSet, aEncounteredSet );
          lDominatingSet = aOriginatingSet;
        }
        else
        {
          join( aEncounteredSet, aOriginatingSet );
          lDominatingSet = aEncounteredSet;
        }
      }
      else if ( !aOriginatingSet->AmFosterChild() && aEncounteredSet->AmFosterChild() ) // originating set isn't fostered.
      {
        join( aOriginatingSet, aEncounteredSet );
        lDominatingSet = aOriginatingSet;
      }
      else // encountered set isn't fostered.
      {
        join ( aEncounteredSet, aOriginatingSet );
        lDominatingSet = aEncounteredSet;
      }
          
      return lDominatingSet;
    }
    void canonicalize()
    {
      // go through all the set offsets and point them to their canonical locations.
      for ( int i = 0; i < HASHES; ++i )
      {
        for ( int j = 0; j < _hash_table_total_bit_counts[i]; ++j )
        {
          if ( _set_offsets[i][j] != 0 ) // there's something in here
          {
            _set_offsets[i][j] = SetOffsetToSet( _set_offsets[i][j] )->GetPrimaryOffset(); 
            
            //assert ( SetOffsetToSet( _set_offsets[i][j] )->AmValid() );
            if ( !( SetOffsetToSet( _set_offsets[i][j] )->AmValid() ) )
            {
              cout << "canonicallize - pointed to set isn't valid" << endl;
              cout << "i=" << i << " j=" << j;
              cout << " pointedset=";
              SetOffsetToSet( _set_offsets[i][j] )->OutputInfo();
              assert( 0 ); // this should work
            }
          }
        }
      }
      
      _released_set_offsets.clear(); // clear this.
      
      for ( int k = 1; k < SETS_SIZE; ++k )
      {
        if ( _sets[ k ] != NULL ) // something in this slot
        {
          SetHandle lSet = *_sets[k];
          
          if ( !lSet->AmValid() )
          {
            cout << "canonicallize - lSet isn't valid" << endl;
            cout << "k=" << k;
            cout << "lSet=";
            lSet->OutputInfo();
            assert( 0 ); // this should work
          }
          
          
          if ( k != lSet->GetPrimaryOffset() )
          {
            _sets[k] = NULL;
            _released_set_offsets.push_back( k ); // it's null! hey!
            //cout << k << endl;
          }
          else // we'll only hit each set once this way.
          {
            //cout << k << endl;
            lSet->BackReferences.clear();
            lSet->BackReferences.push_back( lSet->Self ); // readd the canonical one.
          }
        }
        else
        {
          //cout << k << endl;
          _released_set_offsets.push_back( k ); // it's null! hey!
        }
      }
      
      cout << "released " << _released_set_offsets.size() << " sets." << endl;
    }
    
#define MAX_FOSTERS_ADD 2
#define TOO_MANY_FOSTERS_COUNT MAX_FOSTERS_ADD*2
    
    void reclaim_and_rebalance()
    {
      // now, go through and foster up the tiny sets, and release their set offsets.
      vector<SetOffset> lPotentialFosters;
      vector<SetOffset> lPotentialParents;
      int lFosterChildCount = 0;
      int lFosterParentCount = 0;
      for ( int l = 1; l < SETS_SIZE; ++l )
      {
        if ( _sets[l] != NULL )
        {
          SetHandle lPotential = *_sets[l];
          
          if ( lPotential->AmStoringHashes() && !lPotential->ShouldStopStoringHashes() && lPotential->Fosters.size() == 0 ) // you're small enough, and you're not holding anyone
          {
            lPotentialFosters.push_back( l );
          }
          else if ( lPotential->AmStoringHashes() && lPotential->Fosters.size() < MAX_FOSTERS_ADD ) // you're small enough, and you're not holding too many, but more than one.
          {
            lPotentialParents.push_back( l );
          }
        }
      }
      
      // if we're short on parents, volunteer some fosters.
      if ( lPotentialFosters.size() > 11 && lPotentialParents.size() < lPotentialFosters.size() / MAX_FOSTERS_ADD )
      {
        int lNeededCount = (lPotentialFosters.size() / MAX_FOSTERS_ADD) - lPotentialParents.size();
        for ( int q = 0; q < lNeededCount; ++q )
        {
          lPotentialParents.push_back( lPotentialFosters.back() );
          lPotentialFosters.pop_back();
        }
      }
      
      lFosterChildCount = lPotentialFosters.size();
      lFosterParentCount = lPotentialParents.size();
                        
      // we've crossed a threshold, where we have a parent, and a kid to fill it. This may happen only once, but whatever.
      while ( lPotentialParents.size() > 0 && lPotentialFosters.size() > 0 )
      {
        SetOffset lParentOffset = lPotentialParents.back();
        SetHandle lParent = *_sets[lParentOffset];
        int lSlotsToFill = MAX_FOSTERS_ADD - lParent->Fosters.size();
        
        // while there's room, and kids to fill it
        for ( int m = 0; lPotentialFosters.size() > 0 && m < lSlotsToFill; ++m, lPotentialFosters.pop_back() )
        {
          int lChildOffset = lPotentialFosters.back();            
          SetHandle lChild = *_sets[lChildOffset];
          
          // merge it
          if ( !lParent->AcceptFosterChild( lChild ) )
          {
            cout << "Reclaim and Rebalance -- Failed to accept foster child" << endl;
            cout << "Parent=";
            lParent->OutputInfo();
            cout << "Child=";
            lChild->OutputInfo();
            assert( 0 ); // this should work
          }
          
          if ( !(lParent->AmValid() && lChild->AmValid()) )
          {
            cout << "Reclaim and Rebalance -- either parent or child is invalid" << endl;
            cout << "Parent=";
            lParent->OutputInfo();
            cout << "Child=";
            lChild->OutputInfo();
            assert ( 0 );
          }
          
          // now, point all the kid's hashes to the new parent offset.
          for ( int n = 0; n < HASHES; ++n )
          {
            for ( set<HashIntoType>::iterator lIt = lChild->Hashes.begin(); lIt != lChild->Hashes.end(); ++lIt )
            {
              SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBinCached( *lIt, n ), n);
              _set_offsets[n][lBin] = SetOffsetToSet( _set_offsets[n][lBin] )->GetPrimaryOffset();  // canonicalize
              
              if (!( _set_offsets[n][ lBin ] > 0 ))
              {
                cout << "Reclaim and Rebalance -- Empty bin" << endl;
                cout << "Parent=";
                lParent->OutputInfo();
                cout << "Child=";
                lChild->OutputInfo();
                assert(0);
              }
            }
          }
                    // release it
          _sets[ lChildOffset ] = NULL;
          _released_set_offsets.push_back( lChildOffset );
          
          if ( !(lParent->AmValid() && lChild->AmValid()) )
          {
            cout << "Reclaim and Rebalance -- Failed to accept foster child" << endl;
            cout << "Parent=";
            lParent->OutputInfo();
            cout << "Child=";
            lChild->OutputInfo();
            assert ( 0 );
          }
        }
                
        if ( lParent->Fosters.size() >= MAX_FOSTERS_ADD )
          lPotentialParents.pop_back(); // this guy's filled up
      }    
      
      cout << "Fostered " << lFosterChildCount - lPotentialFosters.size() << " sets to " << lFosterParentCount - lPotentialParents.size() << "-ish parents" << endl;

      canonicalize();
    }
    
    void rehome_too_big_sets()
    {

      if ( _released_set_offsets.size() > 0 ) // no point in doing anything if we have no space to go to.
      {
        // now, go through and release those fosters that have gotten too big
        vector<SetHandle> lSetsThatNeedReleaseBecauseParentIsTooBig;
        vector<SetHandle> lSetsThatNeedReleaseBecauseParentIsTooCrowded;
        vector<SetHandle> lSetsThatNeedReleaseBecauseTheyAreTooBig;
        for ( int o = 1; o < SETS_SIZE; ++o )
       
         {
          if ( _sets[o] != NULL )
          {
            if ( (*(_sets[o]))->Fosters.size() > 0 ) // evaluate the children
            {
              SetHandle lParent = *_sets[o];
              
              if ( !(lParent->AmValid()) )
              {
                cout << "Rehome_too_big_set -- parent invalid" << endl;
                cout << "Parent=";
                lParent->OutputInfo();
                assert ( 0 );
              }

              bool lTooCrowded = false;
              if ( lParent->Fosters.size() > TOO_MANY_FOSTERS_COUNT )
                lTooCrowded = true;
                
              bool lTooBig = false;
              if ( !lParent->AmStoringHashes() )
                lTooBig = true;
              
              for ( int p = 0; p < lParent->Fosters.size(); ++p )
              {
                if ( lParent->Fosters[p]->ShouldStopStoringHashes() )
                  lSetsThatNeedReleaseBecauseTheyAreTooBig.push_back( lParent->Fosters[p] );
                else if ( lTooBig )
                  lSetsThatNeedReleaseBecauseParentIsTooBig.push_back( lParent->Fosters[p] );
                else if ( lTooCrowded )
                  lSetsThatNeedReleaseBecauseParentIsTooCrowded.push_back( lParent->Fosters[p] );
              }            
            }
          }
        }
        
        // now, do the "too big" sets
        int lPromoteCount = 0;
        while ( _released_set_offsets.size() > 0 && lSetsThatNeedReleaseBecauseTheyAreTooBig.size() > 0 ) 
        {
          SetHandle lSet = lSetsThatNeedReleaseBecauseTheyAreTooBig.back();
          lSetsThatNeedReleaseBecauseTheyAreTooBig.pop_back();
          
          promote_set( lSet );
          
          ++lPromoteCount;
        }
        cout << "Promoted " << lPromoteCount << " out of " << lPromoteCount + lSetsThatNeedReleaseBecauseTheyAreTooBig.size() << " too big sets " << endl;

        // if there's space, do the parent is too big sets
        lPromoteCount = 0;
        while ( _released_set_offsets.size() > 0 && lSetsThatNeedReleaseBecauseParentIsTooBig.size() > 0 ) 
        {
          SetHandle lSet = lSetsThatNeedReleaseBecauseParentIsTooBig.back();
          lSetsThatNeedReleaseBecauseParentIsTooBig.pop_back();
          
          promote_set( lSet );
          ++lPromoteCount;

        }
        cout << "Promoted " << lPromoteCount << " out of " << lPromoteCount + lSetsThatNeedReleaseBecauseParentIsTooBig.size() << " parent too big sets " << endl;

        
        // if there's space, do the parent is too crowded sets
        lPromoteCount = 0;
        while ( _released_set_offsets.size() > 0 && lSetsThatNeedReleaseBecauseParentIsTooCrowded.size() > 0 ) 
        {
          SetHandle lSet = lSetsThatNeedReleaseBecauseParentIsTooCrowded.back();
          lSetsThatNeedReleaseBecauseParentIsTooCrowded.pop_back();
          
          promote_set( lSet );
          ++lPromoteCount;

        }
        cout << "Promoted " << lPromoteCount << " out of " << lPromoteCount + lSetsThatNeedReleaseBecauseParentIsTooCrowded.size() << " parent too crowded sets " << endl;        
      }
    }
    
    void promote_set( SetHandle aSet )
    {
      SetOffset lAddress = get_free_address();
      
      if (!( lAddress > 0 ))// this should work      
      {
        cout << "promote set -- address <= 0" << endl;
        assert(0);
      }
    
      if (!( aSet->Parent->EmancipateFosterChild( aSet, lAddress ) ))// this should work      
      {
        cout << "promote set -- failed to emancipate foster child" << endl;
        cout << "aSet=";
        aSet->OutputInfo();
        cout << "Parent=";
        aSet->Parent->OutputInfo();

        assert(0);
      }
      _sets[ lAddress ] = aSet->Self;
      
      // now, point all the kid's hashes to the new offset.
      for ( int n = 0; n < HASHES; ++n )
      {
        for ( set<HashIntoType>::iterator lIt = aSet->Hashes.begin(); lIt != aSet->Hashes.end(); ++lIt )
        {
          SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBinCached( *lIt, n), n);
          _set_offsets[n][ lBin ] = lAddress;
        }
      }
      
      if ( aSet->ShouldStopStoringHashes() )
        aSet->StopStoringHashes();  
      
      if ( !(aSet->AmValid()) )
      {
        cout << "promote_set aSet is invalid" << endl;
        cout << "Set=";
        aSet->OutputInfo();
        assert ( 0 );
      }
    }
    
    // add a hash to an existing set
    void add_to_set( SetHandle aSet, HashIntoType aHash )
    {      
      for ( int i = 0; i < HASHES; ++i ) // go through and make this hash and its bins and set offsets point at this set
      {
        SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBinCached(aHash, i), i );    
        _set_offsets[i][ lBin ] = aSet->GetPrimaryOffset();        
               
        if (!( _set_offsets[i][ lBin ] > 0 ))// somehow we fucked this up.
        {
          cout << "Add_to_set -- bin is empty" << endl;
          cout << "set=";
          aSet->OutputInfo();
          assert(0);
        }
      }
      aSet->AddToSet(aHash);
      aSet->Increment();
      
      if ( !aSet->AmFosterChild() && aSet->ShouldStopStoringHashes() )
      {
        if (!( aSet->StopStoringHashes() ))
        {
          cout << "Add_to_set -- set couldn't stop storing hashes" << endl;
          cout << "set=";
          aSet->OutputInfo();
          assert(0);
        }
      }
      
      //assert( aSet->AmValid() );
      if (!( aSet->AmValid() ))
      {
        cout << "Add_to_set -- aSet isn't valid" << endl;
        cout << "set=";
        aSet->OutputInfo();
        assert(0);
      }
    }
    
    //
    // memory allocation functions
    //
    void populate_hash_table_bit_count_lookups()
    {
      cout << "populate_hash_table_bit_count_lookups" << endl;
      for (int i = 0; i < HASHES; ++i )
      {
        _hash_table_total_bit_counts[i] = _hash_table[i]->CountBits();
        
        unsigned long long lLookupTableSize = (_tablesizes[i] / BIT_COUNT_PARTITION)+1;
        for (unsigned long long j = 0; j < lLookupTableSize; ++j)
        {      
          unsigned long long lSectionStartIndex = j * BIT_COUNT_PARTITION;
          unsigned long long lSectionStopIndex = ((j + 1) * BIT_COUNT_PARTITION) - 1;
          if ( lSectionStopIndex >= _tablesizes[i] )
            lSectionStopIndex = _tablesizes[i] - 1;
          
          // temporarily store the single section count. (this section may be problematic, since ..lookup[i]'s value is a pointer, so what does [j] do?
          _hash_table_bit_counts_lookup[i][j] = _hash_table[i]->CountBits(lSectionStartIndex, lSectionStopIndex);
          
          if ( j > 0 ) // apply the summation
          {
            _hash_table_bit_counts_lookup[i][j] += _hash_table_bit_counts_lookup[i][j-1];
          }
        }
        cout << i << ": " << _hash_table_total_bit_counts[i] << " -- " << ((double)_hash_table_total_bit_counts[i] / (double)_tablesizes[i]) * 100 << "% occupancy" << endl;
      }
    }
    
    void deallocate_hash_table_preliminary()
    {
      for (int i = 0; i < HASHES; ++i)
      {
        delete _hash_table_preliminary[i];
      }      
    }
    
    void allocate_set_offset_table()
    {
      for (int i = 0; i < HASHES; ++i)
      {
        _set_offsets[i] = new unsigned short[ _hash_table_total_bit_counts[i] ];
        memset(_set_offsets[i], 0, _hash_table_total_bit_counts[i] * sizeof(unsigned short));
      }
    }

    // create a set
    // ideally, create each new set at its own address
    // otherwise, foster the sets out
    SetHandle create_set() 
    {
      SetOffset lAddress = get_free_address();
      SetHandle lSet = NULL;
      
      if ( lAddress == 0 )
      {
        // so, we don't have any free addresses. start a round of fostering
        SetHandle lParentSet = get_least_crowded_set();
        
        if ( lParentSet == NULL ) // we're out of potential foster parents...
        {
          // slim down the list, free up some numbers, also de-foster guys that should be free.
          canonicalize();
          
          if ( _released_set_offsets.empty() )
            reclaim_and_rebalance();
          
          if ( !_released_set_offsets.empty() )
            rehome_too_big_sets();
          
          // re-populate the potential foster parents.
          re_sort_sets();
          
          // try one more time for a free address
          lAddress = get_free_address();
          if ( lAddress == 0 ) // if there are no free addresses, foster
          {
            lSet = new CanonicalSet(); // empty foster set.
            //assert(lParentSet->AcceptFosterChild( lSet )); // if this doesn't work, there's smething very wrong
            if (!( lParentSet->AcceptFosterChild( lSet ) ))
            {
              cout << "create_set -- parent set doesn't accept foster child" << endl;
              cout << "lParentSet=";
              lParentSet->OutputInfo();
              cout << "lSet=";
              lSet->OutputInfo();
              assert(0);
            }
            
          }
          else // woohoo, no need to start the fostering round yet.
          {
            lSet = new CanonicalSet( lAddress );
            _sets[ lAddress ] = lSet->Self;
          }
        } 
        else // here's one
        {
          lSet = new CanonicalSet(); // empty foster set.
          //assert (lParentSet->AcceptFosterChild( lSet )); //if this doesn't work, there's something very wrong.
          if (!( lParentSet->AcceptFosterChild( lSet ) ))
          {
            cout << "create_set -- parent set doesn't accept foster child" << endl;
            cout << "lParentSet=";
            lParentSet->OutputInfo();
            cout << "lSet=";
            lSet->OutputInfo();
            assert(0);
          }
        }
      }
      else
      {
        lSet = new CanonicalSet( lAddress );
        _sets[ lAddress ] = lSet->Self;
      }
      
      //assert( lSet->AmValid() );
      if (!( lSet->AmValid() ))
      {
        cout << "create_set -- lSet not valid" << endl;
        cout << "lSet=";
        lSet->OutputInfo();
        assert(0);
      }
      
      return lSet;
    }  
    
    void join ( SetHandle aJoinee, SetHandle aJoiner )
    {
      if ( aJoinee->AmFosterChild() )      
      {
//        assert( aJoiner->AmFosterChild() );
        
        if (!( aJoiner->AmFosterChild() ))
        {
          cout << "join -- joinee is foster child, but joiner isn't" << endl;
          cout << "aJoinee=";
          aJoinee->OutputInfo();
          cout << "aJoiner=";
          aJoiner->OutputInfo();
          assert(0);
        }
      }
      
//      assert ( aJoiner->AmValid() );
      if (!( aJoiner->AmValid() ))
      {
        cout << "join -- aJoiner isn't valid" << endl;
        cout << "aJoiner=";
        aJoiner->OutputInfo();
        assert(0);
      }
      
      //assert ( aJoinee->AmValid() );
      if (!( aJoinee->AmValid() ))
      {
        cout << "join -- aJoiner isn't valid" << endl;
        cout << "aJoinee=";
        aJoinee->OutputInfo();
        assert(0);
      }
      
      
      // move the back-references
      for ( int i = 0; i < aJoiner->BackReferences.size(); ++i )
      {
        *(aJoiner->BackReferences[i]) = aJoinee;
        aJoinee->BackReferences.push_back( aJoiner->BackReferences[i] );
      }

      // add up the stuff
      aJoinee->Increment( aJoiner->GetKmerCount() );
      
      if ( aJoiner->AmStoringHashes() ) // bring hashes over in
      {
        for ( set<HashIntoType>::iterator lIt = aJoiner->Hashes.begin(); lIt != aJoiner->Hashes.end(); ++lIt )
        {
          if ( aJoinee->AmStoringHashes() )
            aJoinee->Hashes.insert( *lIt );
          
          if ( aJoiner->AmFosterChild() )
          {
            for ( int i = 0; i < HASHES; ++i ) // go through and make this hash and its bins and set offsets point at this set
            {
              SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBinCached(*lIt, i), i );    
              _set_offsets[i][ lBin ] = aJoinee->GetPrimaryOffset();        
              //assert ( _set_offsets[i][ lBin ] > 0 ); // somehow we fucked this up.
              
              if ( !(_set_offsets[i][ lBin ] > 0) )// somehow we fucked this up.
              {
                cout << "get_existing_buck - set offset is 0" << endl;
                cout << "i=" << i;
                cout << "lBin=" << lBin;
                assert( 0 ); // this should work
              }
            }
          }
        }
      }
      
      if ( !aJoiner->AmStoringHashes() && aJoinee->AmStoringHashes() )
      {
        aJoinee->StopStoringHashes();
      }
      
      if ( aJoiner->AmFosterChild() )
      {
        //assert( aJoiner->Parent->EmancipateFosterChild( aJoiner, aJoinee->GetPrimaryOffset()) );
        
        if ( !(aJoiner->Parent->EmancipateFosterChild( aJoiner, aJoinee->GetPrimaryOffset())) )// somehow we fucked this up.
        {
          cout << "join - couldn't emancipate foster child" << endl;
          cout << "aJoinee=";
          aJoinee->OutputInfo();
          cout << "aJoiner=";
          aJoiner->OutputInfo();
          assert( 0 ); // this should work
        }
      }
      else
      {
        if ( aJoiner->Fosters.size() > 0 ) // foster children need to be moved to the main
        {
          for ( int i = 0; i < aJoiner->Fosters.size(); ++i )
          {
            aJoinee->TakeFosterChild( aJoiner->Fosters[i] ); 
          }
        }
      }
      
      delete aJoiner;
      
      if ( !aJoinee->AmFosterChild() && aJoinee->ShouldStopStoringHashes() )
      {
        //assert( aJoinee->StopStoringHashes() );
        
        if ( !(aJoinee->StopStoringHashes()) )// somehow we fucked this up.
        {
          cout << "join - couldn't stop storing hashes" << endl;
          cout << "aJoinee=";
          aJoinee->OutputInfo();
          assert( 0 ); // this should work
        }
      }
      
//      assert ( aJoinee->AmValid() );      
      if ( !(aJoinee->AmValid()) )// somehow we fucked this up.
      {
        cout << "join - joinee isn't valid" << endl;
        cout << "aJoinee=";
        aJoinee->OutputInfo();
        assert( 0 ); // this should work
      }
      
    }
 
    SetHandle get_least_crowded_set()
    { 
      if (_sorted_sets.empty() )
        return NULL;
        //re_sort_sets();
      
      SetHandle lSmallestSet = *(_sorted_sets.back());
      _sorted_sets.pop_back();
      
      return lSmallestSet;
    }
    
    void re_sort_sets()
    {
      _sorted_sets.clear();
      for (int i = 0; i < SETS_SIZE; ++i )
      {
        if ( _sets[i] != NULL && (*(_sets[i]))->GetPrimaryOffset() == i ) // I'm in my canonical spot
          _sorted_sets.push_back( _sets[i] );
      }
      
      sort( _sorted_sets.begin(), _sorted_sets.end(), CanonicalSet::CompSet() );
    }

    SetOffset get_free_address()
    {
      if ( !_released_set_offsets.empty() ) // we've got some released ones to go with.
        return get_a_released_offset(); 
      else if ( _last_set_offset < SETS_SIZE ) // no released ones, but we still have room at the head of the list
      {
        return ++_last_set_offset;
      }
      else
      {
        return 0; // we're fucked. gotta start joining sets
      }
    }
    SetOffset get_a_released_offset()
    {
      SetOffset lSet = _released_set_offsets.back();
      _released_set_offsets.pop_back();
      
      return lSet;
    }
    
    //
    // helper functions
    //
    
    // calculate the hash bin from a kmer's hash
    // store the modulo output for future reference.
    HashBin HashToHashBinCached( HashIntoType aHash, int i ) // no idea if this will be faster than moduloing every time.
    {
      for (int j = 0, index = _hash_bin_cache_last_used_index[i]; j < CACHESIZE; ++j, ++index)
      {
        if ( index >= CACHESIZE )
          index = 0;
        
        if ( _hash_bin_cache[i][index].first == aHash )
        {
          _hash_bin_cache_last_used_index[i] = index;
          return _hash_bin_cache[i][index].second; // found it
        }
      }
      
      // didn't find it.
      HashBin lHashBin = HashToHashBin(aHash, i);
      
      // insert it into the cache;
      _hash_bin_cache_last_used_index[i]++;
      if ( _hash_bin_cache_last_used_index[i] >= CACHESIZE )
        _hash_bin_cache_last_used_index[i] = 0;
      
      _hash_bin_cache[i][ _hash_bin_cache_last_used_index[i] ] = pair<HashIntoType,HashBin>(aHash, lHashBin);
      
      return lHashBin;  
    }    
    HashBin HashToHashBin( HashIntoType aHash, int i ) // no idea which is faster.
    {
      return (aHash % _tablesizes[i]);
    }
    
    // calculate the set offset bin from a hash bin
    // store the value for future reference.
    SetOffsetBin HashBinToSetOffsetBinCached( HashBin aBin, int i ) // no idea if this will be faster
    {
      for (int j = 0, index = _set_offset_bin_cache_last_used_index[i]; j < CACHESIZE; ++j, ++index)
      {
        if ( index >= CACHESIZE )
          index = 0;
        
        if ( _set_offset_bin_cache[i][index].first == aBin )
        {
          _set_offset_bin_cache_last_used_index[i] = index;
          return _set_offset_bin_cache[i][index].second; // found it
        }
      }
      
      // didn't find it.
      HashBin lSetOffsetBin = HashBinToSetOffsetBin(aBin, i);
      
      // insert it into the cache;
      _set_offset_bin_cache_last_used_index[i]++;
      if ( _set_offset_bin_cache_last_used_index[i] >= CACHESIZE )
        _set_offset_bin_cache_last_used_index[i] = 0;
      
      _set_offset_bin_cache[i][ _set_offset_bin_cache_last_used_index[i] ] = pair<HashBin,SetOffsetBin>(aBin, lSetOffsetBin);
      
      return lSetOffsetBin;  
    }    
    SetOffsetBin HashBinToSetOffsetBin( HashBin aBin, int i )
    {      
      //assert( aBin < _tablesizes[i] ); // make sure it's a valid bin.
      if ( !(aBin < _tablesizes[i]) )// make sure it's a valid bin.
      {
        cout << "HashBinToSetOffsetBin - invalid bin" << endl;
        cout << "aBin=" << aBin << endl;
        assert( 0 ); // this should work
      }
      
//      assert( _hash_table[i]->Get( aBin ) == true );
      if ( !(_hash_table[i]->Get( aBin ) == true) )// make sure we have a set
      {
        cout << "HashBinToSetOffsetBin - has a set is false" << endl;
        cout << "aBin=" << aBin << endl;
        assert( 0 ); // this should work
      }
      
      
      unsigned long long lBinSectionIndex = (aBin / BIT_COUNT_PARTITION); // the index of the section before the one we're in
      
//      assert( lBinSectionIndex <= (_tablesizes[i] / BIT_COUNT_PARTITION));       
      if ( !(lBinSectionIndex <= (_tablesizes[i] / BIT_COUNT_PARTITION)) )// make sure we have a set
      {
        cout << "HashBinToSetOffsetBin - lBinSection index is too large" << endl;
        cout << "lBinSectionIndex=" << lBinSectionIndex << endl;
        assert( 0 ); // this should work
      }
      
      unsigned long long lSetOffsetBin = _hash_table[i]->CountBits( lBinSectionIndex * BIT_COUNT_PARTITION, aBin );
      
      if ( lBinSectionIndex > 0 )      
        lSetOffsetBin += _hash_table_bit_counts_lookup[i][ lBinSectionIndex - 1 ];
      
//      assert( lSetOffsetBin > 0 );      
      if ( !(lSetOffsetBin > 0) )// make sure we have a set
      {
        cout << "HashBinToSetOffsetBin - lSetOffsetBin == 0" << endl;
        cout << "lSetOffsetBin=" << lSetOffsetBin << endl;
        assert( 0 ); // this should work
      }
      
      lSetOffsetBin -= 1; // to make it an array index rather than a count.
      
      return lSetOffsetBin;
    }
    
    SetHandle SetOffsetToSet( SetOffset aOffset )
    {      
//      assert( aOffset > 0 ); // we should not be poking here if we don't have a set to go to
      if ( !(aOffset > 0) )// // we should not be poking here if we don't have a set to go to
      {
        cout << "SetOffsetToSet - aOffset == 0" << endl;
        cout << "aOffset=" << aOffset << endl;
        assert( 0 ); // this should work
      }
      
      return *_sets[ aOffset ]; // pick up the gateway node      
    }
    
    bool is_prime( unsigned long long aCandidate )
    {
      if ( aCandidate < 2 )
        return false;
      
      if ( aCandidate == 2 )
        return true;
      
      if ( aCandidate % 2 == 0 )
        return false;
      
      for ( unsigned long long i = 3; i < pow((double)aCandidate, 0.5) + 1; i+= 2 )
      {
        if ( aCandidate % i == 0 )
          return false;
      }
      
      return true;  
    }
    
    unsigned long long get_first_prime_below( unsigned long long aNumber )
    {
      unsigned long long i = aNumber - 1;
      
      if ( i % 2 == 0 ) // no even
        --i;
      
      while ( i > 0 )
      {
        if ( is_prime( i ) )
          return i;
        
        i -= 2;
      }
      
      return 0;
    }
    
    unsigned long long get_first_prime_above( unsigned long long aNumber )
    {
      unsigned long long i = aNumber + 1;
      
      if ( i % 2 == 0 ) // no even
        ++i;
      
      while ( true )
      {
        if ( is_prime( i ) )
          return i;
        
        i += 2;
      }
    }
  };
}
    