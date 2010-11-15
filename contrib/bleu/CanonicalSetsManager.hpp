/*
 *  CanonicalSetManager.h
 *  bleu
 *
 *  Created by Rosangela Canino-Koning on 10/8/10.
 *
 */

#include "../../lib/hashtable.hh"
#include "../../lib/parsers.hh"
#include "../external_lib/cBitArray.h"
#include "CanonicalSet.hpp"
#include "SequenceHash.hpp"
#include <algorithm>

#define BIT_COUNT_PARTITION 1000

#define HASHES 8
#define CACHESIZE 10
#define SET_OFFSET_BITS 20
#define SETS_SIZE pow(2, SET_OFFSET_BITS)
#define CANONICALIZATION_THRESHOLD SETS_SIZE * .05

// 8mb cache, 2gb footprint
//#define CACHED_HASH_SEGMENT_SIZE 1000000
//#define CACHED_HASH_SEGMENT_CAPACITY 10000

//// 4mb cache, 2gb footprint
#define CACHED_HASH_SEGMENT_SIZE 500000
#define CACHED_HASH_SEGMENT_CAPACITY 5000

//// 2mb cache, 2gb footprint
//#define CACHED_HASH_SEGMENT_SIZE 250000
//#define CACHED_HASH_SEGMENT_CAPACITY 2500

//// 2mb cache, 4gb footprint
//#define CACHED_HASH_SEGMENT_SIZE 250000
//#define CACHED_HASH_SEGMENT_CAPACITY 5000

// 2mb cache, 8gb footprint
//#define CACHED_HASH_SEGMENT_SIZE 250000
//#define CACHED_HASH_SEGMENT_CAPACITY 10000

//// 1mb cache, 2gb footprint
//#define CACHED_HASH_SEGMENT_SIZE 125000
//#define CACHED_HASH_SEGMENT_CAPACITY 1250

// 512k cache, 2gb footprint
//#define CACHED_HASH_SEGMENT_SIZE 62500
//#define CACHED_HASH_SEGMENT_CAPACITY 625

//// 512k cache, 4gb footprint
//#define CACHED_HASH_SEGMENT_SIZE 62500
//#define CACHED_HASH_SEGMENT_CAPACITY 1250

// 256k cache, 4gb footprint
//#define CACHED_HASH_SEGMENT_SIZE 31250
//#define CACHED_HASH_SEGMENT_CAPACITY 1024






namespace bleu {
  
  using namespace khmer;
  using namespace std;
  
  typedef CanonicalSet * SetHandle;
  typedef SetHandle * SetPointer;
  
  typedef unsigned long long HashBin;
  typedef unsigned long long SetOffsetBin;
  typedef unsigned long long SetOffset;
  
  typedef unsigned long long SetOffsetContainer;
  
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
    SetOffsetContainer * _set_offsets[HASHES]; 
    

    
    
    // caching system so I don't have to re-modulo. Dunno if this is faster than modulo or not. will test.
    pair<SequenceHash, HashBin> _hash_bin_cache[HASHES][CACHESIZE];
    int _hash_bin_cache_last_used_index[HASHES];
    
    // caching system so I don't have to keep re-counting the bits. Dunno if this is faster than counting or not.
    pair<HashBin, SetOffsetBin> _set_offset_bin_cache[HASHES][CACHESIZE];
    int _set_offset_bin_cache_last_used_index[HASHES];
    
    vector<SetOffset> _released_set_offsets;
    unsigned int _releasable_set_offsets_count;
    
    vector<SetPointer> _sorted_sets;
    
    HashBin ** SeenHashes[HASHES];
    unsigned int * SeenHashCounts[HASHES];
  public:
    unsigned long long _average_size_at_sort;
    
    // contents set during (based on the processing of the reads)  
    vector<SetPointer> _sets; // array of sets
    
  public:
    CanonicalSetsManager( unsigned long long aMaxMemory )
    {
      _last_set_offset = 0; // init to zero
      _average_size_at_sort = 0;
      _releasable_set_offsets_count = 0;
      
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
        
        SeenHashes[j] = new HashBin * [(_tablesizes[j] / CACHED_HASH_SEGMENT_SIZE)];
        SeenHashCounts[j] = new unsigned int [(_tablesizes[j] / CACHED_HASH_SEGMENT_SIZE)];
        
        for ( int k = 0; k < (_tablesizes[j] / CACHED_HASH_SEGMENT_SIZE) + 1; ++k )
        {
          SeenHashes[j][k] = new HashBin[ CACHED_HASH_SEGMENT_CAPACITY ];
          SeenHashCounts[j][k] = 0;
        }
      }
      
      
      _sets.resize(SETS_SIZE, NULL);
      cout << SETS_SIZE << endl;
    }
    //
    // first pass through the reads -- determine if reads are interesting 
    //   enough to warrant getting space to call out a set
    //    
    
    // mark a kmer/hash as "interesting" if it appears more than once.
//    void seen_hash( HashBin * aBins, unsigned long long lCount, int lTable )
//    {
//      cBitArray * lPrelim = _hash_table_preliminary[lTable];
//      cBitArray * lFinal = _hash_table[lTable];
//      HashBin lBin = 0;
//      for ( int i = 0; i < lCount; ++i )
//      {
//        lBin = aBins[i];
//        
//        if ( lPrelim->Get(lBin) )
//          lFinal->Set(lBin, true);
//        else
//          lPrelim->Set(lBin, true);
//
//      }
//    }
    
    void seen_hash( SequenceHash aHash )
    {
      for ( int i = 0; i < HASHES; ++i )
      {
        unsigned long long lHashBin = HashToHashBin(aHash, i); 
        unsigned long long lSegmentIndex = lHashBin / CACHED_HASH_SEGMENT_SIZE;
        
        SeenHashes[i][ lSegmentIndex ][ SeenHashCounts[i][lSegmentIndex]++ ] = lHashBin;
        
        if ( SeenHashCounts[i][ lSegmentIndex ] == CACHED_HASH_SEGMENT_CAPACITY )
        {
          //cout << "dumping: " << i << " segment " << lHashBin / CACHED_HASH_SEGMENT_SIZE << endl;
          dump_seen_hashes_into_table( SeenHashes[i][ lSegmentIndex ], SeenHashCounts[i][ lSegmentIndex ], i );
          SeenHashCounts[i][ lSegmentIndex ]=0;
          //cout << " done" << endl;
        }
      }      
    }
    
    void dump_seen_hashes_into_table( HashBin * aHashBins, unsigned int aCount, int i )
    {
      //cout << ".";
      for ( int j = 0; j < aCount; ++j )
      {
        if ( _hash_table_preliminary[i]->Get( aHashBins[j] ) == true )
          _hash_table[i]->Set( aHashBins[j], true );
        else
          _hash_table_preliminary[i]->Set( aHashBins[j], true );
      }
    }
    
    void finalize_seen_hash()
    {
      for ( int i = 0; i < HASHES; ++i )
      {
        for ( int j = 0; j < (_tablesizes[i] / CACHED_HASH_SEGMENT_SIZE) + 1; ++j )
        {
          dump_seen_hashes_into_table( SeenHashes[i][j], SeenHashCounts[i][j], i );
          delete SeenHashes[i][j];
        }
        //delete SeenHashes[i];
        //delete SeenHashCounts[i];                          
      }
    }
    
//    void seen_hash( HashIntoType aHash )
//    {
//      
//      for ( int i = 0; i < HASHES; ++i )
//      {
//        unsigned long long lHashBin = HashToHashBin(aHash, i); 
//
//        if ( _hash_table_preliminary[i]->Get(lHashBin) == true )
//          _hash_table[i]->Set(lHashBin, true);
//        else
//          _hash_table_preliminary[i]->Set(lHashBin, true);
//       } 
//          
//    }
    
    //
    // second pass through the reads -- actually perform the set assignments
    //
    
    // check whether a kmer can have a set. (I keep wanting to type can_has_set) Damned lolcats.
    bool can_have_set( SequenceHash aHash )
    {
      bool lCanHaveASet = false;
      for ( int i = 0; i < HASHES; ++i )
      {
        HashBin lHashBin = HashToHashBin(aHash, i);
        
        if ( i == 0 )
          lCanHaveASet = _hash_table[i]->Get(lHashBin) == true;
        else if ( _hash_table[i]->Get(lHashBin) != lCanHaveASet ) // we have a disagreement.
          return false;
      }

      return lCanHaveASet;
    }    
    // find a set for this hash. if there isn't one, create it.
    SetHandle get_set( SequenceHash aHash )
    { 
      SetHandle lHandle = NULL;
      
      if ( has_existing_set( aHash ) )
        lHandle = get_existing_set( aHash );
      else
      {
        lHandle = create_set_canonicalize_first();
        add_to_set(lHandle, aHash);
      }
      
      return lHandle;
    }
    // check whether a kmer has already been assigned a set.
    bool has_existing_set( SequenceHash aHash )
    {
      for ( int i = 0; i < HASHES; ++i )
      {
        SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBin(aHash, i), i );
        if ( SetOffsetBinToSetOffset( lBin, i ) == 0 )
          return false;
      }
      return true;
    }
    // find a set for this hash that already exists.
    SetHandle get_existing_set( SequenceHash aHash, SetHandle aWorkingSet=NULL )
    {    
      map<SetHandle, int> lRepresented;
      
      for ( int i = 0; i < HASHES; ++i )
      {
        SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBin(aHash, i), i );
        SetOffset lOffset = SetOffsetBinToSetOffset(lBin, i);
        assert(lOffset > 0);        
        lRepresented[ SetOffsetToSet( lOffset ) ]++;
      }
      
      
      // figure out what the consensus set was
      SetHandle lMostRepresented = NULL; 
      for ( map<SetHandle, int>::iterator lSet = lRepresented.begin(); lSet != lRepresented.end(); ++lSet ) // count the votes
      {
        if ( lMostRepresented == NULL || lSet->second > lRepresented[ lMostRepresented ] ) // because maps are sorted, this will end up being the set with the lowest pointer, and the highest count.
          lMostRepresented = lSet->first;
      }
      
      // even if there isn't consensus, well, we just join up the working set with the one who won the election. Good enough.
      return lMostRepresented;      
      
    }

    
    // add a hash to an existing set
    void add_to_set( SetHandle aSet, SequenceHash aHash )
    {   
      for ( int i = 0; i < HASHES; ++i ) // go through and make this hash and its bins and set offsets point at this set
      {
        SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBin(aHash, i), i );    
        AssignSetOffsetToBin( lBin, i, aSet->GetPrimarySetOffset(), false ); // don't overwrite
      }
      aSet->KmerCount++;
    }

    // create a set
    SetHandle create_set_foster_first()
    {
      SetOffset lAddress = get_free_address();
      SetHandle lSet = NULL;
      
      if ( lAddress == 0 ) // damn. do a round of fostering
      {
        lSet = get_least_crowded_set();
        if ( lSet != NULL ) 
        {
          lSet->JoinOfConvenience = true;
        }
        else // damn damn.
        {
          // does it make sense to canonicalize? payback would be more than 5% of total
          if ( _releasable_set_offsets_count > CANONICALIZATION_THRESHOLD )
          {
            canonicalize();
            // one more try
            lAddress = get_free_address();
          }
          else 
          {
            re_sort_sets();
            lSet = get_least_crowded_set();
            lSet->JoinOfConvenience = true;
          }
        }
      }

      if ( lAddress != 0 )
      {
        lSet = new CanonicalSet( lAddress );
        _sets[ lAddress ] = lSet->Self;  
      }
      
      return lSet;
    }
    
    // create a set
    // try for a free address
    // if that fails, try for canonicalization
    // if that fails, attempt a canonicalization
    // if I can't do that, repopulate the set lists, and fetch least crowded and go with that.
    
    SetHandle create_set_canonicalize_first()
    {
      SetOffset lAddress = get_free_address();
      SetHandle lSet = NULL;
      
      if ( lAddress == 0 ) // damn
      {        
        if ( _releasable_set_offsets_count > CANONICALIZATION_THRESHOLD ) // try to canonicalize first
        {
          _sorted_sets.clear(); // reset here
          canonicalize();
          // one more try
          lAddress = get_free_address();
        } 
        else // foster for this one.
        {
          lSet = get_least_crowded_set();
          if ( lSet != NULL ) // did it work?
          {
            lSet->JoinOfConvenience = true;
          }
          else // resort and try again.
          {
            re_sort_sets();
            lSet = get_least_crowded_set();
            lSet->JoinOfConvenience = true;
          }
        }
      }
      
      if ( lAddress != 0 )
      {
        lSet = new CanonicalSet( lAddress );
        _sets[ lAddress ] = lSet->Self;  
      }
      
      return lSet;
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
      assert( aOriginatingSet != aEncounteredSet); // we really really shouldn't be the same set.

      if ( aOriginatingSet->GetPrimarySetOffset() < aEncounteredSet->GetPrimarySetOffset() ) // the bigger one wins, since it almost certainly has more back references.
      {
        join( aOriginatingSet, aEncounteredSet );
        return aOriginatingSet;
      }
      else
      {
        join( aEncounteredSet, aOriginatingSet );
        return aEncounteredSet;
      }
      //return lDominatingSet;
    }
    
    void canonicalize()
    {   
      time_t lStart;
      lStart = time(NULL);
      cout << "Canonicalization starting." << endl;
      cout << "expect to release " << _releasable_set_offsets_count << " offsets." << endl;
      
      if ( _releasable_set_offsets_count <= CANONICALIZATION_THRESHOLD )
      {
        cout << "Canonicalization benefit expectation is less than " << CANONICALIZATION_THRESHOLD << " freed offsets. Not advised." << endl; 
      }


      // go through all the set offsets and point them to their canonical locations.
      for ( int i = 0; i < HASHES; ++i )
      {
        for ( int j = 0; j < _hash_table_total_bit_counts[i]; ++j )
        {
          SetOffset lOffset = SetOffsetBinToSetOffset(j, i);
          if ( lOffset != 0 ) // there's something in here
          {
            AssignSetOffsetToBin( j, i, SetOffsetToSet( lOffset )->GetPrimarySetOffset(), true );
//            _set_offsets[i][j] = SetOffsetToSet( lOffset )->PrimarySetOffset; 
          }
        }
      }
      
      _released_set_offsets.clear(); // clear this.
      
      for ( int k = 1; k < SETS_SIZE; ++k )
      {
        if ( _sets[ k ] != NULL ) // something in this slot
        {
          SetHandle lSet = *_sets[k];
          
          if ( k != lSet->GetPrimarySetOffset() )
          {
            _sets[k] = NULL;
            _released_set_offsets.push_back( k ); // it's null! hey!
          }
          else // we'll only hit each set once this way.
            lSet->BackReferences.clear();

        }
        else
        {
          _released_set_offsets.push_back( k ); // it's null! hey!
        }
      }
      
      cout << "Canonicalizing: Released " << _released_set_offsets.size() << " sets, in " << difftime(time(NULL), lStart) << " seconds." << endl;
      _releasable_set_offsets_count = 0;
    }
    
    void join ( SetHandle aJoinee, SetHandle aJoiner )
    {
      // move the back-references
      *aJoiner->Self = aJoinee;
      aJoinee->BackReferences.push_back( aJoiner->Self );
      for ( int i = 0; i < aJoiner->BackReferences.size(); ++i )
      {
        *(aJoiner->BackReferences[i]) = aJoinee;
        aJoinee->BackReferences.push_back( aJoiner->BackReferences[i] );
      }
      
      // add up the stuff
      aJoinee->KmerCount += aJoiner->KmerCount;
      delete aJoiner;
      
      _releasable_set_offsets_count++;
    }
 
    SetHandle get_least_crowded_set()
    { 
      if (_sorted_sets.empty() )
        return NULL;
      
      SetHandle lSmallestSet = NULL;
      while ( (lSmallestSet == NULL || lSmallestSet->KmerCount >= _average_size_at_sort ) && !_sorted_sets.empty() )
      { 
        lSmallestSet = *(_sorted_sets.back());

        // if the one we just grabbed was too big, pop it, and grab the next fresh one. We don't care about balance here, just efficiency.
        if ( lSmallestSet != NULL && lSmallestSet->KmerCount >= _average_size_at_sort )
        {
          _sorted_sets.pop_back();
        }
      }      
      
      return lSmallestSet; // it'll either have a legit set, or it'll be empty.
    }
    
    void re_sort_sets()
    {
      _sorted_sets.clear();

      unsigned long long lCount = 0;
      for (int i = 0; i < SETS_SIZE; ++i )
      {
        if ( _sets[i] != NULL && (*(_sets[i]))->GetPrimarySetOffset() == i ) // I'm in my canonical spot
        {
          lCount += (*(_sets[i]))->KmerCount;
          _sorted_sets.push_back( _sets[i] );
        }
      }
      int lSize = _sorted_sets.size();
      _average_size_at_sort = lCount / lSize; // re-set it.

      
      sort( _sorted_sets.begin(), _sorted_sets.end(), CanonicalSet::CompSet() );
      int pos = 0;
      for (; pos < lSize; ++pos) // troll through them until you find the edge of the too-bigs
      {
        if ( (*(_sorted_sets[pos]))->KmerCount < _average_size_at_sort )
          break;
      }
      
      _sorted_sets.erase( _sorted_sets.begin(), _sorted_sets.begin() + pos ); // pull them out.
    }

    SetOffset get_free_address()
    {
      if ( !_released_set_offsets.empty() ) // we've got some released ones to go with.
      {
        return get_a_released_offset(); 
      }
      else if ( _last_set_offset < SETS_SIZE-1 ) // no released ones, but we still have room at the head of the list
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
    // memory allocation functions
    //
    void populate_hash_table_bit_count_lookups()
    {
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
          _hash_table_bit_counts_lookup[i][j] = _hash_table[i]->CountBits2(lSectionStartIndex, lSectionStopIndex);
          
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
        
        unsigned long long lContainerCount = ceil( (SET_OFFSET_BITS/(double)(sizeof(SetOffsetContainer)*8)) *_hash_table_total_bit_counts[i]);
      
        _set_offsets[i] = new SetOffsetContainer[ lContainerCount ];
        memset(_set_offsets[i], 0, lContainerCount * sizeof(SetOffsetContainer));
      }
    }
    
    //
    // helper functions
    //
    
//    // calculate the hash bin from a kmer's hash
//    // store the modulo output for future reference.
//    HashBin HashToHashBinCached( HashIntoType aHash, int i ) // no idea if this will be faster than moduloing every time.
//    {
//      for (int j = 0, index = _hash_bin_cache_last_used_index[i]; j < CACHESIZE; ++j, ++index)
//      {
//        if ( index >= CACHESIZE )
//          index = 0;
//        
//        if ( _hash_bin_cache[i][index].first == aHash )
//        {
//          _hash_bin_cache_last_used_index[i] = index;
//          return _hash_bin_cache[i][index].second; // found it
//        }
//      }
//      
//      // didn't find it.
//      HashBin lHashBin = HashToHashBin(aHash, i);
//      
//      // insert it into the cache;
//      _hash_bin_cache_last_used_index[i]++;
//      if ( _hash_bin_cache_last_used_index[i] >= CACHESIZE )
//        _hash_bin_cache_last_used_index[i] = 0;
//      
//      _hash_bin_cache[i][ _hash_bin_cache_last_used_index[i] ] = pair<HashIntoType,HashBin>(aHash, lHashBin);
//      
//      return lHashBin;  
//    }    
    HashBin HashToHashBin( SequenceHash aHash, int i ) // no idea which is faster.
    {
      return (aHash.canonical_hash % _tablesizes[i]);
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
      assert( aBin < _tablesizes[i] ); // make sure it's a valid bin.
//      assert( _hash_table[i]->Get( aBin ) == true );
      
      unsigned long long lBinSectionIndex = (aBin / BIT_COUNT_PARTITION); // the index of the section before the one we're in
      
      assert( lBinSectionIndex <= (_tablesizes[i] / BIT_COUNT_PARTITION));       
      
      unsigned long long lSetOffsetBin = _hash_table[i]->CountBits2(lBinSectionIndex * BIT_COUNT_PARTITION, aBin );
      
      if ( lBinSectionIndex > 0 )      
        lSetOffsetBin += _hash_table_bit_counts_lookup[i][ lBinSectionIndex - 1 ];
      
      assert( lSetOffsetBin > 0 );      
      
      lSetOffsetBin -= 1; // to make it an array index rather than a count.
      
      return lSetOffsetBin;
    }
    
    SetOffset SetOffsetBinToSetOffset( SetOffsetBin aBin, int i )
    {
      unsigned long long lGlobalBit = aBin * SET_OFFSET_BITS;
      unsigned long long lField = lGlobalBit / (sizeof(SetOffsetContainer)*8);
      unsigned long long lLocalBit = lGlobalBit % (sizeof(SetOffsetContainer)*8);
      
      SetOffset lMask = pow(2, SET_OFFSET_BITS)-1;
      SetOffset lOffset = (_set_offsets[i][lField] >> lLocalBit) & lMask;

      int lSpaceUsed = (sizeof(SetOffsetContainer)*8) - lLocalBit;      
      
//      SetOffsetContainer lDELETEMEField1 = _set_offsets[i][lField];
//      SetOffsetContainer lDELETEMEField2 = 0;
      if (lSpaceUsed < SET_OFFSET_BITS)
      {
        lOffset |= (_set_offsets[i][lField+1] << lSpaceUsed) & lMask;
//        lDELETEMEField2 = _set_offsets[i][lField+1];
      }
      
      assert( lOffset < SETS_SIZE );
      
      return lOffset;      
    }
    
    bool AssignSetOffsetToBin( SetOffsetBin aBin, int i, SetOffset aVal, bool aOverwrite )
    {

      if ( !aOverwrite ) // if overwrite is false
      {
        SetOffset lOldOffset = SetOffsetBinToSetOffset(aBin, i);
        if ( lOldOffset > 0 )
          return false;
      }
    
      SetOffsetContainer lWorkingValue = aVal;
    
      unsigned long long lGlobalBit = aBin * SET_OFFSET_BITS;
      unsigned long long lField = lGlobalBit / (sizeof(SetOffsetContainer)*8);
      unsigned long long lLocalBit = lGlobalBit % (sizeof(SetOffsetContainer)*8);

      int lSpaceUsed = (sizeof(SetOffsetContainer)*8) - lLocalBit;
      SetOffsetContainer lAssignmentMask = ~( ((SetOffsetContainer)pow(2, SET_OFFSET_BITS) - 1) << lLocalBit );
      
//      SetOffsetContainer lDELETEMEField1 = _set_offsets[i][lField];
//      SetOffsetContainer lDELETEMEField2 = 0;
//      SetOffsetContainer lDELETEMEField2After = 0;
//      
//      SetOffsetContainer lDELETEMEShiftedValue = (lWorkingValue << lLocalBit);
//      SetOffsetContainer lDELETEMEMaskedOriginalValue = (_set_offsets[i][lField] & lAssignmentMask);
//      SetOffsetContainer lDELETEMEOREDTOGETHER = (_set_offsets[i][lField] & lAssignmentMask) | (lWorkingValue << lLocalBit);
//
//      SetOffsetContainer lDELETEMEShiftedValue2 = 0;
//      SetOffsetContainer lDELETEMEMaskedOriginalValue2 = 0;
//      SetOffsetContainer lDELETEMEOREDTOGETHER2 = 0;

      
      _set_offsets[i][lField] = (_set_offsets[i][lField] & lAssignmentMask) | (lWorkingValue << lLocalBit);
      
//      SetOffsetContainer lDELETEMEField1After = _set_offsets[i][lField];
      
      if ( lSpaceUsed < SET_OFFSET_BITS )
      {
//        lDELETEMEField2 = _set_offsets[i][lField+1];
        int lShiftCount = lSpaceUsed;
        SetOffsetContainer lShiftMask = ~( ((SetOffsetContainer)pow(2, SET_OFFSET_BITS) - 1) >> lShiftCount );
        
//        lDELETEMEShiftedValue2 = (lWorkingValue >> lShiftCount);
//        lDELETEMEMaskedOriginalValue2 = (_set_offsets[i][lField+1] & lShiftMask );
//        lDELETEMEOREDTOGETHER2 = (_set_offsets[i][lField+1] & lShiftMask ) | (lWorkingValue >> lShiftCount);
        
        _set_offsets[i][lField+1] = (_set_offsets[i][lField+1] & lShiftMask ) | (lWorkingValue >> lShiftCount);
        
//        lDELETEMEField2After = _set_offsets[i][lField+1];
      }
  
      SetOffset lCircVal = SetOffsetBinToSetOffset(aBin, i);
      assert( lCircVal == aVal);
      
      return true;
    }
    
    SetHandle SetOffsetToSet( SetOffset aOffset )
    {      
      assert( aOffset > 0 ); // we should not be poking here if we don't have a set to go to
      return *_sets[ aOffset ]; // pick up the gateway node      
    }
    
    //
    // prime calculation functions
    //
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
    
