///*
// *  SetsManager.hpp
// *  bleu
// *
// *  Created by Rosangela Canino-Koning on 10/4/10.
// *
// */
//#include "../../lib/hashtable.hh"
//#include "../../lib/parsers.hh"
//#include "../external_lib/cBitArray.h"
//#include "Set.hpp"
//
//#define BIT_COUNT_PARTITION 1000
//
//#define HASHES 8
//#define BIN_SIZE 50
//#define FOSTER_COUNT_MULTIPLIER 5
//
//#define CACHESIZE 10
//
//namespace bleu {
//  
//  using namespace khmer;
//  using namespace std;
//  
//  typedef Set * SetHandle;
//  
//  typedef unsigned long long HashBin;
//  typedef unsigned long long SetOffsetBin;
//  typedef unsigned short SetOffset;
//  
//  class ChainedSetsManager
//  {
//  private:
//    // sizes set during constructor
//    cBitArray * _hash_table_preliminary[HASHES]; // two-dimensional hash-table
//    cBitArray * _hash_table[HASHES]; // two-dimensional hash-table
//    
//    unsigned long long _tablesizes[HASHES]; // the sizes of the hash tables
//    unsigned long long * _hash_table_bit_counts_lookup[HASHES]; // the number of bits set in the hash table, by every 10k entries
//    unsigned int _last_set_offset;
//    
//    // sizes set during prep 1 (based on processing of _hash_table)
//    unsigned long long _hash_table_total_bit_counts[HASHES];
//    unsigned short * _set_offsets[HASHES]; 
//    
//    // contents set during (based on the processing of the reads)  
//    vector<SetHandle> _sets; // array of sets
//    
//    map<unsigned long long, list<SetHandle*> > _binned_sets; // our sets, sorted by size, in heap form
//    unsigned long long _smallest_set_count;
//    
//    // caching system so I don't have to re-modulo. Dunno if this is faster than modulo or not. will test.
//    pair<HashIntoType, HashBin> _hash_bin_cache[HASHES][CACHESIZE];
//    int _hash_bin_cache_last_used_index[HASHES];
//    
//    // caching system so I don't have to keep re-counting the bits. Dunno if this is faster than counting or not.
//    pair<HashBin, SetOffsetBin> _set_offset_bin_cache[HASHES][CACHESIZE];
//    int _set_offset_bin_cache_last_used_index[HASHES];
//    
//    vector<SetOffset> _released_set_offsets;
//    
//  public:
//    ChainedSetsManager( unsigned long long aMaxMemory )
//    {
//      _last_set_offset = 0; // init to zero
//      _smallest_set_count = 0;
//      _sets.push_back(NULL); // the 0th entry is invalid.
//      
//      _tablesizes[0] = get_first_prime_below( aMaxMemory / HASHES );
//      for ( int i = 1; i < HASHES; ++i )      
//      {
//        _tablesizes[i] = get_first_prime_below(_tablesizes[i-1]); 
//      }
//      
//      for ( int j = 0; j < HASHES; ++j )      
//      {
//        _hash_table_preliminary[j] = new cBitArray( _tablesizes[j] );
//        _hash_table_preliminary[j]->Clear();      
//        
//        _hash_table[j] = new cBitArray( _tablesizes[j] );
//        _hash_table[j]->Clear(); 
//        
//        _hash_table_bit_counts_lookup[j] = new unsigned long long[(_tablesizes[j] / BIT_COUNT_PARTITION)+1];
//        memset(_hash_table_bit_counts_lookup[j], 0, ((_tablesizes[j] / BIT_COUNT_PARTITION)+1) * sizeof(unsigned long long));
//        
//        // THESE WILL GET SET LATER.            
//        _hash_table_total_bit_counts[j] = 0;
//        _set_offsets[j] = NULL; // a table the size of the number of bits in the hash table.
//        
//        _hash_bin_cache_last_used_index[j] = 0;
//        memset(_hash_bin_cache[j], 0, sizeof(pair<HashIntoType,HashBin>) * CACHESIZE);
//        
//        _set_offset_bin_cache_last_used_index[j] = 0;
//        memset(_set_offset_bin_cache[j], 0, sizeof(pair<HashBin,SetOffsetBin>) * CACHESIZE);
//      }
//    }
//    //
//    // first pass through the reads -- determine if reads are interesting 
//    //   enough to warrant getting space to call out a set
//    //    
//    
//    // mark a kmer/hash as "interesting" if it appears more than once.
//    void seen_hash( HashIntoType aHash )
//    {
//      for (int i = 0; i < HASHES; ++i )
//      {
//        unsigned long long lHashBin = HashToHashBinCached(aHash, i); 
//        
//        if ( _hash_table_preliminary[i]->Get(lHashBin) == true )
//          _hash_table[i]->Set(lHashBin, true);
//        else
//          _hash_table_preliminary[i]->Set(lHashBin, true);
//      }
//    }
//    
//    //
//    // second pass through the reads -- actually perform the set assignments
//    //
//    
//    // check whether a kmer can have a set. (I keep wanting to type can_has_set) Damned lolcats.
//    bool can_have_set( HashIntoType aHash )
//    {
//      bool lCanHaveASet = false;
//      for ( int i = 0; i < HASHES; ++i )
//      {
//        HashBin lHashBin = HashToHashBinCached(aHash, i);
//        
//        if ( i == 0 )
//          lCanHaveASet = _hash_table[i]->Get(lHashBin) == true;
//        else if ( _hash_table[i]->Get(lHashBin) != lCanHaveASet ) // we have a disagreement.
//          return false;
//      }
//      
//      return lCanHaveASet;
//    }
//  
//    // find a set for this hash. if there isn't one, create it.
//    SetHandle get_set( HashIntoType aHash )
//    { 
//      SetHandle lHandle = NULL;
//      
//      if ( has_existing_set( aHash ) )
//        lHandle = get_existing_set( aHash );
//      else
//      {
//        lHandle = create_set();
//        add_to_set(lHandle, aHash);
//      }
//        
//        
//      return lHandle;
//    }
//  
//    
//    // check whether a kmer has already been assigned a set.
//    bool has_existing_set( HashIntoType aHash )
//    {
//      for ( int i = 0; i < HASHES; ++i )
//      {
//        SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBinCached(aHash, i), i );
//                
//        if ( _set_offsets[i][ lBin ] == 0 )
//          return false;
//      }
//      
//      return true;
//    }
//    
//    // find a set for this hash that already exists.
//    SetHandle get_existing_set( HashIntoType aHash )
//    {      
//      map<SetHandle, int> lRepresented;
//      
//      for ( int i = 0; i < HASHES; ++i )
//      {
//        SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBinCached(aHash, i), i );
//        
//        assert ( _set_offsets[i][ lBin ] > 0 ); // somehow we fucked this up.
//        
//        lRepresented[ SetOffsetAndHashToSet( _set_offsets[i][ lBin ], aHash ) ]++;
//      }
//      
//      // figure out what the consensus set was
//      SetHandle lMostRepresented = NULL; 
//      for ( map<SetHandle, int>::iterator lSet = lRepresented.begin(); lSet != lRepresented.end(); ++lSet ) // count the votes
//      {
//        if ( lMostRepresented == NULL || lSet->second > lRepresented[ lMostRepresented ] ) // because maps are sorted, this will end up being the set with the lowest pointer, and the highest count.
//          lMostRepresented = lSet->first;
//      }     
//      
//      return lMostRepresented;      
//    }
//    
//    bool sets_are_disconnected( SetHandle aSet1, SetHandle aSet2 )
//    {
//      if ( aSet1->BiologicalPatriarch == aSet2->BiologicalPatriarch )
//        return false; // they're already in the same set.
//      else
//        return true; // one of them is 
//    }
//    
//    // put together two sets.
//    SetHandle bridge_sets( SetHandle aEncounteredSet, SetHandle aOriginatingSet )
//    {
//      assert ( aOriginatingSet != aEncounteredSet ); // we really really shouldn't be the same set.
//      
//      SetHandle lDominatingSet = NULL;
//      
//      //cout << "originating set: "; aOriginatingSet->output_info();
//      //cout << "encountered set: "; aEncounteredSet->output_info();
//      
//      // now, figure out which one should join to which.
//      // first, if neither of them are fosters, we can just have the larger set adopt the smaller set.
//      if ( !aOriginatingSet->AmFostered() &&
//          !aEncounteredSet->AmFostered() ) // we're not foster children!
//      { 
//        if ( aOriginatingSet->GetTotalKmerCount() > aEncounteredSet->GetTotalKmerCount() )
//        {
//          adopt( aOriginatingSet, aEncounteredSet->BiologicalPatriarch );
//          lDominatingSet = aOriginatingSet;
//        }
//        else
//        {
//          adopt( aEncounteredSet, aOriginatingSet->BiologicalPatriarch );
//          lDominatingSet = aEncounteredSet;
//        }     
//      }
//      else
//      {        
//        // since either one, or both are rooted in fosteredes, we don't join the entire sets, but only the branch of the family rooted up to there.
//        if ( !aOriginatingSet->AmFostered() &&
//               aEncounteredSet->AmFostered() ) // originating set isn't fostered, but encountered is.
//        {
//          // send the encountered set's entire biological tree back to his biological family.
//          adopt( aOriginatingSet, aEncounteredSet->BiologicalPatriarch );        
//          lDominatingSet = aOriginatingSet;
//        }
//        else // the other way around, or they're both fostered, so, makes no difference.
//        {
//          adopt( aEncounteredSet, aOriginatingSet->BiologicalPatriarch ); 
//          lDominatingSet = aEncounteredSet;
//        } 
//      }
//
//      //cout << "dominating set post-bridge: "; lDominatingSet->output_info();
//      
//      return lDominatingSet;
//    }
//    
//    
//    // add a hash to an existing set
//    void add_to_set( SetHandle aSet, HashIntoType aHash )
//    {
//      if ( aSet->StoringHashes )
//        aSet->Hashes.push_back( aHash );
//      
//      aSet->IncrementCount();
//      
//      for ( int i = 0; i < HASHES; ++i ) // go through and make this hash and its bins and set offsets point at this set
//      {
//        SetOffsetBin lBin = HashBinToSetOffsetBinCached( HashToHashBinCached(aHash, i), i );        
//        _set_offsets[i][ lBin ] = aSet->GetPrimaryGatewaySetOffset();        
//        assert ( _set_offsets[i][ lBin ] > 0 ); // somehow we fucked this up.
//      }
//      
//      reevaluate_fostering_order_and_availability( aSet );
//    }
//    
//    //
//    // memory allocation functions
//    //
//    void populate_hash_table_bit_count_lookups()
//    {
//      cout << "populate_hash_table_bit_count_lookups" << endl;
//      for (int i = 0; i < HASHES; ++i )
//      {
//        _hash_table_total_bit_counts[i] = _hash_table[i]->CountBits();
//        
//        unsigned long long lLookupTableSize = (_tablesizes[i] / BIT_COUNT_PARTITION)+1;
//        for (unsigned long long j = 0; j < lLookupTableSize; ++j)
//        {      
//          unsigned long long lSectionStartIndex = j * BIT_COUNT_PARTITION;
//          unsigned long long lSectionStopIndex = ((j + 1) * BIT_COUNT_PARTITION) - 1;
//          if ( lSectionStopIndex >= _tablesizes[i] )
//            lSectionStopIndex = _tablesizes[i] - 1;
//          
//          // temporarily store the single section count. (this section may be problematic, since ..lookup[i]'s value is a pointer, so what does [j] do?
//          _hash_table_bit_counts_lookup[i][j] = _hash_table[i]->CountBits(lSectionStartIndex, lSectionStopIndex);
//          
//          if ( j > 0 ) // apply the summation
//          {
//            _hash_table_bit_counts_lookup[i][j] += _hash_table_bit_counts_lookup[i][j-1];
//          }
//        }
//        cout << i << ": " << _hash_table_total_bit_counts[i] << " -- " << ((double)_hash_table_total_bit_counts[i] / (double)_tablesizes[i]) * 100 << "% occupancy" << endl;
//      }
//    }
//    
//    void deallocate_hash_table_preliminary()
//    {
//      for (int i = 0; i < HASHES; ++i)
//      {
//        delete _hash_table_preliminary[i];
//      }      
//    }
//    
//    void allocate_set_offset_table()
//    {
//      for (int i = 0; i < HASHES; ++i)
//      {
//        _set_offsets[i] = new unsigned short[ _hash_table_total_bit_counts[i] ];
//        memset(_set_offsets[i], 0, _hash_table_total_bit_counts[i] * sizeof(unsigned short));
//      }
//    }
//    
//  private:
//
//    
//    //
//    // functions related to set creation, and binning
//    // 
//    
//    
//    // add every new set to the least crowded available set we can find.
//    SetHandle create_set() // the basic starter set, everyone starts out fostered.
//    {
//      SetHandle lRootSet = get_least_crowded_set();
//      
//      if ( lRootSet == NULL || !lRootSet->AmElgibleForFostering() )      
//        lRootSet = create_allocated_address_set(); 
//            
//      SetHandle lSet = new Set();        
//      foster( lRootSet, lSet );
//      
//      return lSet;
//    }   
//    
//    void output_bins()
//    {
//      for ( map<unsigned long long, list<SetHandle*> >::iterator lIt = _binned_sets.begin(); lIt != _binned_sets.end(); ++lIt )
//      {
//        int lll = lIt->first;
//        int size = lIt->second.size();
//        cout << "(" << lll << " " << size << ")";
//      }
//      cout << endl;
//    }
//    
//    SetHandle get_least_crowded_set()
//    { 
//      cout << "get x: ";
//      output_bins();
//      
//      if ( _binned_sets.begin() == _binned_sets.end() )
//        return NULL;
//      
//      return *((_binned_sets.begin())->second.front());
//    }
//
//    // create a real, honest to goodness addressed set.
//    SetHandle create_allocated_address_set()
//    {
//      SetHandle lSet = NULL;
//      
//      if ( !_released_set_offsets.empty() ) // we've got some released ones to go with.
//      {
//        lSet = new Set( get_a_released_offset() ); 
//        add_to_fostering_elgibility( lSet );
//        _sets[ lSet->GetPrimaryGatewaySetOffset() ] = lSet;         
//      }
//      else if ( _last_set_offset < 65535 ) // no released ones, but we still have room at the head of the list
//      {
//        lSet = new Set( ++_last_set_offset ); 
//        
//        add_to_fostering_elgibility( lSet );
//        _sets.push_back( lSet );  
//      }
//      else
//      {
//        assert(0); // this is pretty fucked. how often does this happen?
//                   // in this case, you have to do some shady jiggery pokery.
//      }
//      
//      return lSet;
//    }
//    SetOffset get_a_released_offset()
//    {
//      SetOffset lSet = _released_set_offsets.back();
//      _released_set_offsets.pop_back();
//      
//      return lSet;
//    }
//    void add_to_fostering_elgibility( SetHandle aSet, bool aTargetGiven = false, unsigned long long aTarget = 0 )
//    {
//      cout << "add b: ";
//      output_bins();
//      
//      if ( !aSet->AmElgibleForFostering() )
//        return;
//
//      unsigned long long lBinTarget = aTarget;
//      
//      if ( !aTargetGiven )
//        lBinTarget = calculate_order_bin( aSet );
//      
//      _binned_sets[lBinTarget].push_front( aSet->Self );
//      
// 
//      aSet->CurrentBin = lBinTarget;      
//      aSet->BinLocation = _binned_sets[ lBinTarget ].begin();
//      aSet->AmInBin = true;
//      
//      cout << "add=" << lBinTarget << " " << aSet->CurrentBin << endl;
//      
//      cout << "add a: ";
//      output_bins();
//
//      
//    }
//    void reevaluate_fostering_order_and_availability( SetHandle aSet )
//    {
//        unsigned long long lBinTarget = calculate_order_bin( aSet );
//        if ( aSet->CurrentBin != lBinTarget )
//        {
//          remove_from_fostering_elgibility( aSet );          
//          add_to_fostering_elgibility( aSet, true, lBinTarget );          
//        }
//
//    }
//    void remove_from_fostering_elgibility( SetHandle aSet )
//    {
//      cout << "remove b: " << aSet->CurrentBin;
//      output_bins();
//
//      if ( aSet->AmInBin )
//      {
//        _binned_sets[ aSet->CurrentBin ].erase( aSet->BinLocation ); // remove it
//        aSet->AmInBin = false;
//      
//        //cout << ">" << _binned_sets[ aSet->CurrentBin ].size() << "<";
//        
//        if ( _binned_sets[ aSet->CurrentBin ].size() == 0 )
//        {
//          cout << "remove eb: ";
//          output_bins();
//          
//          cout << "erasing=" << aSet->CurrentBin << endl;
//          _binned_sets.erase( aSet->CurrentBin );
//          
//          cout << "remove ea: ";
//          output_bins();
//        }
//      }
//        
//      
//      cout << "remove a: ";
//      output_bins();
//
//    }
//    unsigned long long calculate_order_bin( SetHandle aSet ) // do whatever calculus we should do, including counting kmers and how many fosters you already have
//    {
//      return (aSet->GetTotalKmerCount() / BIN_SIZE) + (aSet->GetTotalFosterCount() / FOSTER_COUNT_MULTIPLIER) ;
//    }
////    SetHandle create_set_at_offset( SetOffset aTarget )
////    {
////      SetHandle lSet = new Set( aTarget ); 
////      add_to_order( lSet );
////      _sets[ aTarget ] = lSet; 
////      
////      return lSet;      
////    }
//    
//
//    //
//    // Functions related to fostering and adoption
//    //
//    
//    // a set may be fostered from an unborn state, or from a root state, but never from a child state of another set. Also, I should never be passed around from one foster situation to another.
//    void foster( SetHandle aParent, SetHandle aChild )
//    {
//      assert( aParent != aChild ); // you know?
//      
//      SetHandle lRealNewParent = aParent->BiologicalPatriarch; // always add fosters at the root of a set
//
//      assert( aChild->ParentSet == NULL || aChild->AmFosterRoot() ); // I'm my own root. Srsly.      
//      assert( aChild->StoringHashes ); // can't foster a child that's already lost track of what's in it. fo'sho.
//      
//      lRealNewParent->FosterSets.push_back( aChild );      
//      aChild->ParentSet = lRealNewParent;
//      
//      repoint_to_new_root(lRealNewParent->Root, aChild); // point my root node info and that of all my downstream nodes to this new situation.      
//      repoint_to_new_gateway_node( lRealNewParent->Root, aChild ); // re-point my gateway nodes (assuming I have any) to my new root gateway.
//      
//      // release whatever set offset(s) I was using
//      reclaim_set_offsets( aChild );
//      
//      // make sure that any fosters this new foster has get put in the same line as this one
//      move_fosters_up( aChild );
//      
//      reevaluate_fostering_order_and_availability( lRealNewParent ); // you just added a foster, gotta make sure things are up to date.
//    }    
//    void reclaim_set_offsets( SetHandle aSet )
//    {
//      for (int i = 0; i < aSet->SetOffsetsIOccupy.size(); ++i )
//      {
//        if ( aSet->SetOffsetsIOccupy[i] != 0 ) // don't reclaim blank offsets.
//        {          
//          _released_set_offsets.push_back( aSet->SetOffsetsIOccupy[i] );
//          _sets[ aSet->SetOffsetsIOccupy[i] ] = NULL; // wipe so it doesn't point to me, even accidentally
//        }
//      } 
//      aSet->SetOffsetsIOccupy.clear();
//    }
//    
////    void rebalance_system_if_needed( SetHandle aSet )
////    {
////      assert( aSet->AmBiologicalPatriarch() );
////      collapse_if_needed( aSet );
////      //redistribute_foster_children_if_needed( aSet );
////      reorder_if_needed( aSet );
////    }
//    
//    // a set may be adopted after being root, or being fostered, but never from a child state of another set.
//    void adopt( SetHandle aGatewayNode, SetHandle aChild )
//    {
//      SetHandle lRealNewParent = aGatewayNode->BiologicalPatriarch;
//      assert ( aChild->AmRoot() || aChild->AmFosterRoot() ); // I'm either top root, or the root of a foster branch.
//      if ( aChild->AmFosterRoot() ) 
//      {
//        for ( int i = 0; i < aChild->ParentSet->FosterSets.size(); ++i ) // emancipate the kid so you can adopt it
//        {
//          if ( aChild->ParentSet->FosterSets[i] == aChild )
//          {
//            aChild->ParentSet->FosterSets.erase( aChild->ParentSet->FosterSets.begin() + i );
//            break;
//          }          
//        }
//        
//        if ( aChild->GetPrimaryGatewaySetOffset() != lRealNewParent->GetPrimaryGatewaySetOffset() ) // I'm really moving away, change my address.
//        {
//          // re-point my gateway nodes (assuming I have any) to my new patriarch.
//          repoint_to_new_gateway_node( lRealNewParent->Root, aChild );
//        }
//      }
//      else // I'm not fostered
//        remove_from_fostering_elgibility( aChild ); // I'm a child, soon to be rolled up. I can't be fostering anyone
//      
//      lRealNewParent->ChildSets.push_back( aChild );
//      lRealNewParent->BiologicalPatriarch->NodeCount++; // hey, we're adding one somewhere in your tree.
//      lRealNewParent->TotalKmerCount += aChild->TotalKmerCount; // appropriate this information
//      aChild->ParentSet = lRealNewParent;
//      
//      repoint_to_new_root(lRealNewParent->Root, aChild); // propagate root node info to the kids
//      repoint_to_new_biological_patriarch(lRealNewParent->BiologicalPatriarch, aChild); // if we're rooted in a foster, we need to know. you know, for medical reasons.
//      
//      // keep the tree nice and shallow
//      move_children_up( aChild ); // this guy's children now must become his siblings.
//      move_fosters_up( aChild ); // all this guy's fosters must be parented by the root. The real root (not the foster root)
//
//      rebalance_system_if_needed( lRealNewParent ); // if I'm a foster and overflowing, cut me loose, if I'm otherwise too big collapse or divest me of my foster children.
//      
//      reevaluate_fostering_order_and_availability( lRealNewParent );
//    }
//    void move_fosters_up( SetHandle aNode )
//    {
//      for ( vector<SetHandle>::iterator lChildSet = aNode->ChildSets.begin(); lChildSet != aNode->ChildSets.end(); ++lChildSet )
//      {
//        move_fosters_up( *lChildSet );
//      }
//      for ( vector<SetHandle>::iterator lFosterSet = aNode->FosterSets.begin(); lFosterSet != aNode->FosterSets.end(); ++lFosterSet )
//      {
//        move_fosters_up( *lFosterSet );
//      }
//      
//      if ( aNode != aNode->Root ) // THIS IS DELIBERATE. All Independent Fosters should be made to stand at the same level. No sub-fosters allowed.
//      {
//        for ( vector<SetHandle>::iterator lFosterSet = aNode->FosterSets.begin(); lFosterSet != aNode->FosterSets.end(); ++lFosterSet )
//        {          
//          aNode->Root->FosterSets.push_back( *lFosterSet );
//          (*lFosterSet)->ParentSet = aNode->Root;
//        }
//        aNode->FosterSets.clear();
//      }
//    }
//    void move_children_up( SetHandle aNode )
//    {
//      assert( aNode != aNode->BiologicalPatriarch ); // this should only be done to child nodes, not any root.
//      
//      for ( vector<SetHandle>::iterator lChildSet = aNode->ChildSets.begin(); lChildSet != aNode->ChildSets.end(); ++lChildSet )
//      {
//        move_children_up( *lChildSet );
//        aNode->BiologicalPatriarch->ChildSets.push_back( *lChildSet ); // move this child to top
//        (*lChildSet)->ParentSet = aNode->BiologicalPatriarch;
//      }
//      aNode->ChildSets.clear();
//    }    
//    // introduce all the kids to their new ruler
//    void repoint_to_new_root( SetHandle aNewRoot, SetHandle aSet, int lDepth = 0 )
//    {
//      lDepth++;
//      
//      if ( lDepth > 60 )
//        cout << " (" << lDepth << ") ";
//      
//      aSet->Root = aNewRoot;
//      for ( vector<SetHandle>::iterator lChildSet = aSet->ChildSets.begin(); lChildSet != aSet->ChildSets.end(); ++lChildSet )
//      {
//        repoint_to_new_root( aNewRoot, *lChildSet, lDepth );
//      }
//      for ( vector<SetHandle>::iterator lFosterSet = aSet->FosterSets.begin(); lFosterSet != aSet->FosterSets.end(); ++lFosterSet )
//      {
//        repoint_to_new_root( aNewRoot, *lFosterSet, lDepth );
//      }
//    }    
//    void repoint_to_new_gateway_node( SetHandle aNewRoot, SetHandle aSet )
//    {
//      if ( aSet->StoringHashes )
//        for ( vector<HashIntoType>::iterator lHashIt = aSet->Hashes.begin(); lHashIt != aSet->Hashes.end(); ++lHashIt )
//        {
//          for (int i = 0; i < HASHES; ++i)
//          {
//            unsigned long long lHashBin = HashToHashBinCached(*lHashIt, i );
//            unsigned long long lSetOffsetBin = HashBinToSetOffsetBinCached(lHashBin, i);
//            _set_offsets[i][lSetOffsetBin] = aNewRoot->GetPrimaryGatewaySetOffset();
//          }
//        }
//    }    
//    // propagate the information about their biological grandpa that you just found to all your biological descendents.
//    void repoint_to_new_biological_patriarch( SetHandle aNewBioPatriarch, SetHandle aSet )
//    {
//      aSet->BiologicalPatriarch = aNewBioPatriarch; // I found my real parents!
//      for ( vector<SetHandle>::iterator lChildSet = aSet->ChildSets.begin(); lChildSet != aSet->ChildSets.end(); ++lChildSet )
//      {
//        repoint_to_new_biological_patriarch( aNewBioPatriarch, *lChildSet ); // tell the kids!
//      }      
//    }
//
//    void rebalance_system_if_needed( SetHandle lBioPatriarch )
//    {
//      if ( lBioPatriarch->AmFostered() ) 
//      {
//        if ( is_overflowing( lBioPatriarch ) ) // put this fucker out into its own address
//        {
//          SetHandle lSet = create_allocated_address_set(); 
//          adopt( lSet, lBioPatriarch );
//        }
//      }
//      else 
//      {
//        if ( lBioPatriarch->GetTreeSize() > 100 ) // only non-fosters may ever be collapsed.
//        {
//          collapse( lBioPatriarch );
//        }
//        
//        // only non-fosters can have foster children of their own.
//        if ( lBioPatriarch->GetTotalKmerCount() > 1000 ) // too many kmers to have any foster sets in here, pull-em out.
//        {
//          remove_from_fostering_elgibility( lBioPatriarch );
//          redistribute_foster_children( lBioPatriarch );          
//        }        
//      }
//    }
//    // troll through a sub-tree and see if it's overflowing.
//    bool is_overflowing( SetHandle aNode ) // any individual node or it's children are overflowing.
//    {
//      assert( aNode->StoringHashes );
//      
//      if ( aNode->KmerCount > 100 )
//        return true;
//      
//      for ( vector<SetHandle>::iterator lChildSet = aNode->ChildSets.begin(); lChildSet != aNode->ChildSets.end(); ++lChildSet )
//      {
//        bool lVal = is_overflowing( *lChildSet );
//        
//        if ( lVal == true )
//          return true;
//      }
//      
//      return false;
//    }
//    
//    void redistribute_foster_children( SetHandle aSet )
//    {      
//      assert( aSet->AmRoot() ); // we shouldn't be redistributing foster children if we aren't the head of their household
//      
//      SetHandle lPatsySet = get_least_crowded_set(); // get the smallest set we've got.
//      if ( lPatsySet == NULL )
//        lPatsySet = create_allocated_address_set();      
//      assert ( lPatsySet != NULL ); 
//      
//      // and pawn all of the fosters off on it.
//      for ( int i = aSet->FosterSets.size() - 1; i >= 0; --i ) // leave one foster kid in there. We're going backward because the adoption process rips them out of the foster sets array.
//        foster( lPatsySet, aSet->FosterSets[i] );
//    }
//    
////    void collapse_if_needed( SetHandle aSet )
////    {
////      if ( ready_for_collapse( aSet ) )
////        collapse( aSet );
////    }
//    
////    void redistribute_foster_children_if_needed( SetHandle aSet )
////    {
////      if ( ready_for_foster_redistribution( aSet ) )
////        redistribute_foster_children( aSet );
////    }
//    
//
//    
//    
//    //
//    // Functions related to tree collapse
//    //
////    bool ready_for_collapse( SetHandle aSet )
////    {
////      if ( aSet->GetTreeSize() > 100 )
////        return true;
////      
////      return false;
////    }    
//    void collapse( SetHandle aSet )
//    {
//      assert( aSet->AmBiologicalPatriarch() ); // we shouldn't be collapsing things that aren't the heads of their household
//
//      roll_up_all_child_sets( aSet );
//      aSet->NodeCount = 1; // reset the damned thing.
//
//      stop_counting_if_needed( aSet );
//    }
//    void roll_up_all_child_sets( SetHandle aNode )
//    {
//      for ( vector<SetHandle>::iterator lChildSet = aNode->ChildSets.begin(); lChildSet != aNode->ChildSets.end(); ++lChildSet )
//      {
//        roll_up_all_child_sets( *lChildSet );                
//      }
//      aNode->ChildSets.clear();
//      
//      if ( aNode != aNode->BiologicalPatriarch )
//      {
//        if ( aNode->BiologicalPatriarch->StoringHashes ) // fosters will always have this set to true.
//          for ( vector<HashIntoType>::iterator lHash = aNode->Hashes.begin(); lHash != aNode->Hashes.end(); ++lHash )          
//            aNode->BiologicalPatriarch->Hashes.push_back( *lHash );
//        
//        aNode->Hashes.clear();
//        aNode->BiologicalPatriarch->KmerCount += aNode->KmerCount;
//        
//        // in fostered nodes, this is empty
//        for ( vector<SetOffset>::iterator lSetOffset = aNode->SetOffsetsIOccupy.begin(); lSetOffset != aNode->SetOffsetsIOccupy.end(); ++lSetOffset )
//        {
//          _sets[ *lSetOffset ] = aNode->BiologicalPatriarch; // re-point
//          aNode->BiologicalPatriarch->SetOffsetsIOccupy.push_back( *lSetOffset );
//        }
//        
//        delete aNode;
//      }
//    }
//    
////    bool ready_for_foster_redistribution( SetHandle aSet )
////    {
////      if ( aSet->IsTooCrowded() ) // too crowded.
////        return true;
////      
////      return false;
////    }    
//
//    
//    void stop_counting_if_needed( SetHandle aNode )
//    {
//      if ( !aNode->AmFostered() && aNode->KmerCount > 100 && aNode->StoringHashes )
//      {
//        aNode->StoringHashes = false;
//        aNode->Hashes.clear();
//      }
//      
//    }
//    
////    void collapse_fosters( SetHandle aSet )
////    {
////      for ( vector<SetHandle>::iterator lFosterSet = aSet->FosterSets.begin(); lFosterSet != aSet->FosterSets.end(); ++lFosterSet )
////      {
////        collapse( *lFosterSet );
////      }
////    }
//    
//
//    
////    void cut_loose_too_big_fosters( SetHandle aNode )
////    {
////      for ( int i = aNode->FosterSets.size() - 1; i >= 0; --i )
////      {
////        if ( is_overflowing( aNode->FosterSets[i] ) )
////        {
////          cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~move too big foster to the floor~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
////          
////          
////          SetHandle lSetToJoin = create_allocated_address_set(); // if we're lucky, this'll pull us out onto the main floor with our own address, but it might not. at worst, it'll put us in with a less crowded set.
////                    
////          cout << "set created for this purpose: "; lSetToJoin->output_info();
////          
////          assert ( !lSetToJoin->AmFostered() ); // SHIT
////          
////          adopt( lSetToJoin, aNode->FosterSets[i] );
////          
////          cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~move too big foster to the floor~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
////        }
////      }
////    }
//    
//
//    
//    
//    
//    //
//    // helper functions
//    //
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
//          
//          HashBin lDELETEME = _hash_bin_cache[i][index].second;
//          
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
//    HashBin HashToHashBin( HashIntoType aHash, int i ) // no idea which is faster.
//    {
//      return (aHash % _tablesizes[i]);
//    }
//    
//    // calculate the set offset bin from a hash bin
//    // store the value for future reference.
//    SetOffsetBin HashBinToSetOffsetBinCached( HashBin aBin, int i ) // no idea if this will be faster
//    {
//      for (int j = 0, index = _set_offset_bin_cache_last_used_index[i]; j < CACHESIZE; ++j, ++index)
//      {
//        if ( index >= CACHESIZE )
//          index = 0;
//        
//        if ( _set_offset_bin_cache[i][index].first == aBin )
//        {
//          _set_offset_bin_cache_last_used_index[i] = index;
//          return _set_offset_bin_cache[i][index].second; // found it
//        }
//      }
//      
//      // didn't find it.
//      HashBin lSetOffsetBin = HashBinToSetOffsetBin(aBin, i);
//      
//      // insert it into the cache;
//      _set_offset_bin_cache_last_used_index[i]++;
//      if ( _set_offset_bin_cache_last_used_index[i] >= CACHESIZE )
//        _set_offset_bin_cache_last_used_index[i] = 0;
//      
//      _set_offset_bin_cache[i][ _set_offset_bin_cache_last_used_index[i] ] = pair<HashBin,SetOffsetBin>(aBin, lSetOffsetBin);
//      
//      return lSetOffsetBin;  
//    }    
//    SetOffsetBin HashBinToSetOffsetBin( HashBin aBin, int i )
//    {      
//      assert( aBin < _tablesizes[i] ); // make sure it's a valid bin.
//      
//      assert( _hash_table[i]->Get( aBin ) == true );
//      
//      unsigned long long lBinSectionIndex = (aBin / BIT_COUNT_PARTITION); // the index of the section before the one we're in
//      assert( lBinSectionIndex <= (_tablesizes[i] / BIT_COUNT_PARTITION)); 
//      
//      unsigned long long lSetOffsetBin = _hash_table[i]->CountBits( lBinSectionIndex * BIT_COUNT_PARTITION, aBin );
//      
//      if ( lBinSectionIndex > 0 )      
//        lSetOffsetBin += _hash_table_bit_counts_lookup[i][ lBinSectionIndex - 1 ];
//      
//      assert( lSetOffsetBin > 0 );
//      
//      lSetOffsetBin -= 1; // to make it an array index rather than a count.
//      
//      return lSetOffsetBin;
//    }
//    
//    // find the set that a given hash goes with, based on its offset and hash
//    SetHandle SetOffsetAndHashToSet( SetOffset aOffset, HashIntoType aHash )
//    {      
//      assert( aOffset > 0 ); // we should not be poking here if we don't have a set to go to
//      
//      SetHandle lGatewayNode = _sets[ aOffset ]; // pick up the gateway node
//      
//      SetHandle lRoot = lGatewayNode->Root; // go to the top of the tree
//      
//      SetHandle lTargetSet = find_set_harboring_this_kmer( lRoot, aHash ); // start at the top
//      
//      if ( lTargetSet == NULL ) // means we're in the right place, but we've long since lost track of it. This is fine.
//        return lRoot; // the patriarch takes responsibility for it.
//      
//      return lTargetSet;      
//    }
//    
//    SetHandle find_set_harboring_this_kmer( SetHandle aNode, HashIntoType aHash )
//    {
//      // check all the hashes stored in this node.
//      if ( aNode->AmFostered() ) // I'm rooted in a foster tree, so it matters if I'm holding you.
//      {
//        for ( vector<HashIntoType>::iterator lIt = aNode->Hashes.begin(); lIt != aNode->Hashes.end(); ++lIt )
//        {
//          if ( aHash == *lIt )
//            return aNode;
//        }
//      }
//      
//      // now, ask the child/foster nodes to check.
//      SetHandle lSet = NULL;        
//      for ( vector<SetHandle>::iterator lChildSet = aNode->ChildSets.begin(); lChildSet != aNode->ChildSets.end(); ++lChildSet )
//      {
//        lSet = find_set_harboring_this_kmer( *lChildSet, aHash );
//        if ( lSet != NULL )
//          return lSet;
//      }
//      for ( vector<SetHandle>::iterator lFosterSet = aNode->FosterSets.begin(); lFosterSet != aNode->FosterSets.end(); ++lFosterSet )
//      {
//        lSet = find_set_harboring_this_kmer( *lFosterSet, aHash );
//        if ( lSet != NULL )
//          return lSet;
//      }
//      
//      return NULL; // This node is not specifically holding it, and it's not specifically anywhere in any of this node's downstream trees.
//    }
//    
//    bool is_prime( unsigned long long aCandidate )
//    {
//      if ( aCandidate < 2 )
//        return false;
//      
//      if ( aCandidate == 2 )
//        return true;
//      
//      if ( aCandidate % 2 == 0 )
//        return false;
//      
//      for ( unsigned long long i = 3; i < pow((double)aCandidate, 0.5) + 1; i+= 2 )
//      {
//        if ( aCandidate % i == 0 )
//          return false;
//      }
//      
//      return true;  
//    }
//    
//    unsigned long long get_first_prime_below( unsigned long long aNumber )
//    {
//      unsigned long long i = aNumber - 1;
//      
//      if ( i % 2 == 0 ) // no even
//        --i;
//      
//      while ( i > 0 )
//      {
//        if ( is_prime( i ) )
//          return i;
//        
//        i -= 2;
//      }
//      
//      return 0;
//    }
//    
//    unsigned long long get_first_prime_above( unsigned long long aNumber )
//    {
//      unsigned long long i = aNumber + 1;
//      
//      if ( i % 2 == 0 ) // no even
//        ++i;
//      
//      while ( true )
//      {
//        if ( is_prime( i ) )
//          return i;
//        
//        i += 2;
//      }
//    }
//    
//  };  
//}
