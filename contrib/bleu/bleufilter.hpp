#ifndef BLEUFILTER_HPP
#define BLEUFILTER_HPP

#include "../../lib/hashtable.hh"
#include "../../lib/parsers.hh"
#include "../external_lib/cBitArray.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <sstream>
#include <set>
#include <assert.h>
#include <time.h>

#define CALLBACK_PERIOD 10000
#define KMER_BIT_COUNT_PARTITION 10000


namespace bleu {
  
  extern time_t global_span_start;
  extern time_t last_span_time;
  extern time_t last_span_time_before_switch;
  extern time_t last_collapse_time;
  extern unsigned int collapse_threshold;
  extern short basement;
  
  using namespace khmer;
  using namespace std;
  
  typedef unsigned int SetID;
  
  class BleuFilter
  : public Hashtable
  {
    class Set
    {
    public:
      set<SetID> mSetIDs;
      
      static time_t thingy_start;
      
      Set( SetID aSet )      
      {        
        mSetIDs.insert( aSet );
      }
            
      SetID getCurrentPrimarySetID()
      {
        return *(mSetIDs.begin());      
      }
      
      SetID join( SetID & aSetID, map<SetID, Set*> & aAllSetsByID, SetID * aAllSetIDsByBin, HashIntoType aNumberOfUniqueKmers)
      {
        Set * lSet = aAllSetsByID[ aSetID ];
                
        for (set<SetID>::iterator lIt = lSet->mSetIDs.begin(); lIt != lSet->mSetIDs.end(); ++lIt)
        { 
          aAllSetsByID[*lIt] = this; // re-point all the old pointers.
          mSetIDs.insert( *lIt );
        }
        
        delete lSet; // this should work.
        
        if ( mSetIDs.size() > collapse_threshold )
        {
          last_span_time = difftime( time(NULL), global_span_start );
          
          cout << "LAST SPAN TOOK " << last_span_time << " seconds" << endl;
          
          time_t start, end;
          start = time(NULL);
          
          cout << "COLLAPSING " << mSetIDs.size() << endl;
          SetID lCollapsedID = *(mSetIDs.begin());
          
          int lCollapsed = 0;
          int lEntries = 0;
          
          for (HashIntoType i = 0; i < aNumberOfUniqueKmers; ++i)
            //for (map<HashIntoType, SetID>::iterator lIt = aAllSetIDsByBin.begin(); lIt != aAllSetIDsByBin.end(); ++lIt)
          { 
            if ( aAllSetIDsByBin[ i ] > 0 )
            {
              if ( aAllSetIDsByBin[ i ] != lCollapsedID 
                  && mSetIDs.find( aAllSetIDsByBin[ i ] ) != mSetIDs.end() )
                //if ( lIt->second != lCollapsedID && mSetIDs.find( lIt->second ) != mSetIDs.end() )
              {
                
                aAllSetIDsByBin[ i ] = lCollapsedID;
                //lIt->second = lCollapsedID;
                lCollapsed++;                
              }
              lEntries++;
            
              assert( aAllSetsByID[ aAllSetIDsByBin[ i ] ] > 0 );
              //assert ( aAllSetsByID[ lIt->second ] > 0 );
              
            }
            
          }
          
          set<SetID>::iterator lIt2 = mSetIDs.begin();
          lIt2++;
          
          for (; lIt2 != mSetIDs.end(); ++lIt2)
          {
            aAllSetsByID.erase( *lIt2 ); // remove the extra entry in the map for this set.
          }
          
          mSetIDs.erase( mSetIDs.begin(), mSetIDs.end() ); // erase the whole list
          mSetIDs.insert( lCollapsedID );
          
          end = time(NULL);
          
          last_collapse_time = difftime(end, start);
          
          cout << "COLLAPSING " << lCollapsed << " out of " << lEntries <<  " ENTRIES TOOK " << last_collapse_time << " seconds" << std::endl;          
          
          global_span_start = time(NULL); // reset this puppy.
                    
          if ( last_span_time > 10 && last_span_time > (last_collapse_time * 2) && collapse_threshold > basement )
          {
            last_span_time_before_switch = last_span_time;
            last_span_time = 0;
            last_collapse_time = 0;
            collapse_threshold /= 2;
            cout << "SHRINKING COLLAPSE THRESHOLD TO: " << collapse_threshold << endl;
          }
          else if ( last_span_time < last_collapse_time ) 
          {
            if ( last_span_time > 5 && ((last_span_time * 2) + last_collapse_time ) < last_span_time_before_switch )
              cout << "MIGHT WANT TO RECONSIDER EXPANDING THIS ONE." << endl;
            else
            {
              last_span_time_before_switch = last_span_time;
              last_span_time = 0;
              last_collapse_time = 0;
              collapse_threshold *= 2;
              cout << "EXPANDING COLLAPSE THRESHOLD TO: " << collapse_threshold << endl;
            }
          }
          
          
          return lCollapsedID;
          
        }
        else
          return aSetID;
      }
    };
    
  protected:    
    SetID _last_set;
    
    HashIntoType _table2size;

    cBitArray * _hash_table;    
    cBitArray * _hash_table_2;    

    unsigned int * _hash_bit_counts_lookup; // the number of bits set in the hash table, by every 10k entries

    SetID * _set_IDs;
    HashIntoType _total_unique_hash_count;

    map<SetID, Set*> _sets;
    set<Set *> mUniqueSets;
    
  public:
    BleuFilter(WordLength ksize, HashIntoType tablesize)
    : Hashtable(ksize, (tablesize * .66) + 1 )
    {
      _table2size = _tablesize * .33 - 1;
      _last_set = 0; // zero == none in use. any number > 0 is a set. _last_set indicates the last set that was allocated. It will always be the largest set number.
            
      _hash_table = new cBitArray( _tablesize );
      _hash_table->Clear();      
      
      _hash_table_2 = new cBitArray( _table2size );
      _hash_table_2->Clear();
            
      _hash_bit_counts_lookup = new unsigned int[(_tablesize / KMER_BIT_COUNT_PARTITION)+1];
      memset(_hash_bit_counts_lookup, 0, ((_tablesize / KMER_BIT_COUNT_PARTITION)+1) * sizeof(unsigned int));

      _set_IDs = NULL; // THIS WILL GET SET LATER
      _total_unique_hash_count = 0; // THIS WILL GET SET LATER  
      
      _counts = NULL;   
    }
    
    // consume_string: run through every k-mer in the given string, & hash it.
    // overriding the Hashtable version to support my new thang.
    unsigned int consume_string(const std::string &s,
                                HashIntoType lower_bound = 0,
                                HashIntoType upper_bound = 0)
    {
      const char * sp = s.c_str();
      const unsigned int length = s.length();
      unsigned int n_consumed = 0;
      
      HashIntoType forward_hash = 0, reverse_hash = 0;
      
      // generate the hash for the first kmer in the read (fair amount of work)
      HashIntoType hash = _hash(sp, _ksize, forward_hash, reverse_hash);
      HashIntoType lHashBin = hash % _tablesize;
      HashIntoType lHash2Bin = hash % _table2size;
      
      _hash_table->Set(lHashBin, true);
      _hash_table_2->Set(lHash2Bin, true);
      
      n_consumed++;
      
      // now, do the rest of the kmers in this read (add one nt at a time)
      for (unsigned int i = _ksize; i < length; i++) {
        HashIntoType next_hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );        
        lHashBin = next_hash % _tablesize;
        lHash2Bin = next_hash % _table2size;
        
        _hash_table->Set(lHashBin, true);
        _hash_table_2->Set(lHash2Bin, true);
        
        n_consumed++;
      }
      
      return n_consumed;
    }
    
    void prepare_set_arrays()
    {
      _total_unique_hash_count = _hash_table->CountBits();
      
      _set_IDs = new SetID[_total_unique_hash_count];  
      memset(_set_IDs, 0, _total_unique_hash_count * sizeof(SetID));
      
      unsigned long long lLookupTableSize = (_tablesize / KMER_BIT_COUNT_PARTITION)+1;
      
      for (unsigned long long i = 0; i < lLookupTableSize; ++i)
      {      
        unsigned long long lSectionStartIndex = i * KMER_BIT_COUNT_PARTITION;
        unsigned long long lSectionStopIndex = ((i + 1) * KMER_BIT_COUNT_PARTITION) - 1;
        if ( lSectionStopIndex >= _tablesize )
          lSectionStopIndex = _tablesize - 1;
      
        // temporarily store the single section count.
        //unsigned int lDeleteMe = _hash_table->CountBits(lSectionStartIndex, lSectionStopIndex);
        _hash_bit_counts_lookup[i] = _hash_table->CountBits(lSectionStartIndex, lSectionStopIndex);
        
        if ( i > 0 ) // apply the summation
        {
          //unsigned int lDeleteMe2 = lDeleteMe + _hash_bit_counts_lookup[ i - 1 ];
          _hash_bit_counts_lookup[i] += _hash_bit_counts_lookup[i-1];
        }
      }
      cout << "DONE PREP" << endl;
    }
    
    //
    // generate_sets: consume a FASTA file of reads, again, to generate the sets.
    //
    
    void generate_sets(const std::string &filename,
                       unsigned int &total_reads,
                       unsigned long long &n_consumed,
                       HashIntoType lower_bound = 0,
                       HashIntoType upper_bound = 0,
                       ReadMaskTable ** orig_readmask = NULL,
                       bool update_readmask = true,
                       CallbackFn callback = NULL,
                       void * callback_data = NULL)
    {
      total_reads = 0;
      n_consumed = 0;
      
      IParser* parser = IParser::get_parser(filename.c_str());
      Read read;
      
      string currName = "";
      string currSeq = "";
      
      //
      // readmask stuff: were we given one? do we want to update it?
      // 
      
      ReadMaskTable * readmask = NULL;
      std::list<unsigned int> masklist;
      
      if (orig_readmask && *orig_readmask) {
        readmask = *orig_readmask;
      }
      
      //
      // iterate through the FASTA file & consume the reads.
      //
      
      while(!parser->is_complete())  {
        read = parser->get_next_read();
        currSeq = read.seq;
        currName = read.name; 
        
        // do we want to process it?
        if (!readmask || readmask->get(total_reads)) {
          
          // yep! process.
          unsigned int this_n_consumed;
          bool is_valid;
          
          this_n_consumed = check_and_process_read_for_set(currSeq,
                                                   is_valid,
                                                   lower_bound,
                                                   upper_bound);
          
          // was this an invalid sequence -> mark as bad?
          if (!is_valid && update_readmask) {
            if (readmask) {
              readmask->set(total_reads, false);
            } else {
              masklist.push_back(total_reads);
            }
          } else {		// nope -- count it!
            n_consumed += this_n_consumed;
          }
        }
        
        // reset the sequence info, increment read number
        total_reads++;
        
        if (total_reads % 100 == 0)
          cout << total_reads << endl;
        
        // run callback, if specified
        if (total_reads % CALLBACK_PERIOD == 0 && callback) {
          try {
            callback("consume_fasta", callback_data, total_reads, n_consumed);
          } catch (...) {
            throw;
          }
        }
      }
      
      
      //
      // We've either updated the readmask in place, OR we need to create a
      // new one.
      //
      
      if (orig_readmask && update_readmask && readmask == NULL) {
        // allocate, fill in from masklist
        readmask = new ReadMaskTable(total_reads);
        
        list<unsigned int>::const_iterator it;
        for(it = masklist.begin(); it != masklist.end(); ++it) {
          readmask->set(*it, false);
        }
        *orig_readmask = readmask;
      }
    }
    
    //
    // check_and_process_read: checks for non-ACGT characters before consuming
    //
    
    unsigned int check_and_process_read_for_set(const std::string &read,
                                                   bool &is_valid,
                                                   HashIntoType lower_bound,
                                                   HashIntoType upper_bound)
    {
      is_valid = check_read(read);
      
      if (!is_valid) { return 0; }
      
      return consume_string_for_set(read, lower_bound, upper_bound);
    }
    
    
    // consume_string: run through every k-mer in the given string, & hash it.
    // overriding the Hashtable version to support my new thang.
    unsigned int consume_string_for_set(const std::string &s,
                                HashIntoType lower_bound = 0,
                                HashIntoType upper_bound = 0)
    {
      set<HashIntoType> lListOfSetBins;
      
      const char * sp = s.c_str();
      const unsigned int length = s.length();
      unsigned int n_consumed = 0;
      
      HashIntoType forward_hash = 0, reverse_hash = 0;
      
      // generate the hash for the first kmer in the read (fair amount of work)
      HashIntoType hash = _hash(sp, _ksize, forward_hash, reverse_hash);
      HashIntoType lHashBin = hash % _tablesize;
      HashIntoType lHash2Bin = hash % _table2size;
      
      HashIntoType lSetBin = HashBinToSetBin(lHashBin);
      
      lListOfSetBins.insert( lSetBin );
      
      SetID set_ID = _set_IDs[ lSetBin ]; // init the set_pointer-number to no set.
      
      assert ( set_ID <= _last_set );
      
      bool lNewSet = false;
      
      if ( set_ID == 0 // we're new
        || _hash_table_2->Get( lHash2Bin ) == false // or, we're not a real match.
        )
      { // there wasn't a set already, provisionally create a new one.
        set_ID = init_new_set();
        lNewSet = true;
        _set_IDs[ lSetBin ] = set_ID;
      }
      
      if ( _hash_table_2->Get(lHash2Bin) == false )
        cout << "FP1" << endl;
        
      n_consumed++;
      
      assert ( set_ID <= _last_set );
            
      // now, do the rest of the kmers in this read (add one nt at a time)
      for (unsigned int i = _ksize; i < length; i++) {
        HashIntoType next_hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );        
        lHashBin = next_hash % _tablesize;
        lHash2Bin = next_hash % _table2size;
        lSetBin = HashBinToSetBin(lHashBin);
        
        SetID lBleh = _set_IDs[ lSetBin ];
        
        //assert (lBleh <= _last_set );
        
        if ( _set_IDs[ lSetBin ] > 0 
          && _hash_table_2->Get( lHash2Bin ) == true // we have a match in hash2 as well, so truer match
          && _set_IDs[ lSetBin ] != set_ID ) // there was one here, and it's not us.
        {
          if ( lNewSet ) // we're just provisional anyway.
          {
            delete_provisional_set(set_ID);
            for ( set<HashIntoType>::iterator lIt = lListOfSetBins.begin(); lIt != lListOfSetBins.end(); ++lIt )
            {              
              _set_IDs[ *lIt ] = _set_IDs[ lSetBin ]; // repoint the previous kmers in the provisional set
            }
            set_ID = _set_IDs[ lSetBin ];
            lNewSet = false; // now we're with the band.
          }
          else // a bridge is actually required here.
          {
            set_ID = bridge_sets( _set_IDs[ lSetBin ], set_ID );
          }
        }
        else if ( _hash_table_2->Get( lHash2Bin ) == false )
        {
          cout << "FP" << endl;
        }
        else if ( _set_IDs[ lSetBin ] == 0 ) // we don't have a set assigned. 
        {
          _set_IDs[ lSetBin ] = set_ID;
        }
        
        lListOfSetBins.insert( lSetBin );

        n_consumed++;
      }
      
      return n_consumed;
    }
        
    HashIntoType HashBinToSetBin( HashIntoType aBin )
    {
      unsigned long long lBinSectionIndex = (aBin / KMER_BIT_COUNT_PARTITION); // the index of the section before the one we're in
      assert( lBinSectionIndex <= (_tablesize / KMER_BIT_COUNT_PARTITION));

      HashIntoType lJustThisSectionCount = _hash_table->CountBits( lBinSectionIndex * KMER_BIT_COUNT_PARTITION, aBin );
      HashIntoType lSetBin = lJustThisSectionCount;

      HashIntoType lDeleteMeToo = 0;
            
      if ( lBinSectionIndex > 0 )
      {
        lDeleteMeToo = _hash_bit_counts_lookup[ lBinSectionIndex - 1 ];
        lSetBin += _hash_bit_counts_lookup[ lBinSectionIndex - 1 ];
      }
      
      lSetBin -= 1; // to make it an array index rather than a count.
      
      assert( lSetBin < _total_unique_hash_count );
      
      SetID lDeleteMe = _set_IDs[ lSetBin ];
      
      //for (unsigned long long i = 0; i < _total_unique_hash_count; ++i)
      //{
      //  SetID lDeleteMe3 = _set_IDs[ i ];
      //  assert( _set_IDs[ i ] <= _last_set );
      //}
      
      assert ( _set_IDs[ lSetBin ] <= _last_set );
      
      return lSetBin;
    }
    
    SetID init_new_set()
    {
      SetID lNewSetID = ++_last_set;
      
      Set * lNewSet = new Set( lNewSetID );
      _sets.insert( map<unsigned int, Set*>::value_type( lNewSetID, lNewSet ) );
      
      return lNewSetID;
    }
    
    void delete_provisional_set(SetID aSetID)
    {      
      --_last_set;
      _sets.erase( aSetID );      
    }
    
    SetID bridge_sets( HashIntoType aEncounteredSetID, SetID aCurrentlyUsingSetID )
    {       
      // the bin already has a set ID, which points to a different Set, we should bridge.
      
        Set * lOriginatingSet = _sets[ aCurrentlyUsingSetID ];        
        Set * lEncounteredSet = _sets[ aEncounteredSetID ];
        
        if ( lOriginatingSet != lEncounteredSet ) // bridge them.
        {
          return lEncounteredSet->join( aCurrentlyUsingSetID, _sets, _set_IDs, _total_unique_hash_count );
        }
        else // they weren't different after all.
        {
          return aCurrentlyUsingSetID;    
        }
      
      
    }
    
    HashIntoType _move_hash_foward( HashIntoType & aOldForwardHash, HashIntoType & aOldReverseHash, const char & aNextNucleotide )
    {
      aOldForwardHash = aOldForwardHash << 2; // left-shift the previous hash over
      aOldForwardHash |= twobit_repr(aNextNucleotide); // 'or' in the current nucleotide
      aOldForwardHash &= bitmask; // mask off the 2 bits we shifted over.
      
      // now handle reverse complement
      aOldReverseHash = aOldReverseHash >> 2;
      aOldReverseHash |= (twobit_comp(aNextNucleotide) << (_ksize*2 - 2));
      
      // pick the better bin of the forward or reverse hashes
      return uniqify_rc(aOldForwardHash, aOldReverseHash);
    }
    
    void output_sets()
    {
      for ( map<unsigned int, Set*>::iterator lIt = _sets.begin(); lIt != _sets.end(); ++lIt )
      {
        mUniqueSets.insert( lIt->second );
      }
            
      cout << setw(6) << "unique set count: "<< mUniqueSets.size() << endl;
      cout << endl;
    }
    
    virtual unsigned int output_partitioned_file(const std::string infilename,
                                                 const std::string outputfilename,
                                                 CallbackFn callback=0,
                                                 void * callback_data=0)
    {
      IParser* parser = IParser::get_parser(infilename);
      ofstream outfile(outputfilename.c_str());
      
      unsigned int total_reads = 0;
      unsigned int reads_kept = 0;
      
      Read read;
      string seq;
      
      std::string first_kmer;
      HashIntoType forward_hash = 0, reverse_hash = 0;
      
      map<SetID, unsigned int> lReadCounts;
      
      while(!parser->is_complete()) {
        read = parser->get_next_read();
        seq = read.seq;
        
        if (check_read(seq)) {
          first_kmer = seq.substr(0, _ksize);
          
          // generate the hash for the first kmer in the read (fair amount of work)
          HashIntoType hash = _hash(first_kmer.c_str(), _ksize, forward_hash, reverse_hash);
          HashIntoType bin = hash % _tablesize;
          
          HashIntoType lSetBin = HashBinToSetBin( bin );    
          Set * lDeleteMe = _sets[ _set_IDs[ lSetBin ] ];
          
          SetID lActualFinalSetID = _sets[ _set_IDs[ lSetBin ] ]->getCurrentPrimarySetID();
          
          lReadCounts[ lActualFinalSetID ]++;
          
          outfile << ">" << read.name << "\t" << lActualFinalSetID
          << " " << "\n" 
          << seq << "\n";
          
          // reset the sequence info, increment read number
          total_reads++;
          
          // run callback, if specified
          if (total_reads % CALLBACK_PERIOD == 0 && callback) {
            try {
              callback("do_truncated_partition/output", callback_data,
                       total_reads, reads_kept);
            } catch (...) {
              delete parser; parser = NULL;
              outfile.close();
              throw;
            }
          }
        }
      }
      
      for ( map<SetID, unsigned int>::iterator lIt = lReadCounts.begin(); lIt != lReadCounts.end(); ++lIt )
      {
        cout << setw(10) << lIt->first;
        cout << setw(10) << lIt->second << endl;        
      }
      
      cout << setw(6) << "unique set count: "<< lReadCounts.size() << endl;
      cout << endl;

      
      delete parser; parser = NULL;
      
      return mUniqueSets.size();
    }
    
  };
}

#endif //BLEUFILTER_HPP

