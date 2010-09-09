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
                
                SetID lDeleteMe2 = aAllSetIDsByBin[ i ];
                set<SetID>::iterator lDeleteMeIt = mSetIDs.find( lDeleteMe2 );
                SetID lDeleteMe3 = *lDeleteMeIt;
                
                aAllSetIDsByBin[ i ] = lCollapsedID;
                //lIt->second = lCollapsedID;
                lCollapsed++;                
              }
              lEntries++;
              
              SetID lDeleteMe = aAllSetIDsByBin[ i ];
              Set * lDeleteMeSet = aAllSetsByID[ lDeleteMe ];
            
              assert( aAllSetsByID[ aAllSetIDsByBin[ i ] ] > 0 );
              //assert ( aAllSetsByID[ lIt->second ] > 0 );
              
            }
            
          }
          
          set<SetID>::iterator lIt2 = mSetIDs.begin();
          lIt2++;
          
          for (; lIt2 != mSetIDs.end(); ++lIt2)
          {
            SetID lToDel = *lIt2;
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
    // fixed-length array
    SetID * _set_IDs;
    // end
    
    cBitArray _kmer_flags;
    unsigned int * _kmer_bit_counts;
    HashIntoType _unique_kmer_count;
    
    map<SetID, Set*> _sets;
    set<Set *> mUniqueSets;
    
  public:
    BleuFilter(WordLength ksize, HashIntoType tablesize)
    
    : Hashtable(ksize, tablesize)
    {
      _last_set = 0; // zero == none in use. any number > 1 is a set. _last_set indicates the last set that was allocated.
      
      delete _counts; // not using these right now. Will later for presence flags.
      _counts = NULL;
      
//      _set_IDs = new SetID[_tablesize];
//      memset(_set_IDs, 0, _tablesize * sizeof(SetID));
      _set_IDs = NULL; // THIS WILL GET SET LATER
      _unique_kmer_count = 0; // THIS WILL GET SET LATER
      
      _kmer_flags.ResizeClear(tablesize);
      _kmer_bit_counts = new unsigned int[(tablesize / KMER_BIT_COUNT_PARTITION)+1];
      
      
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
      HashIntoType bin = hash % _tablesize;
      
      _kmer_flags.Set(bin, true);
      
      n_consumed++;
      
      // now, do the rest of the kmers in this read (add one nt at a time)
      for (unsigned int i = _ksize; i < length; i++) {
        HashIntoType next_hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );        
        bin = next_hash % _tablesize;
        
        _kmer_flags.Set(bin, true);
        
        n_consumed++;
      }
      
      return n_consumed;
    }
    
    void prepare_set_arrays()
    {
      _unique_kmer_count = _kmer_flags.CountBits();
      
      _set_IDs = new SetID[_unique_kmer_count];  
      memset(_set_IDs, 0, _unique_kmer_count * sizeof(SetID));
      
      for (int i = 0; i < (_tablesize / KMER_BIT_COUNT_PARTITION)+1; ++i)
      {      
        // temporarily store the single section count.
        _kmer_bit_counts[i] = _kmer_flags.CountBits((i * KMER_BIT_COUNT_PARTITION), ((i+1) * KMER_BIT_COUNT_PARTITION) - 1);
        
        if ( i > 0 ) // apply the summation
          _kmer_bit_counts[i] += _kmer_bit_counts[i-1];
      }
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
        
        if (total_reads % 1000 == 0)
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
      const char * sp = s.c_str();
      const unsigned int length = s.length();
      unsigned int n_consumed = 0;
      
      HashIntoType forward_hash = 0, reverse_hash = 0;
      
      // generate the hash for the first kmer in the read (fair amount of work)
      HashIntoType hash = _hash(sp, _ksize, forward_hash, reverse_hash);
      HashIntoType bin = hash % _tablesize;
      
      SetID set_ID = 0; // init the set_pointer-number to no set.
      set_ID = initial_set_fetch_or_assignment( bin );
      n_consumed++;
      
      // now, do the rest of the kmers in this read (add one nt at a time)
      for (unsigned int i = _ksize; i < length; i++) {
        HashIntoType next_hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );        
        bin = next_hash % _tablesize;
        set_ID = assign_or_bridge_sets( bin, set_ID );
        
        n_consumed++;
      }
      
      return n_consumed;
    }
    
    SetID initial_set_fetch_or_assignment( HashIntoType aBin )
    {
      HashIntoType lSetBin = HashBinToSetBin( aBin );
      
      if ( _set_IDs[ lSetBin ] == 0 )
        _set_IDs[lSetBin] = init_new_set();
      
      return _set_IDs[lSetBin];
    }
    
    HashIntoType HashBinToSetBin( HashIntoType aBin )
    {
      unsigned int lBinSectionIndex = (aBin / KMER_BIT_COUNT_PARTITION); // the index of the section before the one we're in
      unsigned int lSectionCountsSize = _tablesize / KMER_BIT_COUNT_PARTITION;
      assert( lBinSectionIndex <= (_tablesize / KMER_BIT_COUNT_PARTITION));
      
      unsigned int lStartBit = lBinSectionIndex * KMER_BIT_COUNT_PARTITION;
            

      HashIntoType lJustThisSectionCount = _kmer_flags.CountBits( lBinSectionIndex * KMER_BIT_COUNT_PARTITION, aBin );
      HashIntoType lSetBin = lJustThisSectionCount;
      
      if ( lBinSectionIndex > 0 )
        lSetBin += _kmer_bit_counts[ lBinSectionIndex - 1 ];
      
      lSetBin -= 1; // to make it an array index rather than a count.
      
      assert( lSetBin < _unique_kmer_count );
      return lSetBin;
    }
    
    SetID init_new_set()
    {
      SetID lNewSetID = ++_last_set;
      
      Set * lNewSet = new Set( lNewSetID );
      _sets.insert( map<unsigned int, Set*>::value_type( lNewSetID, lNewSet ) );
      
      return lNewSetID;
    }
    
    SetID assign_or_bridge_sets( HashIntoType aBin, SetID aSetID )
    {
      
      
      // so, we have a bin, and a setID it should go with. Right, sounds good.
      
      // if the bin is empty, that's easy. put in our set ID, and move on.
      
      // otherwise, if the bin already has a set ID, which points to a different Set them, we should bridge.
      
      HashIntoType lSetBin = HashBinToSetBin( aBin );      
      
      if ( _set_IDs[ lSetBin ] == 0 )
        //if ( _set_IDs.find( aBin ) == _set_IDs.end() ) // no set found for this bin
      {
        _set_IDs[lSetBin] = aSetID; // assign
        return aSetID; // yay!
      }
      else 
      {
        Set * lOriginatingSet = _sets[ aSetID ];
        
        SetID lEncounteredSetID = _set_IDs[ lSetBin ];        
        Set * lEncounteredSet = _sets[ lEncounteredSetID ];
        
        if ( lOriginatingSet != lEncounteredSet ) // bridge them.
        {
          SetID lDominatingSetID = lEncounteredSet->join( aSetID, _sets, _set_IDs, _unique_kmer_count );
          return lDominatingSetID;
        }
        else // they weren't different after all.
        {
          return aSetID;    
        }
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

