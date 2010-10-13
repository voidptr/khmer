#ifndef BLEUFILTER_HPP
#define BLEUFILTER_HPP

#include "../../lib/hashtable.hh"
#include "../../lib/parsers.hh"
#include "CanonicalSetsManager.hpp"
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
#include <algorithm>
#include <pthread.h>

#define CALLBACK_PERIOD 10000


namespace bleu {
  
  using namespace khmer;
  using namespace std;
  
  typedef unsigned long long MaxSize;
  
  class BleuFilter
  : public Hashtable
  {
  private:
    CanonicalSetsManager * _Sets_Manager;
    
    
    template<class T, void(T::*mem_fn)()>
    void* thunk(void* p)
    {
      (static_cast<T*>(p)->*mem_fn)();
      return 0;
    }

  public:
    typedef unsigned int (BleuFilter::*ConsumeStringFN)( string * aReads, int aReadCount );
    
    BleuFilter(WordLength ksize, unsigned long long aMaxMemory)
    : Hashtable(ksize, 0)
    { 
      _Sets_Manager = new CanonicalSetsManager( aMaxMemory );
      
      // housekeeping
      _counts = NULL;   
    }
    
    void populate_hash_table_bit_count_lookups() {
      _Sets_Manager->populate_hash_table_bit_count_lookups();
    }
    
    void deallocate_hash_table_preliminary() {
      _Sets_Manager->deallocate_hash_table_preliminary();
    }
    
    void allocate_set_offset_table() {
      _Sets_Manager->allocate_set_offset_table();      
    }
        
    // consume_string: run through every k-mer in the given string, & hash it.
    // overriding the Hashtable version to support my new thang.
    unsigned int consume_strings_for_hash_table(string * aReads,
                                                int aReadCount)
    {      
      unsigned int n_consumed = 0;

      for ( int i = 0; i < aReadCount; ++i )        
      {
        if ( check_read( aReads[i] ) ) // read's ok
        {
          const char * sp = aReads[i].c_str();
          const unsigned int length = aReads[i].length();
          
          HashIntoType forward_hash = 0, reverse_hash = 0;
          
          HashIntoType hash = 0;
          bool lHashGenerated = false;      
          for (unsigned int i = _ksize - 1; i < length; ++i)
          {
            if ( lHashGenerated )
              hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );
            else
            {
              hash = _hash(sp, _ksize, forward_hash, reverse_hash);
              lHashGenerated = true;
            }
            
            _Sets_Manager->seen_hash( hash );
            ++n_consumed;
            
          }   
        }
      }
      
      return n_consumed;
    }
    
    
    // consume_string: run through every k-mer in the given string, & hash it.
    // overriding the Hashtable version to support my new thang.
    unsigned int consume_strings_for_set(
                                string * reads,
                                int aReadCount)

    {
      unsigned int n_consumed = 0;
      
      for ( int i = 0; i < aReadCount; ++i )
      {
        if ( check_read( reads[i] ) )
        {
          const char * sp = reads[i].c_str();
          const unsigned int length = reads[i].length();
         
          HashIntoType forward_hash = 0, reverse_hash = 0;
          
          SetHandle lWorkingSet = NULL;
          HashIntoType hash = 0;
          bool lHashGenerated = false;
          
          for (unsigned int i = _ksize - 1; i < length; ++i)
          {
            if ( lHashGenerated )
              hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );
            else
            {
              hash = _hash(sp, _ksize, forward_hash, reverse_hash);
              lHashGenerated = true;
            }
            
            ++n_consumed;
             
            if ( _Sets_Manager->can_have_set( hash ) )
            {
              if ( lWorkingSet == NULL )
              {
                lWorkingSet = _Sets_Manager->get_set( hash ); // this'll either find the one it goes in, or create it          
                continue; // move on
              }
              if ( _Sets_Manager->has_existing_set( hash ) )
              {
                SetHandle lExistingSet = _Sets_Manager->get_existing_set( hash );
                if ( _Sets_Manager->sets_are_disconnected( lExistingSet, lWorkingSet ) )
                  lWorkingSet = _Sets_Manager->bridge_sets( lExistingSet, lWorkingSet );            
              }
              else
                _Sets_Manager->add_to_set( lWorkingSet, hash );
            }
          }
        }
      }
      return n_consumed;
    }
    
    HashIntoType _move_hash_foward( HashIntoType & aOldForwardHash, HashIntoType & aOldReverseHash, const char & aNextNucleotide )
    {
      unsigned long long lTwoBit = twobit_repr( aNextNucleotide );
      
      aOldForwardHash = aOldForwardHash << 2; // left-shift the previous hash over
      aOldForwardHash |= lTwoBit; // 'or' in the current nucleotide
      aOldForwardHash &= bitmask; // mask off the 2 bits we shifted over.
      
      // now handle reverse complement
      aOldReverseHash = aOldReverseHash >> 2;
      aOldReverseHash |= (compl_twobit(lTwoBit) << (_ksize*2 - 2));
      
      // pick the better bin of the forward or reverse hashes
      return uniqify_rc(aOldForwardHash, aOldReverseHash);
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
      
      //std::string first_kmer;
      //HashIntoType forward_hash = 0, reverse_hash = 0;
      
      map<unsigned int, unsigned int> lReadCounts;
      map<unsigned int, unsigned int> lFosteredCounts;
      
      while(!parser->is_complete()) {
        read = parser->get_next_read();
        seq = read.seq;
        
        if (check_read(seq)) {
          
          const char * sp = seq.c_str();
          const unsigned int length = seq.length();
          
          
          SetHandle lSet = NULL;
          HashIntoType hash = 0;
          HashIntoType forward_hash = 0, reverse_hash = 0;          
          bool lHashGenerated = false;          
          for (unsigned int i = _ksize - 1; i < length; ++i) // run through the kmers until you find a set
          {
            if ( lHashGenerated )
              hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );
            else
            {
              hash = _hash(sp, _ksize, forward_hash, reverse_hash);
              lHashGenerated = true;
            }
            
            if ( _Sets_Manager->can_have_set( hash ) )
            {
              lSet = _Sets_Manager->get_existing_set( hash );
              break;
            }
                          
          }   
          
          unsigned int lSetID = 0;
          if ( lSet != NULL )
          {
            lSetID = lSet->PrimarySetOffset;
            
            if ( lSet->JoinOfConvenience )
            {
              lFosteredCounts[ lSetID ]++;
            }
          }
          
          lReadCounts[ lSetID ]++;
          
          outfile << ">" << read.name << "\t" 
          << lSetID 
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
      
      for ( map<unsigned int, unsigned int>::iterator lIt = lReadCounts.begin(); lIt != lReadCounts.end(); ++lIt )
      {
        cout << setw(10) << lIt->first;
        cout << setw(10) << lIt->second;
        
        if ( lFosteredCounts.find( lIt->first ) != lFosteredCounts.end() )
          cout << setw(10) << "T";
        
                
        cout << endl;
      }
      
      cout << setw(6) << "unique set count: "<< lReadCounts.size() << endl;
      cout << endl;
      
      //cout << "Hash Entries: " << _total_unique_hash_count << endl;
      
      delete parser; parser = NULL;
      
      return lReadCounts.size();
    }
    


    // fucking duplicate code drives me nuts. I swear I will clean this up once I get some functionality that I'm happy with.
    void consume_reads(const std::string &filename,
                       unsigned int &total_reads,
                       unsigned long long &n_consumed,
                       ConsumeStringFN consume_string_fn,
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
      
      string reads[10000];
      
      while(!parser->is_complete())  {
        
        int lCount = 0;
        for (int i = 0; i < 10000 && !parser->is_complete(); ++i)
        {
          read = parser->get_next_read();
          reads[i] = read.seq;
          lCount = i + 1;
        }
          
        // yep! process.
        
        unsigned int this_n_consumed = (this->*consume_string_fn)(reads, lCount);
        
//        unsigned int this_n_consumed = check_and_process_read(reads, lCount, consume_string_fn);        

        n_consumed += this_n_consumed;
      
        // reset the sequence info, increment read number
        total_reads += lCount;
        
        cout << total_reads << " " << n_consumed << endl;
        
      }
    }
    
    //
    // check_and_process_read: checks for non-ACGT characters before consuming
    //
    
//    unsigned int check_and_process_read(string * reads,
//                                        int aCount,
//                                        ConsumeStringFN consume_string_fn)
//    {
//      
//      for ( int i = 0; i < aCount; ++i )        
//      {
//        if ( !check_read( reads[i] ) )
//          return 0;
//      }
//      
//      return (this->*consume_string_fn)(reads, aCount);
//    }
  };
}



#endif //BLEUFILTER_HPP

