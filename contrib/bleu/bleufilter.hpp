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
#define THREADCT 1

namespace bleu {
  
  using namespace khmer;
  using namespace std;
  
  typedef unsigned long long MaxSize;
  
  void * ThreadStart( void * aArgs );
  
  class BleuFilter
  : public Hashtable
  {
  private:
    CanonicalSetsManager * _Sets_Manager;

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
        
//    // consume_string: run through every k-mer in the given string, & hash it.
//    // overriding the Hashtable version to support my new thang.
//    unsigned int consume_strings_for_hash_table(string * aReads,
//                                                int aReadCount)
//    {      
//      unsigned int n_consumed = 0;
//
//      for ( int i = 0; i < aReadCount; ++i )        
//      {
//        if ( check_read( aReads[i] ) ) // read's ok
//        {
//          const char * sp = aReads[i].c_str();
//          const unsigned int length = aReads[i].length();
//          
//          HashIntoType forward_hash = 0, reverse_hash = 0;
//          
//          HashIntoType hash = 0;
//          bool lHashGenerated = false;      
//          for (unsigned int i = _ksize - 1; i < length; ++i)
//          {
//            if ( lHashGenerated )
//              hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );
//            else
//            {
//              hash = _hash(sp, _ksize, forward_hash, reverse_hash);
//              lHashGenerated = true;
//            }
//            
//            _Sets_Manager->seen_hash( hash );
//            ++n_consumed;
//            
//          }   
//        }
////        if ( i % 100 == 0 )
////          cout << i << endl;
//      }
//      
//      
//      
//      return n_consumed;
//    }
//    
//    
//    // consume_string: run through every k-mer in the given string, & hash it.
//    // overriding the Hashtable version to support my new thang.
//    unsigned int consume_strings_for_set(
//                                string * reads,
//                                int aReadCount)
//
//    {
//      unsigned int n_consumed = 0;
//      
//      for ( int i = 0; i < aReadCount; ++i )
//      {
//        if ( check_read( reads[i] ) )
//        {
//          const char * sp = reads[i].c_str();
//          const unsigned int length = reads[i].length();
//         
//          HashIntoType forward_hash = 0, reverse_hash = 0;
//          
//          SetHandle lWorkingSet = NULL;
//          HashIntoType hash = 0;
//          bool lHashGenerated = false;
//          
//          for (unsigned int i = _ksize - 1; i < length; ++i)
//          {
//            if ( lHashGenerated )
//              hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );
//            else
//            {
//              hash = _hash(sp, _ksize, forward_hash, reverse_hash);
//              lHashGenerated = true;
//            }
//            
//            ++n_consumed;
//             
//            if ( _Sets_Manager->can_have_set( hash ) )
//            {
//              if ( lWorkingSet == NULL )
//              {
//                lWorkingSet = _Sets_Manager->get_set( hash ); // this'll either find the one it goes in, or create it          
//                continue; // move on
//              }
//              if ( _Sets_Manager->has_existing_set( hash ) )
//              {
//                SetHandle lExistingSet = _Sets_Manager->get_existing_set( hash );
//                if ( _Sets_Manager->sets_are_disconnected( lExistingSet, lWorkingSet ) )
//                  lWorkingSet = _Sets_Manager->bridge_sets( lExistingSet, lWorkingSet );            
//              }
//              else
//                _Sets_Manager->add_to_set( lWorkingSet, hash );
//            }
//          }
//        }
//      }
//      return n_consumed;
//    }
    
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
            lSetID = lSet->GetPrimarySetOffset();
            
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
        if ( lIt->first > 0 )
        {
          cout << setw(10) << (*(_Sets_Manager->_sets[ lIt->first ]))->KmerCount;
          
          if ( (*(_Sets_Manager->_sets[ lIt->first ]))->JoinOfConvenience )
          {
            cout << setw(10) << "JoC";
            
            if ( (*(_Sets_Manager->_sets[ lIt->first ]))->KmerCount >= _Sets_Manager->_average_size_at_sort )
              cout << setw(10) << "*";
          }
        }
        
        cout << endl;
      }
      
      cout << setw(6) << "unique set count: "<< lReadCounts.size() << endl;
      cout << endl;
      
      //cout << "Hash Entries: " << _total_unique_hash_count << endl;
      
      delete parser; parser = NULL;
      
      return lReadCounts.size();
    }
    

    // fucking duplicate code drives me nuts. I swear I will clean this up once I get some functionality that I'm happy with.
    void consume_strings_for_hash_table(const std::string &filename)
    {
      time_t start, end;      
      start = time(NULL);
      
      int total_reads = 0;
      
      IParser* parser = IParser::get_parser(filename.c_str());
      Read read;
      
      string currName = "";
      string currSeq = "";
      
      string reads[100000];
            
      while(!parser->is_complete())  {
        
        int lCount = 0;
        //int lKmerCt = 0;
        for (int i = 0; i < 100000 && !parser->is_complete(); ++i)
        {
          read = parser->get_next_read();
          if ( check_read( read.seq ) )
          {
            reads[lCount] = read.seq;
            lCount++;
            
            //lKmerCt += read.seq.length() - _ksize + 1;
          }
        }

        for ( int j = 0; j < lCount; ++j )
        {        
          const char * sp = reads[j].c_str();
          const unsigned int length = reads[j].length();
          
          HashIntoType forward_hash = 0, reverse_hash = 0;
                    
          HashIntoType hash = 0;
          bool lHashGenerated = false;
          
          for (unsigned int k = _ksize - 1; k < length; ++k)
          {
            if ( lHashGenerated )
              hash = _move_hash_foward( forward_hash, reverse_hash, sp[k] );
            else
            {
              hash = _hash(sp, _ksize, forward_hash, reverse_hash);
              lHashGenerated = true;
            }
            
            _Sets_Manager->seen_hash( hash );
          }
        }
        
        total_reads += lCount;        
        
        cout << total_reads << endl;
      }
      
      _Sets_Manager->finalize_seen_hash();
      
      end = time(NULL);        
      std::cout << "READTEST DONE: " << difftime(end, start)<< " seconds" << std::endl;
    }
    
    
    // fucking duplicate code drives me nuts. I swear I will clean this up once I get some functionality that I'm happy with.
    void generate_sets(const std::string &filename)
    {
      time_t start, end;      
      start = time(NULL);
      
      int total_reads = 0;
      
      IParser* parser = IParser::get_parser(filename.c_str());
      Read read;
      
      string currName = "";
      string currSeq = "";
      
      string reads[100000];
      
      while(!parser->is_complete())  {
        
        int lCount = 0;
        for (int i = 0; i < 100000 && !parser->is_complete(); ++i)
        {
          read = parser->get_next_read();
          if ( check_read( read.seq ) )
          {
            reads[lCount] = read.seq;
            lCount++;
          }
        }
        
        for ( int j = 0; j < lCount; ++j )
        {        
          const char * sp = reads[j].c_str();
          const unsigned int length = reads[j].length();
          
          HashIntoType forward_hash = 0, reverse_hash = 0;
          
          SetHandle lWorkingSet = NULL;
          HashIntoType hash = 0;
          bool lHashGenerated = false;
          
          for (unsigned int k = _ksize - 1; k < length; ++k)
          {
            if ( lHashGenerated )
              hash = _move_hash_foward( forward_hash, reverse_hash, sp[k] );
            else
            {
              hash = _hash(sp, _ksize, forward_hash, reverse_hash);
              lHashGenerated = true;
            }
            
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
              {
                
                _Sets_Manager->add_to_set( lWorkingSet, hash );
              }
            }
          }
        }
        total_reads += lCount;        
        
        cout << total_reads << endl;
      }
      
      end = time(NULL);        
      std::cout << "READTEST DONE: " << difftime(end, start)<< " seconds" << std::endl;
    }
    
  };
  
  typedef unsigned int (BleuFilter::*ConsumeStringFN)( string * aReads, int aReadCount );
  
  typedef struct ThreadArgs {
    ConsumeStringFN aMethod;
    string * aReads;
    int aCount;
    BleuFilter * aThis;
    int aChunk;
  };
  
  void * ThreadStart( void * aArgs )
  {
    ThreadArgs * lArgs = (ThreadArgs*) aArgs;

    int lConsumed = (lArgs->aThis->*lArgs->aMethod)( lArgs->aReads, lArgs->aCount );
    cout << "Chunk " << lArgs->aChunk << " consumed: " << lConsumed << endl;
    
    return (void*) lConsumed;
  }
}



#endif //BLEUFILTER_HPP

