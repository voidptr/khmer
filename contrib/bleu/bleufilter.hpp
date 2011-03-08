#ifndef BLEUFILTER_HPP
#define BLEUFILTER_HPP

#include "../../lib/hashtable.hh"
#include "../../lib/parsers.hh"
#include "CanonicalSetsManager.hpp"
//#include "SequenceHash.hpp"
#include "KmerJoins.hpp"
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
    SequenceHashArbitrary_Builder * _hash_builder;

  public:
    typedef unsigned int (BleuFilter::*ConsumeStringFN)( string * aReads, int aReadCount );
    
    BleuFilter(WordLength ksize, unsigned long long aMaxMemory)
    : Hashtable(ksize)
    { 
      _Sets_Manager = new CanonicalSetsManager( aMaxMemory );
      _hash_builder = new SequenceHashArbitrary_Builder();
      
      // housekeeping
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

    // fucking duplicate code drives me nuts. I swear I will clean this up once 
    // I get some functionality that I'm happy with.
    void consume_strings_for_hash_table(const std::string &filename)
    {
      cout << "Populating Hash Table(s)..." << endl;
      
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
            
          }
        }

        for ( int j = 0; j < lCount; ++j )
        { 
          const unsigned int length = reads[j].length();
          
          for ( unsigned int k = 0; k < length - _ksize + 1; ++k )
          {
            string lSeq = reads[j].substr(k, _ksize);  
            SequenceHashArbitrary hash( lSeq, _hash_builder->hash( lSeq ) );
            _Sets_Manager->seen_hash( hash );
          }
        }
        
        total_reads += lCount;        
        cout << total_reads << endl;
      }
      
      _Sets_Manager->finalize_seen_hash();
      end = time(NULL);        
      cout << "Elapsed: " << difftime(end, start)<< " seconds" << endl;
    }
    
    // fucking duplicate code drives me nuts. I swear I will clean this up once 
    // I get some functionality that I'm happy with.
    void generate_sets(const std::string &filename)
    {
      cout << endl << "Generating Sets..." << endl;
      time_t start, end;      
      start = time(NULL);
      
      int total_reads = 0;
      
      IParser* parser = IParser::get_parser(filename.c_str());
      Read read;
      
      string currName = "";
      string currSeq = "";
      
      string reads[100000];
      
      while(!parser->is_complete())  
      {
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
          const unsigned int length = reads[j].length();
          
          SetHandle lWorkingSet = NULL;
          for ( unsigned int k = 0; k < length - _ksize + 1; ++k )
          {
            string lSeq = reads[j].substr(k, _ksize);
            SequenceHashArbitrary lHash( lSeq, _hash_builder->hash( lSeq ) );
            
            lWorkingSet = _Sets_Manager->assign_hash_to_set( lWorkingSet, lHash );
          }
        }
        total_reads += lCount;        
        
        cout << total_reads << endl;
      }
      
      end = time(NULL);        
      std::cout << "Elapsed: " << difftime(end, start)<< " seconds" << std::endl;
    }

    // fucking duplicate code drives me nuts. I swear I will clean this up once 
    // I get some functionality that I'm happy with.
    void output_join_reads(const std::string &filename, const std::string &outputfilename)
    {
      cout << endl << "Outputting Join Reads..." << endl;
    
      ofstream outfile(outputfilename.c_str());
    
      time_t start, end;      
      start = time(NULL);
      
      int total_reads = 0;
      
      IParser* parser = IParser::get_parser(filename.c_str());
      Read read;
      
      string currName = "";
      string currSeq = "";
      
      string reads[100000];
      string names[100000];
      
      while(!parser->is_complete())  
      {
        int lCount = 0;
        for (int i = 0; i < 100000 && !parser->is_complete(); ++i)
        {
          read = parser->get_next_read();
          if ( check_read( read.seq ) )
          {
            reads[lCount] = read.seq;
            names[lCount] = read.name;
            lCount++;
          }
        }
        
        for ( int j = 0; j < lCount; ++j )
        {        
          const unsigned int length = reads[j].length();
          
          SetHandle lWorkingSet = NULL;
          for ( unsigned int k = 0; k < length - _ksize + 1; ++k )
          {
            string lSeq = reads[j].substr(k, _ksize);
            SequenceHashArbitrary lHash( lSeq, _hash_builder->hash( lSeq ) );          

            // WARNING - this section is the content of assign_hash_to_set, but with the output line. :/
            // can this kmer have a set. Is it interesting?
            if ( _Sets_Manager->can_have_set( lHash ) )
            {
              // does this kmer have an existing set?
              if ( _Sets_Manager->has_existing_set( lHash ) )
              {
                SetHandle lFoundSet = _Sets_Manager->get_existing_set( lHash );
                
                if ( lWorkingSet == NULL ) // woohoo!
                  lWorkingSet = lFoundSet;
                else if ( _Sets_Manager->sets_are_disconnected( lFoundSet, lWorkingSet ) )                
                  lWorkingSet = _Sets_Manager->combine_sets( lFoundSet, lWorkingSet );
                
                outfile << ">" << names[j] << "\t"
                << " " << lSeq << " " << "\n" 
                << reads[j] << "\n";
              }
              else // this hash is brand new, never before seen. :)
              {
                if ( lWorkingSet == NULL )
                  lWorkingSet = _Sets_Manager->get_new_set( lHash );
                else
                  _Sets_Manager->add_to_set( lWorkingSet, lHash );
              }              
            }
            // END WARNING
          }
        }
        total_reads += lCount;        
        
        cout << total_reads << endl;
      }
      
      end = time(NULL);        
      std::cout << "Elapsed: " << difftime(end, start)<< " seconds" << std::endl;
    }
    
    virtual unsigned int output_partitioned_file(const std::string infilename,
                                                 const std::string outputfilename,
                                                 CallbackFn callback=0,
                                                 void * callback_data=0)
    {
      cout << endl << "Outputting Partitioned File..." << endl;
      
      IParser* parser = IParser::get_parser(infilename);
      ofstream outfile(outputfilename.c_str());
      
      unsigned int total_reads = 0;
      unsigned int reads_kept = 0;
      
      Read read;
      string seq;
      
      map<unsigned int, unsigned int> lReadCounts;
      map<unsigned int, unsigned int> lFosteredCounts;
      
      while(!parser->is_complete()) {
        read = parser->get_next_read();
        seq = read.seq;
        
        if (check_read(seq)) 
        {
          const unsigned int length = seq.length();          
          
          SetHandle lSet = NULL;
          for ( unsigned int i = 0; i < length - _ksize + 1; ++i )
          {
            string lSeq = seq.substr(i, _ksize);  
            SequenceHashArbitrary hash( lSeq, _hash_builder->hash( lSeq ) );
            
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
      
      cout << endl;
      cout << "SUMMARY:" << endl;
      cout << "Partitions w/ more than 100 reads, and unhomed reads" << endl;
      cout << setw(20) << "Part #" << setw(10) << "Reads" << setw(20) << "Kmers (experimental)" << endl;
      
      int lSummary[100];
      memset(lSummary, 0, sizeof(int)*100);
      
      for ( map<unsigned int, unsigned int>::iterator lIt = lReadCounts.begin(); lIt != lReadCounts.end(); ++lIt )
      {
        if ( lIt->second > 100 || lIt->first == 0) // does this partition have more than 100 reads in it?
        {
          // set number
          if ( lIt->first == 0 )
            cout << setw(20) << "unhomed (no partition)";
          else
            cout << setw(20) << lIt->first;
            
          // read count
          cout << setw(10) << lIt->second;
          
          // if we have a set number, do the rest of this crap.
          if ( lIt->first > 0 )
          {
            cout << setw(20) << (*(_Sets_Manager->_sets[ lIt->first ]))->KmerCount;
            
            if ( (*(_Sets_Manager->_sets[ lIt->first ]))->JoinOfConvenience )
            {
              cout << setw(10) << "JoC";
              
              if ( (*(_Sets_Manager->_sets[ lIt->first ]))->KmerCount >= _Sets_Manager->_average_size_at_sort )
                cout << setw(10) << "*";
            }
          }
          cout << endl;
        }
        else
        {
          lSummary[ lIt->second ]++;
        }
      }
      
      cout << endl;      
      cout << "Summary of smaller partitions" << endl;
      cout << setw(10) << "Read Ct" << setw(10) << "# partitions" << endl;
      for ( int i = 0; i < 100; ++i )
      {
        if ( lSummary[i] > 0 )
          cout << setw(10) << i << setw(10) << lSummary[i] << endl;
      }
      
      cout << endl;
      cout << endl;
      cout << setw(6) << "total unique set count: "<< lReadCounts.size() << endl;
      cout << endl;
      
      delete parser; parser = NULL;
      
      return lReadCounts.size();
    }
    
    ///////////////////////////
    /// ANALYZE JOINED READS   
    /////////////////////////// 

    struct kmer_info {
      unsigned int number_of_joins;
      map<string, string> reads_containing_kmer;
    };
    virtual unsigned int analyze_joined_reads(const std::string infilename)
    {
      IParser* parser = IParser::get_parser(infilename);
      
      Read read;
      string seq;
      string kmer;
      string readname;
      
      typedef map<string, kmer_info> Kmer_with_reads;
      typedef map<unsigned int, Kmer_with_reads> Set_to_Kmer;
      
      Set_to_Kmer lSets; // twoD map, [SetID][Kmer]->(count,reads)
      
      while(!parser->is_complete()) {
        read = parser->get_next_read();
        seq = read.seq;
        kmer = read.name.substr( read.name.length() - _ksize - 1, _ksize ); // extract the kmer from where we stashed it.
        readname = read.name.substr(1, read.name.find(' ') );

        if (check_read(seq)) 
        {
          // loop through the reads in the file, and figure out what set they go in. Then, count it up.
          SequenceHashArbitrary hash( kmer, _hash_builder->hash( kmer ) );
          SetHandle lSet = _Sets_Manager->get_existing_set( hash );

          unsigned int lSetID = lSet->GetPrimarySetOffset();
          
          if ( lSets.find(lSetID) == lSets.end() )
            lSets[lSetID][kmer].number_of_joins = 0;
          else if (lSets[lSetID].find(kmer) == lSets[lSetID].end() )
            lSets[lSetID][kmer].number_of_joins = 0;
          
          lSets[lSetID][kmer].number_of_joins++; 
          lSets[lSetID][kmer].reads_containing_kmer.insert( pair<string,string>(readname, seq) );                        
        }
      }
      cout << "Kmers that Joined sets and are unmatched beyond the kmer" << endl;

      
      // sets
      for ( Set_to_Kmer::iterator lSetIt = lSets.begin(); lSetIt != lSets.end(); ++lSetIt )
      {        
        // if there are more than one kmer in the join
        // or, if there's only one kmer, that it was joined against more than 
        // once.
       
        if ( lSetIt->second.begin()->second.number_of_joins > 1 )
        {
        
          stringstream lWholeOutput;
          bool lFoundBadKmers = false;
                  
          lWholeOutput << "Set " << lSetIt->first << ": " << lSetIt->second.size() << " kmers." << endl;

          
          // loop through the kmers in the set
          // but first, sort them by how many joins they have.
          multimap<unsigned int, string> lKmersSortedByAbundance;
          for ( Kmer_with_reads::iterator lKmerIt0 = lSetIt->second.begin(); lKmerIt0 != lSetIt->second.end(); ++lKmerIt0 )
          {
            lKmersSortedByAbundance.insert( pair<unsigned int, string>( lKmerIt0->second.number_of_joins, lKmerIt0->first));
          }
          
          // now, run throught the automatically sorted multimap
          for ( multimap<unsigned int, string>::reverse_iterator lKmerAbundanceIt = lKmersSortedByAbundance.rbegin(); 
            lKmerAbundanceIt != lKmersSortedByAbundance.rend(); ++lKmerAbundanceIt )
          {
            Kmer_with_reads::iterator lKmerIt = lSetIt->second.find( lKmerAbundanceIt->second );
          
            if ( lKmerIt->second.number_of_joins > 2 ) // number of joins
            {
            
              // pull out the reads and their names.
              vector<string> lNames;
              vector<string> lReads;
              for ( map<string,string>::iterator lReadIt = lKmerIt->second.reads_containing_kmer.begin(); 
                lReadIt != lKmerIt->second.reads_containing_kmer.end(); ++lReadIt )
              {
                lNames.push_back( lReadIt->first );
                lReads.push_back( lReadIt->second );
              }
            
              KmerJoins lKmerJoin( lKmerIt->first, lReads, lNames, lKmerIt->second.number_of_joins );
                          
              if ( lKmerJoin.S() > 0 ) // there are differences.
              {
                if ( lKmerJoin.Pi() >= 0.05 ) // some arbitrary threshold of difference.
                {
                  lWholeOutput << lKmerJoin.OutputKmerInfo() << endl;                                  
                  lFoundBadKmers = true;
                }
              }
            } 
          }
          
          if ( lFoundBadKmers )
          {
            cout << lWholeOutput.str();
            cout << "============================================================================" << endl;
            cout << "============================================================================" << endl;
          }

        }
      }
      
      cout << endl;
      cout << "DONE." << endl;
      delete parser; parser = NULL;
      
      return lSets.size();
    }
  };
  
  typedef unsigned int (BleuFilter::*ConsumeStringFN)
    ( string * aReads, int aReadCount );
  
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

    int lConsumed = (lArgs->aThis->*lArgs->aMethod)
      ( lArgs->aReads, lArgs->aCount );
    cout << "Chunk " << lArgs->aChunk << " consumed: " << lConsumed << endl;
    
    return (void*) lConsumed;
  }
}



#endif //BLEUFILTER_HPP

