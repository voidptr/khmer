#ifndef BLEUFILTER_HPP
#define BLEUFILTER_HPP

#include "../../lib/hashtable.hh"
#include "../../lib/parsers.hh"
#include "CanonicalSetsManager.hpp"
//#include "SequenceHash.hpp"
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
    : Hashtable(ksize, 0)
    { 
      _Sets_Manager = new CanonicalSetsManager( aMaxMemory );
      _hash_builder = new SequenceHashArbitrary_Builder();
      
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

    // fucking duplicate code drives me nuts. I swear I will clean this up once 
    // I get some functionality that I'm happy with.
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
            
          }
        }

        for ( int j = 0; j < lCount; ++j )
        {        
          const unsigned int length = reads[j].length();
          
          for ( unsigned int k = 0; k < length - _ksize + 1; ++k )
          {
            string lSeq = reads[j].substr(k, _ksize);  
            SequenceHashArbitrary hash( lSeq, _hash_builder->hash( lSeq ) );
//            cout << hash.canonical_hash << endl;
            
            _Sets_Manager->seen_hash( hash );
          }
        }
        
        total_reads += lCount;        
        
        cout << total_reads << endl;
      }
      
      _Sets_Manager->finalize_seen_hash();
      
      end = time(NULL);        
      cout << "READTEST DONE: " << difftime(end, start)<< " seconds" << endl;
    }
    
    // fucking duplicate code drives me nuts. I swear I will clean this up once 
    // I get some functionality that I'm happy with.
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
          
            if ( _Sets_Manager->can_have_set( lHash ) )
            {
              if ( lWorkingSet == NULL )
              {
                lWorkingSet = _Sets_Manager->get_set( lHash ); // this'll either find the one it goes in, or create it          
                continue; // move on
              }
              if ( _Sets_Manager->has_existing_set( lHash ) )
              {
                SetHandle lExistingSet = _Sets_Manager->get_existing_set( lHash );
                if ( _Sets_Manager->sets_are_disconnected( lExistingSet, lWorkingSet ) )
                {
                  lWorkingSet = _Sets_Manager->bridge_sets( lExistingSet, lWorkingSet );            
                }
              }
              else
              {
                
                _Sets_Manager->add_to_set( lWorkingSet, lHash );
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

    // fucking duplicate code drives me nuts. I swear I will clean this up once 
    // I get some functionality that I'm happy with.
    void output_join_reads(const std::string &filename, const std::string &outputfilename)
    {
    
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
          
            if ( _Sets_Manager->can_have_set( lHash ) )
            {
              if ( lWorkingSet == NULL )
              {
                // if this hash goes with something that's already there, throw the read in there, since it's a kind of join.
                if ( _Sets_Manager->has_existing_set( lHash ) )
                {
                  outfile << ">" << names[j] << "\t"
                  << " " << lSeq << " " << "\n" 
                  << reads[j] << "\n"; 
                }
                  
                lWorkingSet = _Sets_Manager->get_set( lHash ); // this'll either find the one it goes in, or create it          
                continue; // move on
              }
              if ( _Sets_Manager->has_existing_set( lHash ) )
              {
                SetHandle lExistingSet = _Sets_Manager->get_existing_set( lHash );
                if ( _Sets_Manager->sets_are_disconnected( lExistingSet, lWorkingSet ) )
                {
                  // if you are joining two sets, throw the read in there
                  outfile << ">" << names[j] << "\t"
                  << " " << lSeq << " " << "\n" 
                  << reads[j] << "\n"; 
//                  << lSeq << "\n";
                  lWorkingSet = _Sets_Manager->bridge_sets( lExistingSet, lWorkingSet );            
                }
              }
              else
              {
                
                _Sets_Manager->add_to_set( lWorkingSet, lHash );
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
      
      for ( map<unsigned int, unsigned int>::iterator lIt = lReadCounts.begin(); lIt != lReadCounts.end(); ++lIt )
      {
        if ( lIt->first == 0 )
          cout << setw(10) << "none";
        else
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
      
      delete parser; parser = NULL;
      
      return lReadCounts.size();
    }
    
    virtual unsigned int analyze_joined_reads(const std::string infilename)
    {
      IParser* parser = IParser::get_parser(infilename);
      
      Read read;
      string seq;
      string kmer;
      string readname;
      
      //map<unsigned int, unsigned int> lReadCounts;
      
      //  set ID      ,     kmer  , count
      map<unsigned int, map<string, pair<unsigned int, set<string> > > > lJoinedReads;
      
      while(!parser->is_complete()) {
        read = parser->get_next_read();
        seq = read.seq;
        kmer = read.name.substr( read.name.length() - _ksize - 1, _ksize );
        readname = read.name.substr(1, read.name.find(' ') );
        readname.append(" ");
        readname.append(seq);
        
        if (check_read(seq)) 
        {
          const unsigned int length = seq.length();          
         
          // loop through the reads in the file, and figure out what set they go in. Then, count it up.
          SequenceHashArbitrary hash( kmer, _hash_builder->hash( kmer ) );
          SetHandle lSet = _Sets_Manager->get_existing_set( hash );

          unsigned int lSetID = lSet->GetPrimarySetOffset();
          
          if ( lJoinedReads.find(lSetID) == lJoinedReads.end() )
            lJoinedReads[lSetID][kmer].first = 0;
          else if (lJoinedReads[lSetID].find(kmer) == lJoinedReads[lSetID].end() )
            lJoinedReads[lSetID][kmer].first = 0;
          
          lJoinedReads[lSetID][kmer].first++; 
          lJoinedReads[lSetID][kmer].second.insert(readname);
                        
        }
      }
      
      vector<int> lPercentages;
      
      cout << "Kmers that Joined sets" << endl;
      //          setID             kmer         count           reads
      for ( map<unsigned int, map<string, pair<unsigned int, set<string> > > >::iterator lIt = lJoinedReads.begin(); lIt != lJoinedReads.end(); ++lIt )
      {        
        // if there are more than one kmer in the join
        // or, if there's only one kmer, that it was joined against more than once.
       
        if ( lIt->second.size() > 1 || lIt->second.begin()->second.first > 1 )
        {      
        
          cout << "Set " << lIt->first << ": " << lIt->second.size() << " kmers." << endl;
          
          //         kmer          count             reads
          for ( map<string, pair<unsigned int, set<string> > >::iterator lIt2 = lIt->second.begin(); lIt2 != lIt->second.end(); ++lIt2 )
          {
            if ( lIt2->second.first > 1 ) // number of joins
            {
              cout << endl;
              cout << " " << lIt2->second.first << " join points." << endl; // kmer & count
              
              // output the reads they came from, aligned.
              
              // first, divine the alignment
              int lLongestAlignment = 0;
              //        reads                          reads
              for ( set<string>::iterator lIt3 = lIt2->second.second.begin(); lIt3 != lIt2->second.second.end(); ++lIt3 )
              {
                int lAlignmentStart = lIt3->find( lIt2->first ); //lRead.find( lKmer );
                
                if ( lAlignmentStart > lLongestAlignment )
                  lLongestAlignment = lAlignmentStart;
              }
              // then output the alignment
              string lMarker = lIt2->first;
              lMarker.insert(0, lLongestAlignment + 3, ' ');
              cout << lMarker << endl;
              
              vector<string> lLines;
              int lBareAlignment = 0;
              for ( set<string>::iterator lIt3 = lIt2->second.second.begin(); lIt3 != lIt2->second.second.end(); ++lIt3 )
              {
                int leftpad = ( lLongestAlignment - lIt3->find( lIt2->first ) );
                string lRead = *lIt3;
                lRead.insert( lRead.find(" "), leftpad, ' ');
                lRead.insert(0, 3, ' ');
                cout << lRead << endl;

                string lBareAlignedRead = lIt3->substr(lIt3->rfind(" ") + 1);
                
                if ( lBareAlignment < lBareAlignedRead.find( lIt2->first ) )
                  lBareAlignment = lBareAlignedRead.find( lIt2->first );
                lLines.push_back( lBareAlignedRead );
              }
              
              ///// finally, compare the aligned bits
              
              // resize the lines appropriately.
              int lBareAlignedReadMaxLength = 0;
//              cout << endl;
              for (int i = 0; i < lLines.size(); ++i)
              {
                lLines[i].insert(0, lBareAlignment - lLines[i].find(lIt2->first), ' '); // pad left to align

                if ( lLines[i].length() > lBareAlignedReadMaxLength )
                  lBareAlignedReadMaxLength = lLines[i].length();
                  
//                cout << lLines[i] << endl;
              }
              
              // go over each character and score it appropriately
              int lBareAlignedLocation = lLines[0].find( lIt2->first ); // just grab the first one. the others should already be lined up.
              int lMismatchScore = 0;
              int lTotalPositionsCompared = 0;
              for (int j = 0; j < lBareAlignedReadMaxLength; ++j) // each nucleotide
              {
//                bool lCompared = false;
                set<char> lNucleotides;
                int lComparableLines = 0;
                for ( int k = 0; k < lLines.size(); ++k ) // line by line
                {
                  if ( ( j < lBareAlignedLocation || j >= lBareAlignedLocation + lIt2->first.length() ) && lLines[k].length() > j )
                  {

                    if( lLines[k][j] != ' ' )
                    {
                      lNucleotides.insert( lLines[k][j] );
                      //lCompared = true;
                      lComparableLines++;
                    }
                  }
                }    
                if ( lComparableLines > 1 )
                {
                  lTotalPositionsCompared++;
 //                 cout << "x";
                }
//                else
//                {
//                  cout << ".";
//                }
                if ( lNucleotides.size() > 1 )
                  lMismatchScore++;
              }
              if ( lTotalPositionsCompared > 0 )
              {
                int lPercentCorrect =  (lTotalPositionsCompared - lMismatchScore) / (double)lTotalPositionsCompared * 100;
                cout << " " << lTotalPositionsCompared - lMismatchScore << "/" << lTotalPositionsCompared << " correct. " << lPercentCorrect << "% " <<endl;
                lPercentages.push_back( lPercentCorrect );
              }
              else 
                cout << "NO BASIS FOR COMPARISON" << endl;

            } 
          }
        }
      }
      
      cout << "Machine readable list of join percentages." << endl;
      for (int i = 0; i < lPercentages.size(); ++i)
      {
        cout << lPercentages[i] << endl;
      }
      cout << "DONE." << endl;
      delete parser; parser = NULL;
      
      return lJoinedReads.size();
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

