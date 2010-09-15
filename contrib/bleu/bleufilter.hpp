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

#define HASHES 8

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
    private:
      vector<Set**> BackReferences;
      BleuFilter * Parent;
      
    public: 
      Set ** Self;
      
    public:
      Set( BleuFilter * aParent )      
      {        
        Self = new Set*();
        *Self = this;

        Parent = aParent;
        BackReferences.push_back( Self );
      }
      
      ~Set()
      {
        for ( vector<Set**>::iterator lReference = BackReferences.begin(); lReference != BackReferences.end(); ++lReference )
        {
          assert ( **lReference != this );
        }  
      }
      
      unsigned long long size()
      {
        return BackReferences.size();
      }
      
      void consume( Set * aSet )
      { 
        
        for ( int i = 0; i < aSet->BackReferences.size(); ++i )
        {
          *(aSet->BackReferences[i]) = this;
          BackReferences.push_back( aSet->BackReferences[i] );
        }
        
        delete aSet; // this should work.
      }
     

    };
    
  protected:    
    cBitArray * _hash_table[HASHES]; // two-dimensional hash-table
    unsigned int * _hash_bit_counts_lookup[HASHES]; // the number of bits set in the hash table, by every 10k entries

    Set *** _sets[HASHES]; // two-dimensional array of doublepointer set (really? seriously?)
    HashIntoType _total_unique_hash_count[HASHES];
    
    HashIntoType _tablesizes[HASHES];
    
  public:
    BleuFilter(WordLength ksize, HashIntoType tablesize)
    : Hashtable(ksize, get_first_prime_below( tablesize / HASHES ))
    {       
      _tablesizes[0] = _tablesize;
      for ( int i = 1; i < HASHES; ++i )      
      {
        _tablesizes[i] = get_first_prime_below(_tablesizes[i-1]); 
      }
      
      for ( int j = 0; j < HASHES; ++j )      
      {
        _hash_table[j] = new cBitArray( _tablesizes[j] );
        _hash_table[j]->Clear();      
                    
        _hash_bit_counts_lookup[j] = new unsigned int[(_tablesizes[j] / KMER_BIT_COUNT_PARTITION)+1];
        memset(_hash_bit_counts_lookup[j], 0, ((_tablesizes[j] / KMER_BIT_COUNT_PARTITION)+1) * sizeof(unsigned int));
                    
        _sets[j] = NULL; // THIS WILL GET SET LATER
        _total_unique_hash_count[j] = 0; // THIS WILL GET SET LATER
      }

      // housekeeping
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
      
      for (int i = 0; i < HASHES; ++i )
      {
        HashIntoType lHashBin = hash % _tablesizes[i];
        _hash_table[i]->Set(lHashBin, true);
      }
           
      
      n_consumed++;
      
      // now, do the rest of the kmers in this read (add one nt at a time)
      for (unsigned int i = _ksize; i < length; i++) {
        HashIntoType next_hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );        
        
        for (int i = 0; i < HASHES; ++i )
        {
          HashIntoType lHashBin = next_hash % _tablesizes[i];
          _hash_table[i]->Set(lHashBin, true);
        }
        
        n_consumed++;
      }
      
      return n_consumed;
    }
    
    void prepare_set_arrays()
    {
      
      for (int i = 0; i < HASHES; ++i )
      {
        _total_unique_hash_count[i] = _hash_table[i]->CountBits();
        
        _sets[i] = new Set**[_total_unique_hash_count[i]];  
        memset(_sets[i], 0, _total_unique_hash_count[i] * sizeof(Set**));
        
        
        unsigned long long lLookupTableSize = (_tablesizes[i] / KMER_BIT_COUNT_PARTITION)+1;
        
        for (unsigned long long j = 0; j < lLookupTableSize; ++j)
        {      
          unsigned long long lSectionStartIndex = j * KMER_BIT_COUNT_PARTITION;
          unsigned long long lSectionStopIndex = ((j + 1) * KMER_BIT_COUNT_PARTITION) - 1;
          if ( lSectionStopIndex >= _tablesizes[i] )
            lSectionStopIndex = _tablesizes[i] - 1;
          
          // temporarily store the single section count. (this section may be problematic, since ..lookup[i]'s value is a pointer, so what does [j] do?
          _hash_bit_counts_lookup[i][j] = _hash_table[i]->CountBits(lSectionStartIndex, lSectionStopIndex);
          
          if ( i > 0 ) // apply the summation
          {
            _hash_bit_counts_lookup[i][j] += _hash_bit_counts_lookup[i][j-1];
          }
        }  
  
      }
      
            
      
            
       
      
      cout << "DONE PREP" << endl;
    }
    
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
      //cout << endl << s << " " << s.length() << endl;
      
      const char * sp = s.c_str();
      const unsigned int length = s.length();
      unsigned int n_consumed = 0;
      
      HashIntoType forward_hash = 0, reverse_hash = 0;
      
      // generate the hash for the first kmer in the read (fair amount of work)
      HashIntoType hash = _hash(sp, _ksize, forward_hash, reverse_hash);
      
      // for the first hash
      
      Set * lWorkingSet = SelectOrCreateProperSet( hash );
      
      ++n_consumed;
      
      // for the rest of the string
      // now, do the rest of the kmers in this read (add one nt at a time)      
      for (unsigned int i = _ksize; i < length; i++) 
      {
        HashIntoType next_hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] ); 
        
        lWorkingSet = AssignOrBridgeToProperSet( next_hash, lWorkingSet );
        
        ++n_consumed;
      }      
      
      return n_consumed;
    }
    
    Set * SelectOrCreateProperSet( HashIntoType aHash )
    {
      Set * lProspectiveSets[HASHES];
      HashIntoType lProspectiveSetBins[HASHES];
      for (int i = 0; i < HASHES; ++i)
      {
        lProspectiveSetBins[i] = HashBinToSetsBin( HashToHashBin( aHash, i ), i );
        lProspectiveSets[i] = SetsBinToSet( lProspectiveSetBins[i], i );
      }
      
      Set * lSet = NULL;
      for (int i = 0; i < HASHES; ++i)
      {
        if ( lProspectiveSets[i] == NULL ) // we found an empty slot! Woohoo!
        {
          if ( lSet == NULL ) // put a new set there.
            lSet = init_new_set();
          
          _sets[i][ lProspectiveSetBins[i] ] = lSet->Self;
        }
      }
      if ( lSet != NULL )
        return lSet;
      
      // So, none of them were empty. Perform pairwise comparisons of the hashes to see which agree.
      map<Set *, int> lRepresented;
      for (int i = 0; i < HASHES; ++i)
      {
        lRepresented[ lProspectiveSets[i] ]++;
      }
      
      Set * lMostRepresented = NULL; // because maps are sorted, this will be the set with the lowest pointer, and the highest count.
      for ( map<Set*, int>::iterator lSet = lRepresented.begin(); lSet != lRepresented.end(); ++lSet )
      {
        if ( lMostRepresented == NULL || lSet->second > lRepresented[ lMostRepresented ] )
          lMostRepresented = lSet->first;
      }
      
      int lMaxCount = lRepresented[ lMostRepresented ];
      
      if ( lRepresented[ lMostRepresented ] > 1 ) // woohoo, good enough.
        ;//cout << "FOUND AN ACTUAL MATCH: " << lMaxCount << endl;
      else 
        cout << "GOOD ENOUGH. FINE. JOINING A FP SET." << endl;

      
      return lMostRepresented;
    }
    
    Set * AssignOrBridgeToProperSet( HashIntoType aHash, Set * aWorkingSet )
    {
      Set * lProspectiveSets[HASHES];
      HashIntoType lProspectiveSetBins[HASHES];
      for (int i = 0; i < HASHES; ++i)
      {
        lProspectiveSetBins[i] = HashBinToSetsBin( HashToHashBin( aHash, i ), i );
        lProspectiveSets[i] = SetsBinToSet( lProspectiveSetBins[i], i );
      }
      
      Set * lSet = NULL;
      for (int i = 0; i < HASHES; ++i)
      {
        if ( lProspectiveSets[i] == NULL ) // we found an empty slot, so we can't possibly be like the others that might be around! Woohoo!
        {
          if ( lSet == NULL ) // put a new set there.
            lSet = aWorkingSet;          
          _sets[i][ lProspectiveSetBins[i] ] = lSet->Self;
        }
      }
      if ( lSet != NULL )
        return lSet;
      
      // so, none of them were empty. were any of them us?
      map<Set *, int> lRepresented;
      
      lRepresented.insert( map<Set*, int>::value_type(aWorkingSet, 1) );
      
      for (int i = 0; i < HASHES; ++i)
      {
        lRepresented[ lProspectiveSets[i] ]++;
      }
      Set * lMostRepresented = NULL; // because maps are sorted, this will be the set with the lowest pointer, and the highest count.
      for ( map<Set*, int>::iterator lSet = lRepresented.begin(); lSet != lRepresented.end(); ++lSet )
      {
        if ( lMostRepresented == NULL || lSet->second > lRepresented[ lMostRepresented ] )
          lMostRepresented = lSet->first;
      }
      
      int lMaxCount = lRepresented[ lMostRepresented ];
      
      if ( lRepresented[ lMostRepresented ] > 1 ) // woohoo, good enough.
        ;//cout << "FOUND AN ACTUAL MATCH: " << lMaxCount << endl;
      else 
        cout << "GOOD ENOUGH. FINE. JOINING A FP SET." << endl;
      
      if ( lMostRepresented != aWorkingSet )
        return bridge_sets( lMostRepresented, aWorkingSet);
      else
        return aWorkingSet;
    } 
    HashIntoType HashToHashBin ( HashIntoType aHash, int i )
    {
      return (aHash % _tablesizes[i]);
    }
   
    HashIntoType HashBinToSetsBin( HashIntoType aBin, int i )
    {
      unsigned long long lBinSectionIndex = (aBin / KMER_BIT_COUNT_PARTITION); // the index of the section before the one we're in
      assert( lBinSectionIndex <= (_tablesizes[i] / KMER_BIT_COUNT_PARTITION));
      
      HashIntoType lSetsBin = _hash_table[i]->CountBits( lBinSectionIndex * KMER_BIT_COUNT_PARTITION, aBin );      
      
      if ( lBinSectionIndex > 0 )      
        lSetsBin += _hash_bit_counts_lookup[i][ lBinSectionIndex - 1 ];
            
      lSetsBin -= 1; // to make it an array index rather than a count.
      
      return lSetsBin;
    }
    
    Set * SetsBinToSet( HashIntoType aSetsBin, int i )
    {
      Set** lSetDoublePointer = _sets[i][ aSetsBin ];
      
      if ( lSetDoublePointer == NULL )
        return NULL;
      
      return *lSetDoublePointer;
    }
    
    Set * init_new_set()
    {
      return new Set( this );            
    }
    
    Set * bridge_sets( Set * aEncounteredSet, Set * aOriginatingSet )
    { 
      assert ( aOriginatingSet != aEncounteredSet );

      if ( aOriginatingSet->size() > aEncounteredSet->size() )
      {
        aOriginatingSet->consume( aEncounteredSet );
        return aOriginatingSet;
      }
      else
      {
        aEncounteredSet->consume( aOriginatingSet );
        return aEncounteredSet;
      }
    }
    
    HashIntoType _move_hash_foward( HashIntoType & aOldForwardHash, HashIntoType & aOldReverseHash, const char & aNextNucleotide )
    {
      //cout << " -sF- " << aOldForwardHash << " " << _revhash(aOldForwardHash, 32) << endl;
      //cout << " -sR- " << aOldReverseHash << " " << _revhash(aOldReverseHash, 32) << endl;
      
      HashIntoType lTwoBit = twobit_repr( aNextNucleotide );
      
      aOldForwardHash = aOldForwardHash << 2; // left-shift the previous hash over
      aOldForwardHash |= lTwoBit; // 'or' in the current nucleotide
      aOldForwardHash &= bitmask; // mask off the 2 bits we shifted over.
      
      // now handle reverse complement
      aOldReverseHash = aOldReverseHash >> 2;
      aOldReverseHash |= (compl_twobit(lTwoBit) << (_ksize*2 - 2));
      
      //cout << " -eF- " << aOldForwardHash << " " << _revhash(aOldForwardHash, 32) << endl;
      //cout << " -eR- " << aOldReverseHash << " " << _revhash(aOldReverseHash, 32) << endl;      
      
      //cout << endl;
      
      // pick the better bin of the forward or reverse hashes
      return uniqify_rc(aOldForwardHash, aOldReverseHash);
    }
    
//    void output_sets()
//    {
//      for ( map<unsigned int, Set*>::iterator lIt = _sets.begin(); lIt != _sets.end(); ++lIt )
//      {
//        mUniqueSets.insert( lIt->second );
//      }
//            
//      cout << setw(6) << "unique set count: "<< mUniqueSets.size() << endl;
//      cout << endl;
//    }
    
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
      
      map<Set*, unsigned int> lReadCounts;
      
      while(!parser->is_complete()) {
        read = parser->get_next_read();
        seq = read.seq;
        
        if (check_read(seq)) {
          first_kmer = seq.substr(0, _ksize);
          
          // generate the hash for the first kmer in the read (fair amount of work)
          HashIntoType hash = _hash(first_kmer.c_str(), _ksize, forward_hash, reverse_hash);

          /// TODO FIX PER WHATEVER THE ABOVE ALGO IS SUPPOSED TO BE.
          //assert (0);
          Set * lSet = SelectOrCreateProperSet( hash );
          //Set * lSet = SetsBinToSet( HashBinToSetsBin(HashToHashBin(hash)));
          
          lReadCounts[ lSet ]++;
          
          outfile << ">" << read.name << "\t" 
          //<< lSet FIX ME
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
      
      for ( map<Set *, unsigned int>::iterator lIt = lReadCounts.begin(); lIt != lReadCounts.end(); ++lIt )
      {
        cout << setw(10) << lIt->first;
        cout << setw(10) << lIt->second << endl;        
      }
      
      cout << setw(6) << "unique set count: "<< lReadCounts.size() << endl;
      cout << endl;
      
      cout << "Hash Entries: " << _total_unique_hash_count << endl;
      
      delete parser; parser = NULL;
      
      return lReadCounts.size();
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

#endif //BLEUFILTER_HPP

