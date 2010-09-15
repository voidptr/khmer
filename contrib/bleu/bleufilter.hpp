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
    private:
      vector<Set*> BackReferences;
      BleuFilter * Parent;
      Set ** Self;
      
    public:
      Set( BleuFilter * aParent )      
      {        
        Parent = aParent;
        BackReferences.push_back( this );
        Self = new Set*();
        *Self = this;
      }
      
      ~Set()
      {
        for ( vector<Set*>::iterator lReference = BackReferences.begin(); lReference != BackReferences.end(); ++lReference )
        {
          assert ( *lReference != this );
          (*lReference) = this;
        }  
      }
      
      unsigned long long size()
      {
        return BackReferences.size();
      }
      
      void consume( Set * aSet )
      { 
        for ( vector<Set*>::iterator lReference = aSet->BackReferences.begin(); lReference != aSet->BackReferences.end(); ++lReference )
        {
          (*lReference) = this;
          BackReferences.push_back( *lReference );
        }       
        
        delete aSet; // this should work.
      }
      
      void Add( HashIntoType aHash )
      {
        
        Parent->_sets[ Parent->HashBinToSetsBin( Parent->HashToHashBin( aHash ) ) ] = Self;
      }
    };
    
  protected:    
    SetID _last_set;

    cBitArray * _hash_table;    
    unsigned int * _hash_bit_counts_lookup; // the number of bits set in the hash table, by every 10k entries

    Set *** _sets; // array of doublepointer set
    HashIntoType _total_unique_hash_count;
    
  public:
    BleuFilter(WordLength ksize, HashIntoType tablesize)
    : Hashtable(ksize, get_first_prime_below( tablesize ))
    { 
      _last_set = 0; // zero == none in use. any number > 0 is a set. _last_set indicates the last set that was allocated. It will always be the largest set number.
            
      _hash_table = new cBitArray( _tablesize );
      _hash_table->Clear();      
      
      _hash_bit_counts_lookup = new unsigned int[(_tablesize / KMER_BIT_COUNT_PARTITION)+1];
      memset(_hash_bit_counts_lookup, 0, ((_tablesize / KMER_BIT_COUNT_PARTITION)+1) * sizeof(unsigned int));

      _sets = NULL; // THIS WILL GET SET LATER
      _total_unique_hash_count = 0; // THIS WILL GET SET LATER  

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
      HashIntoType lHashBin = hash % _tablesize;
     
      _hash_table->Set(lHashBin, true);
      
      n_consumed++;
      
      // now, do the rest of the kmers in this read (add one nt at a time)
      for (unsigned int i = _ksize; i < length; i++) {
        HashIntoType next_hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] );        
        lHashBin = next_hash % _tablesize;
       
        _hash_table->Set(lHashBin, true);
        
        n_consumed++;
      }
      
      return n_consumed;
    }
    
    void prepare_set_arrays()
    {
      _total_unique_hash_count = _hash_table->CountBits();
      
      _sets = new Set**[_total_unique_hash_count];  
      memset(_sets, 0, _total_unique_hash_count * sizeof(Set**));
            
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
      
      
      const char * sp = s.c_str();
      const unsigned int length = s.length();
      unsigned int n_consumed = 0;
      
      HashIntoType forward_hash = 0, reverse_hash = 0;
      
      // generate the hash for the first kmer in the read (fair amount of work)
      HashIntoType hash = _hash(sp, _ksize, forward_hash, reverse_hash);

      HashIntoType lHashBinDeleteMe = HashToHashBin( hash );
      HashIntoType lSetsBinDeleteMe = HashBinToSetsBin( lHashBinDeleteMe );
      Set * lSetDeleteMe = SetsBinToSet( lSetsBinDeleteMe );
      
      // for the first hash
      Set * lWorkingSet = SetsBinToSet( HashBinToSetsBin( HashToHashBin( hash ) ) );

      if ( lWorkingSet != NULL ) // if there's a match //, and it's not a false positive
      {
        // no need to apply SetID. It's already there.
      }
      else // create a new set for us to hold on to 
      {
        lWorkingSet = init_new_set();// otherwise, create a new set, and use that going forward.
        lWorkingSet->Add( hash );
      }
      
      ++n_consumed;
      
      // for the rest of the string
      // now, do the rest of the kmers in this read (add one nt at a time)      
      for (unsigned int i = _ksize; i < length; i++) 
      {
        HashIntoType next_hash = _move_hash_foward( forward_hash, reverse_hash, sp[i] ); 
        
        Set * lEncounteredSet = SetsBinToSet( HashBinToSetsBin( HashToHashBin( next_hash ) ) );
        
        if ( lEncounteredSet != NULL ) // there's a match, and it's not a FP
        {          
          if ( lWorkingSet != lEncounteredSet ) // it's not us!
          {             
            // no need to add to the old set's bins, because this hash is already in the encountered one.
            lWorkingSet = bridge_sets( lEncounteredSet, lWorkingSet );            
          } 
          // if it's us, we can safely skip it (we've been seen before)          
        }
        else // empty set here.
        {
          lWorkingSet->Add( next_hash );
        }
        ++n_consumed;
      }      
      
      return n_consumed;
    }
    
//    bool NotFalsePositive( HashIntoType aHash )
//    {
//      SetID lSetID1 = HashToSetID( aHash );
//      
//      HashIntoType lSetBin2 = HashToSetBin2(aHash);
//      
//      if ( lSetID1 == 0 && _set_IDs_2[ lSetBin2 ] == NULL )
//        return true;
//      
//      if ( _set_IDs_2[ lSetBin2 ] != NULL && 
//          _set_IDs_2[ lSetBin2 ]->find( lSetID1 ) != _set_IDs_2[ lSetBin2 ]->end() )
//        return true;
//      
//      Set * lSet1 = _sets[ lSetID1 ];
//      
//      if ( _set_IDs_2[ lSetBin2 ] != NULL )
//      {
//        for ( set<SetID>::iterator lMatch = _set_IDs_2[ lSetBin2 ]->begin();
//            lMatch != _set_IDs_2[ lSetBin2 ]->end();
//            ++lMatch )
//        {
//          if ( lSet1 == _sets[ *lMatch ] )
//            return true;
//        }
//      }
//      
//      return false;                    
//    }
    
    HashIntoType HashToHashBin ( HashIntoType aHash )
    {
      return (aHash % _tablesize);
    }
   
    HashIntoType HashBinToSetsBin( HashIntoType aBin )
    {
      unsigned long long lBinSectionIndex = (aBin / KMER_BIT_COUNT_PARTITION); // the index of the section before the one we're in
      assert( lBinSectionIndex <= (_tablesize / KMER_BIT_COUNT_PARTITION));
      
      HashIntoType lSetsBin = _hash_table->CountBits( lBinSectionIndex * KMER_BIT_COUNT_PARTITION, aBin );      
      if ( lBinSectionIndex > 0 )      
        lSetsBin += _hash_bit_counts_lookup[ lBinSectionIndex - 1 ];
            
      lSetsBin -= 1; // to make it an array index rather than a count.
      
      assert( lSetsBin < _total_unique_hash_count );
      
      return lSetsBin;
    }
    
    Set * SetsBinToSet( HashIntoType aSetsBin )
    {
      Set** lSetDoublePointer = _sets[ aSetsBin ];
      
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
        return aEncounteredSet;
        aEncounteredSet->consume( aOriginatingSet );
      }
    }
    
    HashIntoType _move_hash_foward( HashIntoType & aOldForwardHash, HashIntoType & aOldReverseHash, const char & aNextNucleotide )
    {
      char lTwoBit = twobit_repr( aNextNucleotide );
      
      aOldForwardHash = aOldForwardHash << 2; // left-shift the previous hash over
      aOldForwardHash |= lTwoBit; // 'or' in the current nucleotide
      aOldForwardHash &= bitmask; // mask off the 2 bits we shifted over.
      
      // now handle reverse complement
      aOldReverseHash = aOldReverseHash >> 2;
      aOldReverseHash |= (compl_twobit(lTwoBit) << (_ksize*2 - 2));
      
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

          Set * lSet = SetsBinToSet( HashBinToSetsBin(HashToHashBin(hash)));
          
          lReadCounts[ lSet ]++;
          
          outfile << ">" << read.name << "\t" << lSet
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

