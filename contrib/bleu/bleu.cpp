// bleu.cpp : main project file.

#include "../../lib/khmer.hh"
#include "../../lib/hashtable.hh"
#include "bleufilter.hpp"
//#include "hashtable.hh"

#include <iostream>
#include <time.h>

using namespace std;

int main(int argc, char *argv[])
{
//  unsigned int total_reads;
//  unsigned long long n_consumed;
//
  if (argc < 5)
  {
    cout << argv[0] << " inputfile.fa ksize memoryfootprint outpufile.fa" << endl;
    return 1;
  }
  time_t start, end;

  start = time(NULL);

  bleu::BleuFilter bf(atoi(argv[2]), atoll(argv[3]));

  bf.consume_strings_for_hash_table(argv[1]);
  
  // populate the hash table
  //  bf.consume_reads(argv[1], total_reads, n_consumed, &bleu::BleuFilter::consume_strings_for_hash_table);
  bf.deallocate_hash_table_preliminary();

  // allocate valid permutation table and has_set table
  bf.populate_hash_table_bit_count_lookups();

  bf.allocate_set_offset_table();

  // generate the sets
  bf.generate_sets(argv[1]);
  //  bf.consume_reads(argv[1], total_reads, n_consumed, &bleu::BleuFilter::consume_strings_for_set);

  bf.output_partitioned_file(argv[1], argv[4]);
  end = time(NULL);



  std::cout << "DONE: " << difftime(end, start)<< " seconds" << std::endl;

  return 0;
}



