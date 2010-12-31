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
    cout << argv[0] << " inputfile.fa ksize memoryfootprint outpufile.fa [readsthatjoin.fa]" << endl;
    return 1;
  }
  time_t start, end;

  start = time(NULL);

  bleu::BleuFilter bf(atoi(argv[2]), atoll(argv[3]));

  bf.consume_strings_for_hash_table(argv[1]);
  
  // populate the hash table
  bf.deallocate_hash_table_preliminary();

  // allocate valid permutation table and has_set table
  bf.populate_hash_table_bit_count_lookups();

  bf.allocate_set_offset_table();

  // generate the sets
  if ( argc == 5 )
    bf.generate_sets(argv[1]);
  else // or, also output the join reads.
    bf.output_join_reads(argv[1], argv[5]);

  bf.output_partitioned_file(argv[1], argv[4]);
  end = time(NULL);



  std::cout << "DONE: " << difftime(end, start)<< " seconds" << std::endl;

  return 0;
}



