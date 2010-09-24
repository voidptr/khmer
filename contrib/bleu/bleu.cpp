// bleu.cpp : main project file.

#include "../../lib/khmer.hh"
#include "../../lib/hashtable.hh"
#include "bleufilter.hpp"

#include <iostream>
#include <time.h>

int main(int argc, char *argv[])
{
  unsigned int total_reads;
  unsigned long long n_consumed;

  time_t start, end;
  
  start = time(NULL);
        
  bleu::BleuFilter bf(atoi(argv[2]), atoll(argv[3]));

  bf.consume_fasta(argv[1], total_reads, n_consumed);
  bf.populate_hash_table_bit_count_lookups();
  bf.allocate_has_set_table();

  bf.characterize_reads(argv[1], total_reads, n_consumed);
  bf.populate_has_set_bit_count_lookups();
  bf.allocate_set_offset_table();

  bf.generate_sets(argv[1], total_reads, n_consumed);
  
  bf.output_partitioned_file(argv[1], argv[4]);  
  end = time(NULL);
  
  

  std::cout << "DONE: " << difftime(end, start)<< " seconds" << std::endl;

  return 0;
}



