// cerulean.cpp : main project file.

#include "../../lib/khmer.hh"
#include "../../lib/hashtable.hh"
#include "ceruleanfilter.hpp"

#include <iostream>
#include <time.h>

time_t cerulean::global_span_start;
time_t cerulean::last_span_time;
time_t cerulean::last_span_time_before_switch;
time_t cerulean::last_collapse_time;
unsigned int cerulean::collapse_threshold;
short cerulean::basement;

int main(int argc, char *argv[])
{
  unsigned int total_reads;
  unsigned long long n_consumed;

  time_t start, end;
  
  start = time(NULL);
  
  cerulean::global_span_start = time(NULL);
  cerulean::last_span_time = 0;
  cerulean::last_span_time_before_switch = 0;
  cerulean::last_collapse_time = 0;
  cerulean::collapse_threshold = 5000;
  cerulean::basement = 1000;
    
  cerulean::CeruleanFilter bf(atoi(argv[2]), atoi(argv[3]));

  bf.consume_fasta(argv[1], total_reads, n_consumed);

  bf.output_sets();
  
  bf.output_partitioned_file(argv[1], argv[4]);  
  end = time(NULL);
  
  

  std::cout << "DONE: " << difftime(end, start)<< " seconds" << std::endl;

  return 0;
}
