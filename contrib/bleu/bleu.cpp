// bleu.cpp : main project file.

#include "../../lib/khmer.hh"
#include "../../lib/hashtable.hh"
#include "bleufilter.hpp"

#include <iostream>
#include <time.h>

time_t bleu::global_span_start;
time_t bleu::last_span_time;
time_t bleu::last_span_time_before_switch;
time_t bleu::last_collapse_time;
unsigned int bleu::collapse_threshold;
short bleu::basement;

int main(int argc, char *argv[])
{
  unsigned int total_reads;
  unsigned long long n_consumed;

  time_t start, end;
  
  start = time(NULL);
  
  bleu::global_span_start = time(NULL);
  bleu::last_span_time = 0;
  bleu::last_span_time_before_switch = 0;
  bleu::last_collapse_time = 0;
  bleu::collapse_threshold = 5000;
  bleu::basement = 1000;
    
  //char * lBleh = argv[3];
  //unsigned int lBleh2 = atol( lBleh );
  
  //unsigned long long lBleh3 = atoll( lBleh );
  
 // unsigned long long lVal = 32000000007;
 // unsigned long long lVal =  3200000007;
 //   unsigned long long lVal = 16000000007;
  
  bleu::BleuFilter bf(atoi(argv[2]), atoll(argv[3]));

  bf.consume_fasta(argv[1], total_reads, n_consumed);
  bf.prepare_set_arrays();
  bf.generate_sets(argv[1], total_reads, n_consumed);

  //bf.forecast_memory_consumption();

  
//  bf.output_sets();
  
  bf.output_partitioned_file(argv[1], argv[4]);  
  end = time(NULL);
  
  

  std::cout << "DONE: " << difftime(end, start)<< " seconds" << std::endl;

  return 0;
}



