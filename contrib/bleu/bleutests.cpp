// bleutests.cpp : main project file.

#include "../../lib/khmer.hh"
#include "../../lib/hashtable.hh"
#include "bleufilter.hpp"
//#include "hashtable.hh"

#include <iostream>
#include <time.h>

using namespace std;

int main(int argc, char *argv[])
{

  cout << "Tests!!" << endl;

  cout << "K=16, 4000000007 Tests, " << endl;
  cout << "Symmetric (AAAA <-> TTTT; CCCC <-> GGGG)" << endl;
  cout << "Expected Result - 0 unhomed, 2 sets: " << endl;
  
  bleu::BleuFilter * bf16_1 = new bleu::BleuFilter(16, 4000000007ULL);
  bf16_1->consume_strings_for_hash_table("symmetric.fa");
  bf16_1->deallocate_hash_table_preliminary();
  bf16_1->populate_hash_table_bit_count_lookups();
  bf16_1->allocate_set_offset_table();
  bf16_1->generate_sets("symmetric.fa");
  bf16_1->output_partitioned_file("symmetric.fa", "symmetric_k16_out.fa");
  delete bf16_1;
  cout << endl << endl;

  cout << "K=32, 4000000007 Tests, " << endl;
  cout << "Symmetric (AAAA <-> TTTT; CCCC <-> GGGG)" << endl;
  cout << "Expected Result - 0 unhomed, 2 sets: " << endl;
  
  bleu::BleuFilter * bf16_2 = new bleu::BleuFilter(32, 4000000007ULL);
  bf16_2->consume_strings_for_hash_table("symmetric.fa");
  bf16_2->deallocate_hash_table_preliminary();
  bf16_2->populate_hash_table_bit_count_lookups();
  bf16_2->allocate_set_offset_table();
  bf16_2->generate_sets("symmetric.fa");
  bf16_2->output_partitioned_file("symmetric.fa", "symmetric_k32_out.fa");
  delete bf16_2;
  cout << endl << endl;

  cout << "K=16, 4000000007 Tests, " << endl;
  cout << "Symmetric to K=16 @ k16" << endl;
  cout << "Expected Result - 0 unhomed, 1 set: " << endl;
  
  bleu::BleuFilter * bf16_3 = new bleu::BleuFilter(16, 4000000007ULL);
  bf16_3->consume_strings_for_hash_table("symmetric-to-k16.fa");
  bf16_3->deallocate_hash_table_preliminary();
  bf16_3->populate_hash_table_bit_count_lookups();
  bf16_3->allocate_set_offset_table();
  bf16_3->generate_sets("symmetric-to-k16.fa");
  bf16_3->output_partitioned_file("symmetric-to-k16.fa", "symmetric-to-k16_k16_out.fa");
  delete bf16_3;
  cout << endl << endl;
  
  cout << "K=17, 4000000007 Tests, " << endl;
  cout << "Symmetric to K=16 - @ k17" << endl;
  cout << "Expected Result - 1 unhomed, 1 set: " << endl;
  
  bleu::BleuFilter * bf16_4 = new bleu::BleuFilter(17, 4000000007ULL);
  bf16_4->consume_strings_for_hash_table("symmetric-to-k16.fa");
  bf16_4->deallocate_hash_table_preliminary();
  bf16_4->populate_hash_table_bit_count_lookups();
  bf16_4->allocate_set_offset_table();
  bf16_4->generate_sets("symmetric-to-k16.fa");
  bf16_4->output_partitioned_file("symmetric-to-k16.fa", "symmetric-to-k16_k17_out.fa");
  delete bf16_4;
  cout << endl << endl;

  cout << "K=19, 4000000007 Tests, " << endl;
  cout << "Random A, K=19" << endl;
  cout << "Expected Result - 0 unhomed, 1 set: " << endl;
  
  bleu::BleuFilter * bf19_A = new bleu::BleuFilter(19, 4000000007ULL);
  bf19_A->consume_strings_for_hash_table("random-20-a.fa");
  bf19_A->deallocate_hash_table_preliminary();
  bf19_A->populate_hash_table_bit_count_lookups();
  bf19_A->allocate_set_offset_table();
  bf19_A->generate_sets("random-20-a.fa");
  bf19_A->output_partitioned_file("random-20-a.fa", "random-20-a_k19_out.fa");
  delete bf19_A;
  cout << endl << endl;


  cout << "K=19, 4000000007 Tests, " << endl;
  cout << "Random B, K=19" << endl;
  cout << "Expected Result - 0 unhomed, 1 set: " << endl;
  
  bleu::BleuFilter * bf19_B = new bleu::BleuFilter(19, 4000000007ULL);
  bf19_B->consume_strings_for_hash_table("random-20-b.fa");
  bf19_B->deallocate_hash_table_preliminary();
  bf19_B->populate_hash_table_bit_count_lookups();
  bf19_B->allocate_set_offset_table();
  bf19_B->generate_sets("random-20-b.fa");
  bf19_B->output_partitioned_file("random-20-b.fa", "random-20-b_k19_out.fa");
  delete bf19_B;
  cout << endl << endl;


  cout << "K=20, 4000000007 Tests, " << endl;
  cout << "Random A, K=20" << endl;
  cout << "Expected Result - 99 unhomed, 0 sets: " << endl;
  
  bleu::BleuFilter * bf20_A = new bleu::BleuFilter(20, 4000000007ULL);
  bf20_A->consume_strings_for_hash_table("random-20-a.fa");
  bf20_A->deallocate_hash_table_preliminary();
  bf20_A->populate_hash_table_bit_count_lookups();
  bf20_A->allocate_set_offset_table();
  bf20_A->generate_sets("random-20-a.fa");
  bf20_A->output_partitioned_file("random-20-a.fa", "random-20-a_k20_out.fa");
  delete bf20_A;
  cout << endl << endl;


  cout << "K=20, 4000000007 Tests, " << endl;
  cout << "Random B, K=20" << endl;
  cout << "Expected Result - 99 unhomed, 0 sets: " << endl;
  
  bleu::BleuFilter * bf20_B = new bleu::BleuFilter(20, 4000000007ULL);
  bf20_B->consume_strings_for_hash_table("random-20-b.fa");
  bf20_B->deallocate_hash_table_preliminary();
  bf20_B->populate_hash_table_bit_count_lookups();
  bf20_B->allocate_set_offset_table();
  bf20_B->generate_sets("random-20-b.fa");
  bf20_B->output_partitioned_file("random-20-b.fa", "random-20-b_k20_out.fa");
  delete bf20_B;
  cout << endl << endl;
  

  return 0;
}



