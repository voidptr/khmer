cdef extern from "bleufilter.hpp":
    ## bleufilter class
    
    ctypedef void (*consume_fasta)(char *, int, int)
    ctypedef void (*populate_hash_table_bit_count_lookups)()
    ctypedef void (*allocate_has_set_table)()
    ctypedef void (*characterize_reads)(char *, int, int)
    ctypedef void (*populate_has_set_bit_count_lookups)()
    ctypedef void (*allocate_set_offset_table)()
    ctypedef void (*generate_sets)(char *, int, int)  
    ctypedef void (*output_partitioned_file)(char *, char *)  
    ctypedef struct c_bleufilter "BleuFilter":
        consume_fasta consume_fasta
        populate_hash_table_bit_count_lookups populate_hash_table_bit_count_lookups
        allocate_has_set_table allocate_has_set_table
        characterize_reads characterize_reads
        populate_has_set_bit_count_lookups populate_has_set_bit_count_lookups
        allocate_set_offset_table allocate_set_offset_table
        generate_sets generate_sets
        output_partitioned_file output_partitioned_file
        
    c_bleufilter * new_bleu_filter "new BleuFilter" (int kmer_len, long long total_table_size)

## END Pyrex Definitions

cdef class bleufilter_container:
    cdef c_bleufilter *this_bleufilter
    
    cdef public int reads_consumed, kmers_consumed
  
    def __cinit__(self, kmer_length_string, total_table_size_string):
        self.this_bleufilter = NULL
        self.this_bleufilter = new_bleu_filter( PyInt_FromString(kmer_length_string), PyLong_FromString(total_table_size_string) )
        
    def consume_fasta(self, path):
        self.this_bleufilter.consume_fasta( path, self.reads_consumed, self.kmers_consumed )
        
## END Python Classes
      
  