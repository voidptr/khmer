#define VERSION "0.3"

#define MAX_COUNT 255

#include "cBitArray.h"

#include "hashtable.hh"

namespace khmer {
  // largest number we can count up to, exactly. (8 bytes)
  typedef unsigned long long int ExactCounterType;

  // largest number we're going to hash into. (8 bytes/64 bits/32 nt)
  typedef unsigned long long HashIntoType;
  typedef cBitArrayBytes HashIntoType_Big;
  
  // the base type to represent an individual nucleotide. 
  // Since we can't use a pair of bits, this is as close as it gets.
  typedef unsigned char NucleotideType;

  // largest size 'k' value for k-mer calculations.  (1 byte/255)
  typedef unsigned char WordLength;

  // largest number we can count up to, approximately. (8 bytes/127).
  // see MAX_COUNT, above.
  typedef unsigned char BoundedCounterType;

  typedef void (*CallbackFn)(const char * info, void * callback_data,
			     unsigned int n_reads, unsigned long long other);

};
