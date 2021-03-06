#! /usr/bin/env python
import sys, khmer

K = 32
N_HT=4
HT_SIZE = int(5e8)

###

filename = sys.argv[1]
output = sys.argv[2]

ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

ht.consume_fasta(filename)

print 'preparing hist...'
z = ht.abundance_distribution(filename)
fp = open(output, 'w')

for n, i in enumerate(z[1:]):
    print >>fp, n + 1, i

#ht.fasta_dump_kmers_by_abundance(filename, 255)
