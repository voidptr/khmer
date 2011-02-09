The test files in this directory should behave as follows:

random-20-a.fa & random-20-b.fa - randomly generated genomes that are all one set at k < 20.

At k=19, there is one partition with 99 reads
At k=20, there are 99 unhomed reads, 0 sets.

symmetric.fa - tests AAAA = TTTT; CCCC = GGGG

at k=16, 0 unhomed, 2 sets
at k=32, 0 unhomed, 2 sets.

symmetric-to-k16.fa - tests AAAA = TTTT; CCCC = GGGG, but only until k=16

at k=16, 0 unhomed, 1 set
at k=17, 1 unhomed, 1 set.

200kreads-100genomes-1klength-70bpreads-generated.fa - 200k reads, sampled from 100 (1kb length)
 genomes. 2000 reads per genome. Each read is 70bp, and there are zero sampling errors.

any k, 0 unhomed reads.

k<=12, 1 set
k=13, 28 sets
k=14, 67 sets
k=15, 91 sets
k=16, 98 sets, 0 unhomed.
k>=17, 0 unhomed, 100 sets, 2000 reads per set.
