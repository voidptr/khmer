#!/usr/bin/perl

print "0\n";
`./dataset_generator 1 1000 900 70 2000 0 > ../../../Test\\ Data/200kreads-100genomes-1klength-70bpreads-generated.fa`;

for (1..99)
{
  print "$_\n";
  print `./dataset_generator 1 1000 900 70 2000 $_ >> ../../../Test\\ Data/200kreads-100genomes-1klength-70bpreads-generated.fa`;
}
print "done.\n";
