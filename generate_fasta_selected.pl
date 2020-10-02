#!/usr/bin/perl

$head = `head -1 covid_variants_as_genotypes.tdt`;
chomp $head;

($id, @header) = split(/\t/, $head);


foreach $var(@header)
{
print "$var\n";

$vars{$var}++;}

open(OUT, ">representative_genomes.fa");

open(FASTA, "cog_2020-05-08_alignment.fasta");
while(<FASTA>)
{
if(/>(.*)/){$id = $1; $print = 0; $id =~ s/[^A-Z\-a-z0-9]/_/g; if($vars{$id}>0){$print = 1;}}
if($print ==1){print OUT;}
}



