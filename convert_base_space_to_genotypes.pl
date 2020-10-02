#!/usr/bin/perl

$file = $ARGV[0];
chomp $file;
open(IN, "${file}_covid_variants_duplicates_removed.tdt");
open(OUT, ">${file}_covid_variants_as_genotypes.tdt");
$head = <IN>;
chomp $head; 
($var,@head) = split(/\t/, $head); 
while(<IN>)
{
chomp;
$j = 0;
($var, @data) = split(/\t/, $_);
push @vars, $var;
foreach $cell(@data){$pos = $head[$j]; $genotype[$j].= $cell; $j++;}
}

close IN;
#`rm ${file}_covid_variants_duplicates_removed.tdt`;

$i = -1;
$head = join("\t", @vars);
print OUT "\t$head\n";
foreach $genotype(@genotype)
{
$i++;
$var = $head[$i];
@data = split(//, $genotype);
$data = join("\t",@data);
if($var =~ / ([CATG])\d+([CATG])\d+/){$major = $1; $minor = $2; $data =~ s/$major/0/g; $data =~ s/$minor/2/g;}
$var =~ s/\s/_/g;
print OUT "$var\t$data\n";
}
