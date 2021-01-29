#!/usr/bin/perl

#The maximum number of differences between two lineages for them to be considred identical and duplicates removed.
$max_diff = 0;

$file = $ARGV[0];

open(LOG, ">>$file.pipeline.log");
$date = `date`; chomp $date;

chomp $file;
print LOG "$date\nRunning process_covid_sequences.pl ${file}_covid_variants.tdt\n";


open(IN, "${file}_covid_variants.tdt");
open(OUT, ">${file}_covid_variants_duplicates_removed.tdt");
$head = <IN>;
chomp $head;
($id, @vars) = split(/\t/, $head);

#Read the header from the input file and set up an index for header names to position
$i = 0; foreach $var(@vars){$vars{$var}++; $var2pos{$var} = $i; $pos2var{$i} = $var; $i++;}

$n = keys(%vars);
print LOG "$n varieties present in header of non-dereplicated input file\n";

while(<IN>)
{
$i++;
chomp;

($id, $bases) = split(/\t/, $_);
@bases = split(//, $bases);
$len = @bases;

#Set up an array with sequence data for each variety column $j.
$j = 0; foreach $base(@bases){$array[$j].= $base; $j++;}
push @posinfo, $id;
}
$total = $i * $len;

print LOG "Loaded $i genotype rows with $len columns of variety data ($total datapoints)\n";
close IN;

#`rm ${file}_covid_variants.tdt`;

#The following section compares sequences and removes duplicates.  Only proper basecalls are considered in the comparison: gaps are ignored.
$i = 0; $j = 0;
foreach $coord(@array)
{
print "Processing variety $i of $len\n"; 
$len = @array;
$seq = $array[$i];
@sequence1 = split(//, $seq);
$j = $i+1;
 while($j < $len)
  {
  $seq2 = $array[$j]; $j++; @sequence2 = split(//, $seq2);
  $m = 0;
  $x = @sequence1; $diff = 0;
  while($m < $x){if($sequence1[$m] ne $sequence2[$m] && $sequence1[$m] =~ /[CATG]/ && $sequence2[$m] =~ /[CATG]/){$diff++; last;}$m++;}
  $var = $vars[$j]; if($diff <= $max_diff){delete $vars{$var};}
  }
$i++;
}

$len = keys %vars;
print LOG "$len varieties left after removal of duplicates\n";

$i = 0;
foreach $base(@posinfo){print OUT "\t$base";}
print OUT "\n";
foreach $var(keys %vars)
{
$coord = $var2pos{$var};
$seq = $array[$coord];
@seq = split(//, $seq);
$seq = join("\t", @seq);
print OUT "$var\t$seq\n";
print "$i $var\n"; $i++;
}

