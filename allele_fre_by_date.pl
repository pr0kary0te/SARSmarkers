#!/usr/bin/perl


`mkdir week_analysis`;


#primer_designs_for_cog.tdt
$markers = $ARGV[0];
chomp $markers;

open(IN, $markers);

#$head = <IN>;
while(<IN>)
{
chomp;
($marker, $seq) = split(/\t/, $_);
if($marker =~/CoV_(\d+)/){$pos = $1;}
if($seq =~ /\[([A-Z])\/([A-Z])\]/){$lookup{$pos}{$1}++; $lookup{$pos}{$2}++; print "$pos $1 $2\n";}
}


$meta = $ARGV[1];

open(META, $meta);
#sequence_name,country,adm1,sample_date,epi_week,lineage,lineage_support

$head = <META>;
while(<META>)
{
#COG changed the file format - Jan 2021 - updated this REGEX
#($name,$country,$adm,$date,$week,@other) = split(/\,/, $_);
($name,$country,$adm,$truefalse,$date,$week,@other) = split(/\,/, $_);

$name2week{$name} = $week;
}
close META;
print "Loaded metadata\n";

$seq = $ARGV[2];
open(SEQ, "$seq");
while($line =<SEQ>)
{
chomp $line;
if($line =~ />(.*)/){$name = $1; $week = $name2week{$name};}
else
{

foreach $pos(sort keys %lookup)
{
$actual = $pos -1;
$match = substr($line,$actual,1);
#Make sure this is a valid SNP allele:
if($lookup{$pos}{$match} >0){$final{$pos}{$week}{$match}++; $weeks{$week}++; $overall_total{$pos}{$match}++; $genotypes{$name}.=$match;}


}

}
}
$n = keys %genotypes;

print "Loaded Sequence data for $n samples\n";


foreach $pos(keys %final)
{
open(OUT, ">week_analysis/CoV$pos.txt");
print OUT "Week";
$ref = $lookup{$pos};
%hash = %$ref;
foreach $base(sort keys %hash){if($base =~ /^[CATG]$/){print OUT "\t$base";}}
print OUT "\n";

foreach $week (sort {$a<=>$b} keys %weeks)
{
if($week =~ /\d+/){print OUT "$week";}
$total = 0;
foreach $base(sort keys %hash){if($base =~ /^[CATG]$/){print OUT "\t$final{$pos}{$week}{$base}"; $total+= $final{$pos}{$week}{$base};}}
$max = 0; $major_allele = "";
foreach $base(sort keys %hash)
 {
 $n = $final{$pos}{$week}{$base}; $prop = $n/$total;
 $basetotal = $overall_total{$pos}{$base}; 
 if($basetotal>$max){$max = $basetotal; $major_allele = $base; $maxprop = $prop; }print OUT "\t$prop";
 }
$grouped{$pos}{$week}=$maxprop; $major{$pos}=$major_allele;
print OUT "\n";
}
close OUT;
}


open(OUT, ">grouped_week_analysis.txt");
print OUT "Position";
foreach $week (sort {$a<=>$b} keys %weeks){print OUT "\t$week";}
print OUT "\n";

foreach $pos(sort {$a<=>$b} keys %grouped)
{
print OUT "$pos $major{$pos}";
foreach $week (sort {$a<=>$b} keys %weeks){print OUT "\t$grouped{$pos}{$week}";}
print OUT "\n";
}



$n = keys %lookup;

foreach $name(keys %genotypes)
{
$week = $name2week{$name};
$genotype = $genotypes{$name};
$l = length($genotype);
if($n == $l)
 {
 $complete_genotypes{$week}{$genotype}++;
 }
}

print "
To test the discriminatory power of the selected markers over time, the genotypes observed in each week are placed in random order\n
and then adjacent genotypes compared to see if they are differnt. For example, if there are 50 samples sequenced in week 5, and the\n
pipeline produced 24 markers, we find the 24-marker haplotypes present in this week and put them in random order. We then compare the\n
1st versus 2nd, 2nd verus 3rd, 3rd versus 4th haplotype and so on until we reach the 49th versus 50th. There will be n-1 comparisons\n\n\n";  
print "\nWeek\tDifferences\tComparisons\tPercent comparisons different\n";
foreach $week(sort {$a<=>$b} keys %complete_genotypes)
{
$ref = $complete_genotypes{$week};
%hash = %$ref;
@complete_genotypes = ();
$n = keys %hash;
foreach $genotype(keys %hash)
  {
  $cases = $hash{$genotype};
  while($cases >0) { push @complete_genotypes, $genotype; $cases--;}
  }

shuffle (\@complete_genotypes);
$l = @complete_genotypes;
$i = 0;
$j = 1;
$diff = 0;
$comparisons = $l -1;
 while ($j < $l)
 {
 if($complete_genotypes[$i] ne $complete_genotypes[$j]){$diff++;}
 $i++; $j++;
 }
if($week =~ /\d/ && $comparisons >0){$prop = int(100*($diff/$comparisons)); print "$week\t$diff\t$comparisons\t$prop\n";}
}



sub shuffle
{
my $array = shift;
my $i;
for ($i = @$array; --$i;) {
my $j = int(rand($i+1));
next if $i == $j;
@$array[$i,$j] = @$array[$j,$i];
}
}
