#!/usr/bin/perl

$upstream = 50; $downstream = 50;

$sequence = $ARGV[0];
$markers = $sequence;
$markers =~ s/\.fasta/_minimal_markers.txt/;
$markers =~ s/\.fa/_minimal_markers.txt/;




#Read metadata and split input sequence into regions.
$meta = $ARGV[1];
$startweek = $ARGV[2];
$endweek = $ARGV[3];

open(META, "$meta");

$head = <META>;

while($line = <META>)
{
chomp;

#Header columns are:
#sequence_name,country,adm1,sample_date,epi_week,lineage,lineage_support
#e.g. England/BRIS-121A64/2020,UK,UK-ENG,2020-03-13,11,B,90.0

#COG changed the file format - Jan 2021 - updated this REGEX
($seq_name, $country,$adm1,$TrueFalse,$sample_date,$epi_week,$lineage,$lineage_support) = split(/\,/, $line);
if($seq_name =~ /[^\/]+\/([^\/]+)\-[^\/]+.*/){$city = $1;}
 
  elsif ($seq_name =~ /[^\/]+\/([^\/]+)[^\/]+.*/){$city = $1;} 
 else{die "Sequence name $seq_name did not match REGEX - check on file $meta\n";}
if($lineage !~ /[A-Za-z0-9]/){$lineage = "Unknown";}
if($sample_date =~ /(20\d\d-\d\d)-\d\d/){$month = $1;} else{die "Can'd parse date $sample_date in\n$line\n"; }
$directory = "${city}-$month";
 $seq_name =~ s/[^A-Za-z0-9\-]/_/g;
 $id2lineage{$seq_name} = $lineage;
 $id2week{$seq_name} = $epi_week;
}



open(MARKERS, "$markers");
open(OUT, ">annotated_$markers");
$head = <MARKERS>;
chomp $head;
($power, $id, @head) = split(/\t/, $head);
$head = join("\t", @head);
foreach $cell(@head){$lineage = $id2lineage{$cell}; push @lineage, $lineage; $lines{$lineage}++;}
$lineages = join("\t", @lineage);

print OUT "$Id\t$head\n$id\t$lineages\n";

 while(<MARKERS>)
 {
 chomp;
 ($power, $id, @data) = split(/\t/, $_);
 push @ids,$id;
 $data = join("\t", @data);
 if($id =~ /Base_(\d+)_([CATG])\d+([CATG])\d+/){$usedbases{$1}++; $call_lookup{$1}{$2}= 0; $call_lookup{$1}{$3}=2;}
 $i++;
 $j=0;
 foreach $cell(@data){$product[$j].= $cell; $j++; }
 print OUT "$id\t$data\n";
 }
close OUT;
close MARKERS;



open(GEN, ">genotypes_$markers");
print GEN "Sample\tLineage";
foreach $base(sort{$a<=>$b} keys %usedbases){print GEN "\t$base";}
print GEN "\n";

open(SEQ, "$sequence");
while(<SEQ>)
{
chomp;
if(/>(.*)/){$id = $1; $id =~ s/[^0-9A-Z-a-z]/_/g;} else{$id2seq{$id}.=$_;}
}

foreach $id(keys %id2seq)
{
if($id =~ /_([A-Z]+)-/){$region = $1;}
$lineage = $id2lineage{$id};
if($lineage !~ /Unknown/i && $lineage =~ /\D/)
{
$seq = $id2seq{$id};

foreach $pos(keys %usedbases)
{
$start = $pos - $upstream;
$start--;
$end = $pos+ $downstream;
$end--;
$len = $end-$start;
$subseq = substr($seq,$start,$len);
$variantstring{$pos}{$subseq}++;
}


$i = 1;
@seq = split(//, $seq);
$genotype = ""; $binary = ""; 
foreach $base(@seq){if($usedbases{$i}>0){$genotype.= $base; $binary.= $call_lookup{$i}{$base};} $i++;}
if($genotype =~ /^[CATG]+$/)
{
print GEN "$id\t$lineage\t$genotype\t$binary\n";
$genotype2lineage{$genotype}{$lineage}++;
$genotype2region{$genotype}{$region}++;
$regions{$region}++;

$genotype2binary{$genotype} = $binary;
$genotype_count{$genotype}++;
$week = $id2week{$id};
if($week >= $startweek && $week <= $endweek){push @genotype_array, $genotype;}
}
}
}

%id2seq =();

#Shuffle the array of all haploypes (of used SNPs) and take pairs (they will be random due to shuffle) to see if they're identical
open(DIFF, ">proportion_samples_different_at_selected_SNP_loc.txt");

shuffle(\@genotype_array);

$len = @genotype_array/2;
$i = -2; 
$identical = 0;
$different = 0;
while($i <$len)
{
$i+=2;
$j = $i+1;
if($genotype_array[$i] eq $genotype_array[$j]){$identical++;}
}
$prop_different = $identical/$len;
print DIFF "$identical of $len random sample pairs ($prop_different) were identical\n"; 


open(SUMMARY, ">summary_$markers");
foreach $pos(sort {$a<=>$b} keys %usedbases){print SUMMARY "\t$pos";}

print SUMMARY "Genotpye\tbinary\tfrequency\tn_lineages\tLineages\n";

foreach $genotype(keys %genotype2lineage)
{
$binary = $genotype2binary{$genotype};
$freq = $genotype_count{$genotype};
print SUMMARY "$genotype\t$binary\t$freq";
$ref = $genotype2lineage{$genotype};
%hash = %$ref;
$n = keys %hash; print SUMMARY "\t$n\t";
foreach $lineage(sort keys %hash){print SUMMARY "$lineage ";}
print SUMMARY "\n";
}


open(REGION, ">regional_$markers");

foreach $region(sort keys %regions){print REGION "\t$region"; push @regions, $region;}
print REGION "\n";

foreach $genotype( sort keys %genotype2region)
{
print REGION "$genotype";

foreach $region(@regions){$sum = 0; $sum+= $genotype2region{$genotype}{$region};print REGION "\t$sum";}
print REGION "\t";
foreach $lineage (sort keys %lines)
  {
  $n = $genotype2lineage{$genotype}{$lineage};
  if($n >0){  print REGION " $lineage";}
  }
print REGION "\n";
}

`rm subsequence*.txt`;
open(AMINO, ">aminoacids_$markers");
foreach $pos(keys %variantstring)
{
$max = 0;
$maxstring = "";
open(VAR, ">subsequence_${pos}_$markers");
$ref = $variantstring{$pos};
%hash = %$ref;
foreach $substr(sort keys %hash)
{
$edited = $substr;
$ref2 = $call_lookup{$pos};
%hash2 = %$ref2;
$variant = ""; 
foreach $base(sort keys %hash2){$variant.= $base;}  
$variant =~ s/([A-Z])([A-Z])/\[$1\/$2\]/;
substr($edited,$upstream,1,$variant);
#print VAR "$old\n"; debug line

$i = $variantstring{$pos}{$substr};
 if($i >1)
 {
 if($i >$max){$max = $i; $maxstring = $edited;}
 print VAR "$edited\t$i\n";
 }
}
print "Getting amino acid variants for $maxstring\n"; 
$seq = $maxstring;
$seq2 = $maxstring;
$seq =~ s/\[([A-Z])\/[A-Z]\]/$1/;
$seq2 =~ s/\[[A-Z]\/([A-Z])\]/$1/;
print "$seq\n$seq2\n\n";
$amino1 = `./ORF_find.pl $seq`;
$amino2 = `./ORF_find.pl $seq2`;
chomp $amino1; chomp $amino2;
print AMINO "$pos\n$seq\n$amino1\n$seq2\n$amino2\n\n";
close VAR;
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
