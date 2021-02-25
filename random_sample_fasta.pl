#!/usr/bin/perl


$min = 10;  #The minumum number of sequences below which the script will die with an error. It can be any value but is present as an input  file check
$file = $ARGV[0];
$sample_size = $ARGV[1];

#Check $file looks like a FASTA file

$id = `head -1 $file`; chomp $id;
$seq =`head -2 $file|tail -1`; chomp $seq;
$seqs = `grep -c ">" $file`; chomp $lines;

$usage = "Usage: ./random_sample_fasta.pl <filename> [number of sequences to get]\n";

if($sample_size !~ /^\d+$/){die "\"$sample_size\" is not a valid sample size, please enter a number\n$usage";}
if($sample_size >= $seqs){die "You asked for $sample_size random sequences but $file only seems to have $seqs\n$usage";}
if($id !~ /^>/){die "$id doesn't look like a FASTA header, check the input file $file\n$usage";} 
if($seq !~ /^[A-Z\-]+$/){die "$seq is not a valid FASTA sequence line, check input file $file\n$usage";}
if($seqs <$min){die "There are only $seqs sequences in $file, so not really worth random sampling\n$usage";}

$i = 0;
while($i < $seqs){$ary[$i] = $i; $i++;}
shuffle(\@ary);
$i = 0; while($i <$sample_size){$linenum = $ary[$i]; $selected{$linenum}++; $i++;}


open(IN, "$file" ||die "Cant open $file");
open(OUT, ">random_${sample_size}_$file" ||die "Can't open outfile random_${sample_size}_$file");

$i = -1;
while($line = <IN>)
{
if($line =~ />/){$i++; if($selected{$i} >0){$print = 1;} else{$print = 0;}}
if($print ==1){print OUT "$line";}
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
