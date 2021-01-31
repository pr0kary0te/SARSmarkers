#!/usr/bin/perl

open(VARS, "priority_variants_pos.txt");

$head = <VARS>;
chomp $head;
@head = split(/\t/, $head);


while(<VARS>)
{
#Variant Gene    subs    weight
#B.1.1.7 nsp2    C913T   0
chomp;
($variant, $gene, $call, $weight) = split(/\t/, $_);
if($call =~ /\D(\d+)(\D)/)
   {
   $snp = "$2"; $pos = $1; if($variant =~/\D/ && $pos =~ /\d/ && $snp =~ /[A-Z]/){$lookup{$variant}{$pos}= $snp;} else{die "Bad line in input file $_";}
   }
elsif($call -~ /(\d+)-(\d+)/)
   {
   $startdel = $1;
   $enddel = $2;
   while($startdel <= $enddel){$lookup{$variant}{$startdel}= "-"; $startdel++;}
   }


}


close VARS;

open(REF, "wuhan_hu-1.fa");
$head = <REF>;
while(<REF>){chomp; $wuhan_reference.=$_;}
close REF;

foreach $strain(keys %lookup)
{
open(OUT, ">$strain.fa");
$seq = $wuhan_reference;

$ref = $lookup{$strain};
%hash = %$ref;
foreach $pos(sort {$a<=>$b} keys %hash)
   {
   $snp = $hash{$pos};
   substr($seq, $pos-1, 1) = $snp;      
   }
print OUT ">$strain\n$seq\n";
close OUT;
}
