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
if($marker =~/CoV_(\d+)/){$base = $1;}
if($seq =~ /\[([A-Z])\/([A-Z])\]/){$lookup{$base}{$1}++; $lookup{$base}{$2}++; print "$base $1 $2\n";}
}


$meta = $ARGV[1];

open(META, $meta);
#sequence_name,country,adm1,sample_date,epi_week,lineage,lineage_support

$head = <META>;
while(<META>)
{
($name,$country,$adm,$date,$week,@other) = split(/\,/, $_);
$name2week{$name} = $week;
}
close META;

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
$match = substr($line,$pos,1);
if($lookup{$pos}{$match} >0){$final{$pos}{$week}{$match}++; $weeks{$week}++; $overall_total{$pos}{$match}++;}
}

}
}



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
