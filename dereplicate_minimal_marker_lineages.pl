#!/usr/bin/perl



$markers = $ARGV[0];
$markers =~ s/\.fasta/_minimal_markers.txt/;
$markers = "annotated_$markers";

chomp $markers;
open(IN, "$markers");
open(OUT, ">dereplicated_$markers");

$head = <IN>;
chomp $head;
($id, @head) = split(/\t/, $head);
$lineage = <IN>;
chomp $lineage;
($id, @lineage) = split(/\t/, $lineage);


while(<IN>)
{
chomp;
($marker, @data) = split(/\t/, $_);
$i = 0;
push @markers, $marker;
 foreach $cell(@data){$transpose[$i].=$cell; $i++;}
}


$n = @markers;
print "$n marker rows\n";

print OUT "ID\tLineage";
foreach $marker(@markers){print OUT "\t$marker";}
print OUT "\n";

$i = 0;
foreach $seq(@transpose)
{
$rseq = $seq; $rseq =~ s/x/\\D/g;
$sample = $head[$i];
$lineage = $lineage[$i]; 
$replicate{$seq}{$lineage}++;
$x = 0; foreach $seq2(@transpose){$lineage2 = $lineage[$x]; $x++; $rseq2 = $seq2; $rseq2 =~ s/x/\\D/g; if ($seq ne $seq2 && ($seq =~ /$rseq2/ || $rseq =~ /$seq2/)){$replicate{$seq2}{$lineage2}++;}}
 if($replicate{$seq}{$lineage} ==1)
 {
 print OUT "$sample\t$lineage";
 @seq = split(//, $seq);
 $n = @seq; print "now $n characters\n";
 $seq2 = join("\t", @seq);
 print OUT "\t$seq2\n";
 }
else{print "Removed $sample $lineage $seq as a replicate\n";}
$i++;
}
