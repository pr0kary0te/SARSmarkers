#!/usr/bin/perl

my %lookup =();


#Usage:  ./ORFfinder2.pl input.fasta [min_length]
# Where min_length is the minimum length ORF that you want to consider an ORF

$seq = $ARGV[0];



$min = 10;
$rev = reverse($seq);
$rev =~ tr/CATG/GTAC/;


$i = 0; while($i <6){if($i ==4){$seq = $rev; }$frame = translate($seq); push @possibles, $frame;  $seq =~ s/^.//; $i++;}

$max = 0;
$i = 0;
foreach $possible(@possibles)
{
$i++; 
@data = split(/X/, $possible);
foreach $string(@data){$l = length($string); if($l >$max){$max = $l; $best = $string; $bestframe = $i;}}
}

if($max >$min)
{
$frame = $peptide2frame{$best};
print "$best\t$bestframe\n";
}

sub translate
{
$lookup{AAC}="N";
$lookup{AAG}="K";
$lookup{ACA}="T";
$lookup{ACC}="T";
$lookup{ACG}="T";
$lookup{AGA}="R";
$lookup{AGC}="S";
$lookup{AGG}="R";
$lookup{AGT}="S";
$lookup{ATA}="I";
$lookup{ATC}="I";
$lookup{CAC}="H";
$lookup{CAG}="Q";
$lookup{CCA}="P";
$lookup{CCC}="P";
$lookup{CCG}="P";
$lookup{CGA}="R";
$lookup{CGC}="R";
$lookup{CGG}="R";
$lookup{CTA}="L";
$lookup{CTC}="L";
$lookup{CTG}="L";
$lookup{GAC}="D";
$lookup{GAG}="E";
$lookup{GCA}="A";
$lookup{GCC}="A";
$lookup{GCG}="A";
$lookup{GGA}="G";
$lookup{GGC}="G";
$lookup{GGG}="G";
$lookup{GTA}="V";
$lookup{GTC}="V";
$lookup{GTG}="V";
$lookup{TAC}="Y";
$lookup{TCA}="S";
$lookup{TCC}="S";
$lookup{TCG}="S";
$lookup{TGC}="C";
$lookup{TTA}="L";
$lookup{TTC}="F";
$lookup{TTG}="L";
$lookup{AAA}="K";
$lookup{AAT}="N";
$lookup{ACT}="T";
$lookup{ATG}="M";
$lookup{ATT}="I";
$lookup{CAA}="Q";
$lookup{CAT}="H";
$lookup{CCT}="P";
$lookup{CGT}="R";
$lookup{CTT}="L";
$lookup{GAA}="E";
$lookup{GAT}="D";
$lookup{GCT}="A";
$lookup{GGT}="G";
$lookup{GTT}="V";
$lookup{TAT}="Y";
$lookup{TCT}="S";
$lookup{TGG}="W";
$lookup{TGT}="C";
$lookup{TTT}="F";


 $lookup{TAG}="X";
 $lookup{TAA}="X";
 $lookup{TGA}="X";

my $protein ="";
my $input = shift;
chomp $input;
@codons = split(/(\D{3})/, $input);
#print "\n@codons\n\n";
foreach $codon(@codons){$aa = $lookup{$codon}; $protein.= $aa;}
return($protein);
}
