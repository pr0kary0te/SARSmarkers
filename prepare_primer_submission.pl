#!/usr/bin/perl

$min_maf = 0.005;

#Example input file types
#subsequence_1058_cog_2020-05-08_alignment_minimal_markers-mincallrate-0.5-minmaf-0.001.txt

$study = $ARGV[0];
chomp $study;
$study =~ s/([^_]+).*/$1/;

open(OUT, ">primer_designs_for_$study.tdt");

@files = `ls subsequence_*_${study}*`;
foreach $file(@files)
{
chomp $file;
print "$file\n";
if($file =~ /subsequence_([0-9]+)_${study}/)
{
$pos = $1; 
%first_bases =(); %last_bases =();
open(IN, "$file");
while(<IN>)
{
chomp;
($seq, $count) = split(/\t/, $_);
($first, $ambiguity, $last)  = split(/[\[\]]/, $seq);

$i = 0; while($first =~ /([A-Z-])/g){$first_bases{$i}{$1}+=$count; $i++;}
$i = 0; while($last =~ /([A-Z-])/g){$last_bases{$i}{$1}+=$count; $i++;}

}


$i = 0;
$newfirstseq ="";
foreach $i(sort {$a<=>$b} keys %first_bases){$n = 0; $top = 0; $next = 0; $ref = $first_bases{$i}; %hash = %$ref; 
foreach $base(keys %hash){$n = $hash{$base}; 
if($n >$top){$top = $n; $topbase = $base;} elsif ($n >$next){$next = $n; $nextbase = $base;}}
if($nextbase !~ /[N\-]/ && $top >0 && $next > 0 && $next/$top > $min_maf)
{ $pair = "$topbase$nextbase"; $base = ambiguity($pair); } 
else{$base = $topbase;}
$newfirstseq.=$base;
}

$i = 0;
$newlastseq ="";
foreach $i(sort {$a<=>$b} keys %last_bases){$n = 0; $top = 0; $next = 0; $ref = $last_bases{$i}; %hash = %$ref;
foreach $base(keys %hash){$n = $hash{$base};
if($n >$top){$top = $n; $topbase = $base;} elsif ($n >$next){$next = $n; $nextbase = $base;}}
if($nextbase !~ /[N\-]/ && $top >0 && $next > 0 && $next/$top > $min_maf)
{ $pair = "$topbase$nextbase"; $base = ambiguity($pair
); }
else{$base = $topbase;}
$newlastseq.=$base;
}

print OUT "BRIS_CoV_$pos\t$newfirstseq\[$ambiguity\]$newlastseq\n";




close IN;
}
else{die "File format is incorrect: expecting something like \"subsequence_2890_cog_\"\n";}

}

sub ambiguity
{
my $input = shift;

@final = split(//, $input);;
@final = sort(@final);
$final = join("", @final);

if($final eq "AG"){$final = "R";}
elsif($final =~ /R/){$final = "R";}

elsif($final eq "CT"){$final = "Y";}
elsif($final =~ /Y/){$final = "Y";}

elsif($final eq "GT"){$final = "K";}
elsif($final =~ /K/){$final = "K";}


elsif($final eq "AC"){$final = "M";}
elsif($final =~ /M/){$final = "M";}

elsif($final eq "CG"){$final = "S";}
elsif($final =~ /S/){$final = "S";}

elsif($final eq "AT"){$final = "W";}
elsif($final =~ /W/){$final = "W";}

return $final;
}


