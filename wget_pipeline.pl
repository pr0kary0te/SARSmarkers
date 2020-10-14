#!/usr/bin/perl

`rm index.htm*`;
`rm cog*alignment.fast*`;
`rm cog*metadata.cs*`;
`rm subsequence*.txt`;

$wget = `wget https://www.cogconsortium.uk/data/`;
print "$wget";
$alignment_url = `grep alignment.fasta index.html`;

#example url format:
#https://cog-uk.s3.climb.ac.uk/2020-09-03/cog_2020-09-03_alignment.fasta

if($alignment_url =~ /(https:\/\/cog-uk.s3.climb.ac.uk\/.*\/cog_.*_alignment.fasta)/){$alignment_url = $1; }
else{die "Can't find a valid alignment file url in the web page source https://www.cogconsortium.uk/data/ 
- check the parsing REGEX in line 10 of this script against the page source\n";}

$metadata_url = `grep _metadata.csv index.html`;

#example url format:
#https://cog-uk.s3.climb.ac.uk/2020-09-03/cog_2020-09-03_metadata.csv
if($metadata_url =~ /(https:\/\/cog-uk.s3.climb.ac.uk\/.*\/cog_.*_metadata.csv)/){$metadata_url = $1; }
else{die "Can't find a valid metadata file url in the web page source https://www.cogconsortium.uk/data/ 
- check the parsing REGEX in line 20 of this script against the page source\n";}


if($alignment_url =~ /(cog_.*_alignment.fasta)/){$file = $1;} else{die "alignment file name doesn't match cog_.*_alignment.fasta"; }
if($metadata_url =~ /(cog_.*_metadata.csv)/){$meta = $1;} else{die "metadata file name doesn't match cog_.*_metadata.csv"; }

print "Getting $alignment_url with WGET\n";
`wget $alignment_url`;

print "Getting $metadata_url with WGET\n";

`wget $metadata_url`;

if(-e $file){} else{die "could not locate expected alignment $file - check for complete download\n";}
if(-e $meta){} else{die "could not locate expected metadata $meta - check for complete download\n";}



chomp $file;
chomp $meta;

#Set this variable to be the minimum minor allele frequency at a position considered real.  E.g. 0.01 means min of 1/100 reads
$min_maf = 0.01;

#Need to have an absoulte minimum as well, otherwise sngletons will be viewed as potential SNPs in small datasets.
$abs_min = 4;

#Set this variable to the minumum proportion of bases at a position with a CATG call (as opposed to ambiguous or missing).
$min_call_rate = 0.9;

#Set the weeks of data to be considered for analysis
$startweek = 28;
$endweek = 1000;




`./pre-process_covid_sequences.pl $file $min_maf $abs_min $min_call_rate $meta $startweek $endweek`;
`./process_covid_sequences.pl $file`;
`./convert_base_space_to_genotypes.pl $file`;
`./select_minimal_markers_v2.pl $file $min_maf $min_call_rate`;
`./check_output.pl $file $meta $startweek $endweek`;
`./dereplicate_minimal_marker_lineages.pl $file`;
`./prepare_primer_submission.pl $file`;
`./allele_fre_by_date.pl primer_designs_for_cog.tdt $meta $file`;


#Cleanup 
`rm subsequence*_alignment_minimal_markers.txt`;
