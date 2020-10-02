#!/usr/bin/perl


#Important note: all base poistions used in the files generated in this pipiline are in zero indexed format, i.e. base 0 is the first base in the alignment.
#Bear this in mind when comparing poistions in these files to externally refererenced ones, which are likley to begin with the first base as position 1.


$file = "$ARGV[0]";
$meta = "$ARGV[1]";
chomp $file;
chomp $meta;

#Set this variable to be the minimum minor allele frequency at a position considered real.  E.g. 0.01 means min of 1/100 reads
$min_maf = 0.01;

#Need to have an absoulte minimum as well, otherwise sngletons will be viewed as potential SNPs in small datasets.
$abs_min = 4;

#Set this variable to the minumum proportion of bases at a position with a CATG call (as opposed to ambiguous or missing).
$min_call_rate = 0.9;

#Set the weeks of data to be considered for analysis
$startweek = 0;
$endweek = 1000;
$maxmarkers = 24;



`./pre-process_covid_sequences.pl $file $min_maf $abs_min $min_call_rate $meta $startweek $endweek`;
`./process_covid_sequences.pl $file`;
`./convert_base_space_to_genotypes.pl $file`;
`./select_minimal_markers_v2.pl $file $maxmarkers`;
`./check_output.pl $file $meta $startweek $endweek`;
`./dereplicate_minimal_marker_lineages.pl $file`;
`./prepare_primer_submission.pl $file`;
`./allele_fre_by_date.pl primer_designs_for_cog.tdt $meta $file`;


