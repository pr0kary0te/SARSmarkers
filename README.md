# SARSmarkers

This is a pipeline written in PERL to produce a minimal set of SNP markers capable of discriminating circulating lineages of SARS-CoV-2. It could equally be applied to a sequence alignent of any other organism but has only been tested with ~40,000 accessions of SARS-CoV-2 (~30 kilobases) and may not scale well to larger datasets.

To run:
1) Download all of the PERL scripts in this repositary along with the example alignment and metadata files into a single directory.
2) For sample data, download a copy of the May COG sequence alignment from: https://cog-uk.s3.climb.ac.uk/2020-05-08/cog_2020-05-08_alignment.fasta
and metadata https://cog-uk.s3.climb.ac.uk/2020-05-08/cog_2020-05-08_metadata.csv

3) Edit the pipeline.pl file to specify the outbreak week range you want to design primers for. The default is 0 -1000. be aware that restricting the week range to a small interval may result in insufficient data to get a good primer design. You can also change the minimum minor allele frequency (default 0.01) minumum call-rate (to avoid positions with lots of N's - default 0.9) and absolute minimum allele call rate - default 4, in case you are running a small dataset and don't want SNPs with less than 4 reads to count as real. Finally you can change the value of $maxmarkers - the maximum marker panel size the software will create.  Once this number of markers is reached, no more will be added, even if they add further lineage discrimination.  We added this feature as we found that as more mutations accumulate in the COG data, the pipeline could keep adding SNPs even when thay add the ability to discriminate only extremely rare lineages, and thus have add little value.

In our analysis of the 2020-09-03 COG release, we found that 24 markers was sifficient to discriminate >95% of randomly selected sample pairs, with the tradeoff of being a small marker panel. Note that the pipeline adds the markers with the highest discriminatory power first, so truncating at 24 will always select the best 24 markers.


4) Run the pipeline with ./pipeline.pl cog_2020-05-08_alignment.fasta cog_2020-05-08_metadata.csv (or see point 6 below).

5) If using a different dataset, edit the filenames and be sure to use the matching aligment and metadata files.

6) Experimental option: use ./wget_pipeline.pl : this should download the latest COG data via WGET (if that's on your system) and run the pipeline on these automatically. THis was tested with the 2020-09-03 dataset and worked fine, but future changes to the COG website could break this so use with caution.

Additional analysis
The script qc_genotype_data.pl was used to group samples into groups of identical genotypes after analyisis with PACE markers for publication.  This used the data in the file genotyping_calls_PHE.tdt.  Just run ./qc_genotype_data.pl 






