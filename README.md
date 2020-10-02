# SARSmarkers

This is a pipeline written in PERL to produce a minimal set of SNP markers capable of discriminating circulating lineages of SARS-CoV-2. It could equally be applied to a sequence alignent of any other organism but has only been tested with ~40,000 accessions of SARS-CoV-2 (~30 kilobases) and may not scale well to larger datasets.

To run:
1) Download all of the PERL scripts in this repositary along with the example alignment and metadata files into a single directory.
2) For sample data, download a copy of the May COG sequence alignment from: https://cog-uk.s3.climb.ac.uk/2020-05-08/cog_2020-05-08_alignment.fasta
and metadata https://cog-uk.s3.climb.ac.uk/2020-05-08/cog_2020-05-08_metadata.csv

3) Edit the pipeline.pl file to specify the outbreak week range you want to design primers for. The default is 0 -100. be aware that restricting the week range to a small interval may result in insufficient data to get a good primer design. You can also change the minimum minor allele frequency (default 0.01) minumum call-rate (to avoid positions with lots of N's - default 0.9) and absolute minimum allele call rate - default 4, in case you are running a small dataset and don't want SNPs with less than 4 reads to count as real. 

4) Run the pipeline with ./pipeline.pl cog_2020-05-08_alignment.fasta cpg_2020-05-08_metadata.csv






