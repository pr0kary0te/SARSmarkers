#!/usr/bin/perl


$file = $ARGV[0];
chomp $file;
open(LOG, ">>$file.pipeline.log");

$date = `date`;
chomp $date;

$meta = $ARGV[4]; 
if($meta !~ /\D/){$meta = "cog_2020-09-03_metadata.csv";}
$startweek = $ARGV[5];
$endweek = $ARGV[6];


open(META, $meta);
#Was originally: sequence_name,country,adm1,sample_date,epi_week,lineage,lineage_support
#Jan 2021 COG now use: sequence_name,country,adm1,TrueFalse,sample_date,epi_week,lineage,lineage_support
$head = <META>;
while(<META>)
{
#Jan 2021 - changed REGEX as COG have introduced a new column in the metadata!
#($name,$country,$adm,$date,$week,@other) = split(/\,/, $_);
($name,$country,$adm,$truefalse,$date,$week,@other) = split(/\,/, $_);
$name2week{$name} = $week;
$weeks{$week}++;
}
close META;
  print LOG "$date\nMetadata $meta\nWeek\tSample count\n";
foreach $week (sort {$a<=>$b} keys %weeks){print LOG "$week\t$weeks{$week}\n";}

$keys = keys %name2week;
print LOG "Loaded week data for $keys samples\n";

#$file = "cog_2020-05-08_alignment.fasta";

print LOG "\n\npre-process_covid_sequences.pl seqfile $ARGV[0] minMAF $ARGV[1] ABS-min $ARGV[2] min call rate $ARGV[3] metadata $ARGV[4] startweek $ARGV[5] endweek $ARGV[6]";

#Set this variable to be the minimum minor allele frequency at a position considered real.  E.g. 0.01 means min of 1/100 reads
$min_maf = $ARGV[1]; if($min_maf ==0){$min_maf = 0.01;}

#Need to have an absoulte minimum as well, otherwise sngletons will be viewed as potential SNPs in small datasets.
$abs_min = $ARGV[2]; if($abs_min ==0){$abs_min = 4;}

#Set this variable to the minumum proportion of bases at a position with a CATG call (as opposed to ambiguous or missing).
$min_call_rate = $ARGV[3]; if($min_call_rate ==0){$min_call_rate = 0.9;}





open(OUT, ">${file}_covid_variants.tdt");

chomp $file;
open(IN, "$file");

while(<IN>)
  {
  chomp;
  if(/>(.*)/){$id = $1; $name = $id; $id =~ s/[^0-9A-Z-a-z]/_/g; $week = $name2week{$name}; 
  
#Changed comparison operators to >= and <= (were < and > in error)
if ($week >= $startweek && $week <= $endweek){push @ids, $id;}}

  else
    {
    $week = $name2week{$name};
    if($week >= $startweek && $week <= $endweek)
      {
      $goodweekcriteria++;
      $pos = 0;
      @sequence = split(//, $_); $n = @sequence;
 
     #Split the sequences into column positions and store bases at that position in an array, basically rotating the sequence alignment by 90 degrees
     #and viewing a row of base calls at a poition now.: 
      foreach $base(@sequence){$genotype[$pos].=$base; $pos++;}
      }
    }
  }
close IN;
print LOG "\n$goodweekcriteria sequences fell between weeks $startweek and $endweek\n";


$types = @ids;
$min = $types * $min_maf;
if($min < $abs_min){$min = $abs_min;}

#Prepare the header for output sequence IDs
foreach $id(@ids)
{
print OUT "\t$id"; 
}
print OUT "\n"; 



print LOG "\nLoaded $types sequences\n";

$i = 0;
foreach $genotype(@genotype)
   {
   $i++;
   $bad=0;
   %bases = ();
   %orderbycount = ();
   @counts = ();
   
   @bases = split(//, $genotype);
   foreach $base(@bases){$base =~ s/[^CATG]/N/; if($base =~ /[CATG]/){$bases{$base}++;} else{$bad++;}}
   
   %orderbycount = ();
   foreach $base(keys %bases)
     {
     $n = $bases{$base};
     push @{$orderbycount{$n}}, $base;
     }
   $round = 0; $displaytext = "";
   foreach $n (sort {$b<=>$a} keys %orderbycount)
     {
     $ref = $orderbycount{$n};
     @ary = @$ref;
     foreach $base(@ary)
       {
       $round++; 
       if($round ==1){$displaytext = "${base}$n";}
       if($round == 2 && $n > $min && (1-($bad/$types)) >$min_call_rate )
       {
       #Only the two most abundant bases will be printed out as we're interested in biallelic SNP calling. Third and fourth most common bases (if present) are ignored.
       $displaytext .= "${base}$n";
       print OUT "Base $i $displaytext\t$genotype\n";
       }
     }
   }
}

