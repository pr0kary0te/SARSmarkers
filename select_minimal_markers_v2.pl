#!/usr/bin/perl



$maxmarkers =50000;
#it was set at 50000

#This is currently set to more than the number of markers in the input file (~ 35000) but if it is set lower then markers are prioritised.
#E.g. if set to 5000, then the top 5000 markers by MAF score will be used and the rest ignored.


# ./select_minimal_markers.pl $ARGV[0] $ARGV[1] $ARGV[2]
# ./select_minimal_markers.pl data.txt 0.001 0.3

$min_maf = $ARGV[1];
#it was set at 0.001
#MAF is minor allele frequency.  This is set to a low level to include as many markers as possible but exclude rare error calls.
#It needs optimising.
if($min_maf !~ /\d/){$min_maf = 0.001;}


$min_call_rate = $ARGV[2];
chomp $min_call_rate;
#it was set at 0.3
#Ignore markers with less than this proportion of good (0, 1 of 2 ) calls.  Probably needs optimising.

if($min_call_rate !~ /\d/){$min_call_rate = 0.5;}

#Get the input file name from the command line - first argument:  e.g. ./select_minimal_markers.pl All_marker_data.txt
$infile = $ARGV[0];


chomp $infile;


$infile = "${infile}_covid_variants_as_genotypes.tdt";



#Initiate some hashes.
my %matrix = ();
my %testmatrix =();


#Open the input data file handle - the input file having been specified on the command line
open(IN, "$infile" ||die "Can't opemn $infile\n");
$oldfile = $infile;
$infile =~ s/\..*/_minimal_markers.txt/;
open(OUT, ">$infile");


#Parse the file header before going through all the data lines
#The file format for the header is "Marker name -> Variety1 name-> Variety2 name-> Variety3 name...."
$head = <IN>;
($id, @header) = split(/\t/, $head);
$head = join("\t", @header);
$hlen = @header;
print OUT "Power\t$marker\t$head";

#Start reading the data here
#The file format for the data is "Marker name -> Variety1 score  -> Variety2 score-> Variety3 score...."
while(<IN>)
{
chomp;
($id, @data) = split(/\t/, $_); 
%alleles = ();
foreach $cell(@data)
  {
  #Replace any cells which aren't 0, 1 or 2 with "x"
  if($cell !~ /^[012]/) {$cell = "x";} 

  #Make a hash list of alleles observed in this current row  which aren't bad "x" calls
  if($cell !~ /x/){$alleles{$cell}++;}
  }

$thislen = @data;
#Check that the header and each data row have the same number of cells, in case of file corruption. Die if not.
if($hlen != $thislen){print "$id has length of $thislen which doesn't match header ($hlen)\n"; die;}

$data = join("", @data);
$pattern2id{$data} = $id;
}


$n = keys %pattern2id;

print LOG "Loaded marker data for $n distinct patterns\n";
close IN;




#set up n x n matrix for data
my $l = @data;

$x = 0; $y = 0;

while($x < $l){$y = 0; while($y < $l){$matrix{$x}{$y} = 0; 
#print "Matrix $l build $x $y\n"; 
$y++;} $x++;
}

print "Built initial $l x $l scoring matrix for lineages\n";


$currentscore = 1;
while($currentscore > 0)
{
%testmatrix =();
$bestscore = 0;
$first = 0;
$date = `date`;
$iteration++; print "\n$date Iteration $iteration ";

foreach $pattern1(sort keys %pattern2id)
   {
   $first++;
   $score = 0;
   $id = $pattern2id{$pattern1};
   %testmatrix =();
   @pattern1 = split(//, $pattern1);
   #print "$id\t$pattern1\n";
   $i = 0; $j = 0;
   $len = @pattern1;
   while($i < $len)
      {
      $ichar = $pattern1[$i];
      $j = $i +1;
      while($j <$len)
         {
         $jchar = $pattern1[$j];
         if($matrix{$i}{$j} ==0 && $ichar ne "x" && $jchar ne "x" && $ichar ne $jchar) 
             {
             $testmatrix{$i}{$j}++; $score++;  
             }  
         $j++;
         }
      $i++;
      }
 

  if($score >$bestscore)
      {
      $bestscore = $score; 
      %bestmatrix = %testmatrix; 
      $bestpattern = $pattern1;     
      }
   }

       $maf = $pattern2maf{$bestpattern}; $id = $pattern2id{$bestpattern}; 
      if($bestscore >0)
         {
         print "$bestscore\t$id\t$maf\t$bestpattern\n"; 
         $print = $bestpattern;
         @print = split(//, $print);
         $print = join("\t", @print);
         print OUT "$bestscore\t$id\t$print\n"; 
         }
        $currentscore = $bestscore;
       $i = 0; $j = 0;
         @pattern1 = split(//, $bestpattern);
       $len = @pattern1;
       $improved =0;
        
        while($i < $len)
         {
         $j = $i +1;
         while($j <$len)
            {
            $matrix{$i}{$j} += $bestmatrix{$i}{$j};  
            $improved += $bestmatrix{$i}{$j};
            $j++;
            }
         $i++;
         }
$bestscore =0;
%bestmatrix =();
%sortbyscore =();
}







