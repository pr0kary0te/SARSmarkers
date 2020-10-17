#!/usr/bin/perl

open(IN, "genotyping_calls_PHE.tdt");
$min_call_rate = 0.5;

$head = <IN>;
chomp $head;
($id, @header) = split(/\t/, $head);


while($line = <IN>)
{
$rows++;
chomp $line;
($id, @data) = split(/\t/, $line);
$data = join("\t", @data);
$id2data{$id} = $data;
$i = 0;
$cols = @data;
foreach $cell(@data)
  {
  $lookup{$id}{$i} = $cell;
  if($cell =~ /[CATG]/){$good_id{$id}++; $good_marker{$i}++;} 
  if($cell =~/([CATG])\:([CATG])/){if($1 eq $2){$hom_marker{$i}++;} elsif($1 ne $2){$het_marker{$i}++;}}
  $i++; 
 }
if(($good_id{$id} / @data) < $min_call_rate){delete $lookup{$id}; delete $good_id{$id}; 
}
}

print "Call rates\n";
foreach $id(sort keys %good_id)
{
$n = $good_id{$id};
push @{$sort_by_good{$n}}, $id;
}

print "\nMarker\tnum calls\thom calls\thet calls\n";
foreach $i(sort keys %good_marker)
{
$id = $header[$i];
print "$id\t$good_marker{$i}\t$hom_marker{$i}\t$het_marker{$i}\n";

}
print "\n\n\n";

foreach $id(sort keys %good_id)
 {
 foreach $id2(sort keys %good_id)
  {
  if($id ne $id2)
   {
   $i = 0;
   while($i < $rows)
    {
    $call1 = $lookup{$id}{$i};
    $call2 = $lookup{$id2}{$i};
    if($call1 =~ /[CATG]:[CATG]/ && $call2 =~ /[CATG]:[CATG]/ && $call1 ne $call2){$final{$id}{$id2}++; $final{$id2}{$id}++;}
    $i++;
    }
   }
  }
 }

print "Sample difference matrix\n";
foreach $id (sort keys %good_id){print "\t$id";}
print "\n";
foreach $id(sort keys %good_id)
{
print "$id";
foreach $id2(sort keys %good_id)
{
if($id ne $id2)
 {
 $n = $final{$id}{$id2};
 if($final{$id}{$id2} !~ /\d/){$n = 0;}
 print "\t$n";
 }

else{print "\t-";}
}
print "\n";
}

$i = 0;

foreach $n(sort {$b<=>$a} keys %sort_by_good)
{
$i = 0;
$ref = $sort_by_good{$n};
@ary = @$ref;

foreach $id(@ary)
 {
 push @bigary, $id;
 }
}


foreach $id(@bigary)
{
$row++;
$k++; 
$col = $row;
 foreach $id2(@bigary)
 {
 $col++;
 if($col >= $row)
 {
 if($final{$id}{$id2} ==0 || $final{$id2}{$id} ==0)
  {
  #Has $id2 been assigned to a group previously as group member?  If yes, ignore it, if not, add it to the members hash and add it to $id's group
  if($second{$id2} <1 && $second{$id} <1){$groups{$id}{$id2}++;  if($id ne $id2){$second{$id2}++;} } 
  }
}
 }

}





print "\n\n\nGroup analysis\n\n$head\tgroup\tn\tmembers\n";
foreach $id(keys %groups)
{
$groupnum++;
#print "$id";
print "$id\t$id2data{$id}\t";
$ref = $groups{$id}; 
%hash = %$ref;
$j = keys %hash; 
foreach $id2(sort keys %hash){print "${id2},"; push @lottery, $groupnum;}
print "\n";
}


$i = 0;
while($i <1000)
{
$i++;
shuffle(\@lottery);
if($lottery[0] == $lottery[1]){$same++;}
}
print "\n\n$same random pairs of 1000 samples had the same genotype\n";


sub shuffle
{
my $array = shift;
my $i;
for ($i = @$array; --$i;) {
my $j = int(rand($i+1));
next if $i == $j;
@$array[$i,$j] = @$array[$j,$i];
}
}
