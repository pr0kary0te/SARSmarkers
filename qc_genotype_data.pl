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
if(($good_id{$id} / @data) < $min_call_rate){delete $lookup{$id}; delete $good_id{$id};}
}


print "ID\tnum calls\n";
foreach $id(sort keys %good_id)
{
print "$id\t$good_id{$id}\n";
$n = $good_id{$id};

push @{$sort_by_good{$n}}, $id;
}

print "Marker\tnum calls\thom calls\t$het calls\n";
foreach $i(sort keys %good_marker)
{
$id = $header[$i];
print "$id\t$good_marker{$i}\t$hom_marker{$i}\t$het_marker{$i}\n";

}


foreach $id(sort keys %lookup)
{
foreach $id2(sort keys %lookup)
{
if($id ne $id2)
{
$i = 0;
while($i < $rows)
{
$call1 = $lookup{$id}{$i};
$call2 = $lookup{$id2}{$i};
if($call1 =~ /[CATG]:[CATG]/ && $call2 =~ /[CATG]:[CATG]/ && $call1 ne $call2){$final{$id}{$id2}++;}
$i++;
}
}
}
}
print "\n\n\n";

foreach $id (sort keys %lookup){print "\t$id";}
print "\n";
foreach $id (sort keys %lookup)
{
print "$id";

foreach $id2 (sort keys %lookup)
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
$ref = $sort_by_good{$n};
@ary = @$ref;
print "$n";
foreach $id(@ary)
{
print "\t$id";
$i++; if($i ==1){$groups{$id}{$id}++; $used{$id}++;}
foreach $id2(keys %lookup){if($used{$id2} ==0 && $final{$id}{$id2} ==0){$groups{$id}{$id2}++; $used{$id2}++;}}
if($used{$id} ==0)
 {
 $groups{$id}{$id}++; 
 foreach $id2(keys %lookup){if($used{$id2} ==0 && $final{$id}{$id2} ==0){$groups{$id}{$id2}++; $used{$id2}++;}}

 
 }

}
print "\n";
}

print "Groups\n$head\tgroup\tn\tmembers\n";
foreach $id(keys %groups)
{
$groupnum++;
#print "$id";
print "$id\t$id2data{$id}\t";
$ref = $groups{$id}; %hash = %$ref;
$j = keys %hash; 
#print "$groupnum\t$j\t";
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
print "$same random pairs of 1000 samples had the same genotype\n";


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
