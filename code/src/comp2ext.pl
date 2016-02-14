#!/usr/bin/perl

# Scott Hazelhurst
# Program to convert a compressed cluster table into a 
# an extended one

    use strict;

my ($line,@nums, $n, @link, @cluster, $max, $root, $prev);

$max = 0;
while ($line = <STDIN>) {
    $line =~ s/[\n\.]//g;
    @nums = split(" ",$line);
    $root = $nums[0];
    $prev = $root;
    foreach $n (@nums) {
	if ($n > $max) { $max = $n; }
	$cluster[$n]=$root;
	$link[$prev] = $n;
	$prev = $n;
    }
    $link[$prev] = -1;
}

for ($n=0; $n<=$max; $n++) {
    print "$n $cluster[$n] $link[$n] 0 -1\n";
}

