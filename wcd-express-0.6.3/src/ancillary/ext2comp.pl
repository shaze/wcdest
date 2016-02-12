#!/usr/bin/perl

#    Scott Hazelhurst
#    Program to convert an extended cluster table into a compressed
#    table


    use strict;

my ($line,@nums, $n, @link, @cluster, $k, $max,$sep);

$max = 0;
while ($line = <STDIN>) {
    if ($line =~ m/Index\s+Cluster/) { next;}
    if ($line =~ m/^[ \t]*$/) { next;}
    $n = ($line =~ /\s*([0-9]+)\s+([-0-9]+)\s+([-0-9]+)\s+([-0-9]+)\s+([-0-9]+)/);
    unless ($n) {die "*Error in input at line $.: $line\n"};

    $n = $1;
    $cluster[$n]=$2;
    if ($cluster[$n]==-1) { $cluster[$n]=$n; }
    $link[$n]=$3;
    $max++;
}

for ($n=0; $n<=$max; $n++) {
    if ($n == $cluster[$n]) {
	$sep = "";
	for($k=$n; $k>=0; $k=$link[$k]) {
	    print $sep.$k;
	    $sep = " ";
	}
	print ".\n";
    }
}

