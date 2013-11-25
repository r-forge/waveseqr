#!usr/bin/perl -w

# This script reads a graph file and plugs in zeros
# for windows without any reads. Outputs padded graph
# file to STDOUT

# Apratim Mitra 2011-2012

my $infile = "";
my $winsize = 200;

foreach $i (0..scalar(@ARGV)-1){
    if($ARGV[$i] eq "-i"){
	$infile = $ARGV[$i+1];
    }
    elsif($ARGV[$i] eq "-w"){
	$winsize = $ARGV[$i+1];
    }
}

if(scalar(@ARGV) == 0 || $infile eq ""){
    die "Program Usage:\tperl pad_graph.pl -i infile -w winsize > outfile\n";
}

# check input file format
if($infile !~ /\.graph/){
	die "Input file has to be graph format!\n";
}

open(IN,$infile) or die "Could not open $infile\n";

$current = 0;
$currchr = "";

while(<IN>){
    $line = $_;
    $line =~ s/\s+$//;

    if($line =~ /chr(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/){
	($chr, $start, $end, $reads) = ($1, $2, $3, $4);

	if($chr eq $currchr){
	    if($start > $current){
		while(($current + $winsize - 1) < $start){
		    print "chr",$chr,"\t",$current,"\t",$current+$winsize-1,"\t0\n";
		    $current += $winsize;
		}
		print "chr",$chr,"\t",$start,"\t",$end,"\t",$reads,"\n";
		$current = $end + 1;
	    }
	    else{
		print "chr",$chr,"\t",$start,"\t",$end,"\t",$reads,"\n";
		$current = $end + 1;
	    }
	}
	else{
	    $currchr = $chr;
	    $current = 0;

	    if($start > $current){
		while(($current + $winsize - 1) < $start){
		    print "chr",$chr,"\t",$current,"\t",$current+$winsize-1,"\t0\n";
		    $current += $winsize;
		}
		print "chr",$chr,"\t",$start,"\t",$end,"\t",$reads,"\n";
		$current = $end + 1;
	    }
	    else{
		print "chr",$chr,"\t",$start,"\t",$end,"\t",$reads,"\n";
		$current = $end + 1;
	    }
	}
    }
}
close IN;