#!usr/bin/perl -w

# This script converts a bed file to a bedgraph file.
# BED must be sortted before use with this  script.

# Apratim Mitra 2011-2012

# default parameters
my $winsize = 200;
my $fragsize = 150;
my $infile = "";

foreach $i (0..scalar(@ARGV)-1){
    if($ARGV[$i] eq "-i"){
	$infile = $ARGV[$i+1];
    }
    elsif($ARGV[$i] eq "-w"){
	$winsize = $ARGV[$i+1];
    }
    elsif($ARGV[$i] eq "-f"){
	$fragsize = $ARGV[$i+1];
    }
}

if($infile eq ""){
    die "\nProgram usage:\n\tperl bed2graph.pl -i bedfile -w winsize(def=200) -f fragsize(def=150)>\n\n";
}
elsif($infile !~ /\.bed$/){
	die "Input file has to be in BED format!\n";
}
open(IN,$infile) or die "Could not open $infile\n";

# chromosome flag
$currchr = "";

while(<IN>){
    $line = $_;
    $line =~ s/\s+$//;

	# bed file line format: <chr> <start> <end> <text> <length> <strand>
    if($line =~ /(chr\S+)\s+(\d+)\s+(\d+)\s+.+(\+|\-)/){
	($chr,$start,$end,$strand) = ($1,$2,$3,$4);	
    }    

    # Check if on the same chromosome
    if($currchr ne $chr){
	if($currchr ne ""){
	    # pass hash to subroutine	    
	    @chrkeys = keys %chrhash;	    
	    &countReads($currchr,$winsize,\@chrkeys);	    
	}
	$currchr = $chr;
	%chrhash = ();
    }

    # Extend the strand both ways
    if($strand eq "+"){
	$mid = $start + $fragsize/2;
    }
    else{
	$mid = $end - $fragsize/2;
    }

    # read into hash indexed by midpt
    if(exists $chrhash{$mid}){
	$chrhash{$mid}++;
    }
    else{
	$chrhash{$mid} = 1;
    }
}
# for the last chromosome
@chrkeys = keys %chrhash;
&countReads($currchr, $winsize, \@chrkeys);

# subroutine to count the reads using quick sort
sub countReads{
    my ($currchr,$winsize,$arrayptr) = ($_[0],$_[1],$_[2]);
    my ($currcount, $currstart, $currend);    

    &quicksort($arrayptr,0,scalar(@$arrayptr)-1);

    $sortkeys = $arrayptr;    

    $currstart = 0;
    $currend = $currstart + $winsize - 1;
    $currcount = 0;    

    $arraylen = scalar(@$arrayptr);
    foreach $i (0..$arraylen-1){
	if($$sortkeys[$i] >= $currstart && $$sortkeys[$i] <= $currend){
	    $currcount += $chrhash{$$sortkeys[$i]};
	}
	else{
	    if($currcount > 0){
		print $currchr,"\t",$currstart,"\t",$currend,"\t",$currcount,"\n";
	    }

	    while($$sortkeys[$i] > $currend){
		$currstart += $winsize;
		$currend += $winsize;		
	    }
	    $currcount = $chrhash{$$sortkeys[$i]};
	}
    }
    return;
}

# partition subroutine for quicksort
sub partition{
    my ($arrayref, $left, $right, $pivotind) = ($_[0],$_[1],$_[2],$_[3]);    
    my $storeind = $left;
    my $pivot = $$arrayref[$pivotind];

    #swap pivot and right element
    ($$arrayref[$pivotind],$$arrayref[$right]) = ($$arrayref[$right],$$arrayref[$pivotind]);
    for $i ($left..$right-1){
	if($$arrayref[$i] <= $pivot){
	    ($$arrayref[$i],$$arrayref[$storeind]) = ($$arrayref[$storeind],$$arrayref[$i]);
	    $storeind++;
	}
    }
    ($$arrayref[$storeind],$$arrayref[$right]) = ($$arrayref[$right],$$arrayref[$storeind]);    
    return $storeind;
}

# quicksort subroutine
sub quicksort{
    my ($arrayref, $left, $right) = ($_[0],$_[1],$_[2]);
    my $pivotind = 0;

    if($right > $left){
	#select pivot
	$pivotind = ($left+$right)/2;
	$pivotnewIndex = &partition($arrayref, $left, $right, $pivotind);
	&quicksort($arrayref, $left, $pivotnewIndex - 1);
	&quicksort($arrayref,$pivotnewIndex + 1, $right);
    }    
    return;
}
