#!usr/bin/perl -w

# This script removes duplicate reads and outputs to STDOUT

# Apratim Mitra 2011-2012

my $infile = "";
my $sort = 0;
my $system = "";

foreach $i (0..scalar(@ARGV)-1){
    if($ARGV[$i] eq "-i"){
	$infile = $ARGV[$i+1];
    }	
	elsif($ARGV[$i] eq "-s"){
	$sort = $ARGV[$i+1]; 
	}
}

if($infile eq ""){
    die "No files input!\n";
}

if($infile !~ /(.+)(\.bed)$/){
	die "File has to be in bed format!\n";
}

# Expected fields:
# chr  start  end  txt  len  +/-

# sort file by chromosome and start (if needed)
if($sort == 1){
	$cmd = "sort -k 1,1 -k 2,2n ".$infile."-o ".$infile;
	$system = `$cmd`;
}

open(IN,$infile) or die "Could not open $infile\n";

$plusstart = 0;
$plusend = 0;

$minusstart = 0;
$minusend = 0;

$total = 0;
$unique = 0;

while(<IN>){
	$line = $_;
	$line =~ s/\s+$//;
	
    @cat = split(/\t/,$line);	
	$total++;
	if($cat[5] eq "+"){		
		if($cat[1] == $plusstart){	
		next;
		}
		else{	
		$unique++;

		$plusstart = $cat[1];
		$plusend = $cat[2];	
		print $_;
		}
	}
	elsif($cat[5] eq "-"){		
		if($cat[1] == $minusstart){	
		next;
		}
		else{	
		$unique++;

		$minusstart = $cat[1];
		$minusend = $cat[2];	
		print $_;
		}
	}
}

print "\tTotal reads = $total\n";
print "\tUnique reads = $unique\n";
