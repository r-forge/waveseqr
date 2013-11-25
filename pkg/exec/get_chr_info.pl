#!usr/bin/perl

# This program counts the number of lines
# for each chromosome in an input file
# Also output totalreads in last line.
# (bed or bedgraph)

my $infile = "";
my $outfile = "";

foreach $i (0..scalar(@ARGV)-1){
    if($ARGV[$i] eq "-i"){
	$infile = $ARGV[$i+1];
    }
	elsif($ARGV[$i] eq "-o"){
	$outfile = $ARGV[$i+1];
    }
}

if($infile eq "" || $outfile eq ""){
    die "No files input!\nProgram Usage:\tperl get_chr_info.pl -i infile -o outfile\n";
}
elsif($infile =~ /\.bed$/){
    $type = 1;
}
elsif($infile =~ /\.graph$/){
    $type = 2;
}
else{
    die "File has to bed or bedgraph format!\n";
}

my $currchr = "";
my $currlines = 0;
my $totalreads = 0;
open(IN,$infile) or die "Could not open $infile\n";
open(OUT,">$outfile");
while(<IN>){
    my $line = $_;
    $line =~ s/\s+$//;

    my @col = split(/\t/,$line);
    # count number of lines for each chromosome
    if($currchr eq $col[0]){
	$currlines++;
    }
    else{
		if($currchr ne ""){
			print OUT $currchr,"\t",$currlines,"\n";
		}	
	$currchr = $col[0];
	$currlines = 1;
    }

    # count total reads
    if($type == 1){
	$totalreads++;
    }
    elsif($type == 2){
	$totalreads += $col[3];
    }
}
close IN;
#print last line
#if($currchr !~ /random/i){
	print OUT $currchr,"\t",$currlines."\n";
#}

#foreach $chr (sort keys %chrhash){
#    print $chr,"\t",$chrhash{$chr},"\n";
#}
print OUT "Total\t",$totalreads,"\n";
close OUT;

