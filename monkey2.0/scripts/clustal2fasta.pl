#!/usr/bin/perl
die "usage: [.aln file] ...\n" unless (@ARGV);
while (@ARGV) {
my $file=shift (@ARGV);
my %ALN;
open (IN, "<", $file) or die "$file not found";
while (my $line=<IN>) { #print $line;
	next if ($line=~m/CLUSTAL/);
	if ($line =~ m/(\S+)\s+([-A-Z]+)/) {
		#print "$1\t$2\n";
		chomp $line;
		my $sp=$1;
		my $seq=$2;
		$ALN{$sp} .= $seq;
	}
	#elsif ($line) { print $line; }
}
close IN;

open (OUT, ">", $file.".fasta");
for my $sp (keys %ALN) {
	print OUT ">$sp\n$ALN{$sp}\n";
}
close OUT;
}
