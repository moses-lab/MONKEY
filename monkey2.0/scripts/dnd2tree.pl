#!/usr/bin/perl
die "usage: [.dnd file] ...\n" unless (@ARGV);
while (@ARGV) {
my $file = shift @ARGV;
my $treestring;
open (FILE, "<", $file);
while (my $line = <FILE>) {
	$treestring=$line if ($line =~m/;/);
}
close FILE;

my %tree;
die "couldn't parse $treestring\n" unless (parse_rooted_tree_with_branchlengths ($treestring, \%tree));

print tree_string(\%tree,root(\%tree)), " found\n";

my $treefilename=$file.".tree";
print "wrote $treefilename\n" unless (print_tree_file(\%tree,$treefilename));
}

sub print_tree_file {
	my ($t, $fname) = @_;
	open (OUT, ">", $fname);
	print OUT "t ";
	print OUT $$t{Nnodes};
	print OUT "\n";
	
	for (my $n=$$t{Nnodes}-1;$n>=0;$n--) {
		next if ($$t{$n}{leaf}==1);
		if ($$t{$n}{root}==1) { print OUT "r $n"; }
		else { print OUT "n $n"; }
			
		for my $br ("left","right") { 
			my $chn=$$t{$n}{$br};
			if ($$t{$chn}{leaf}==1) { print OUT "\tl $chn"; }
			else { print OUT "\tn $chn"; }
		}
		print OUT "\n";
	}
	print OUT "b\n";
	for (my $n=0;$n<$$t{Nnodes};$n++) {
		next if ($$t{$n}{leaf}==1);
		print OUT "$n ";			
		for my $br ("left","right") { 
			my $chn=$$t{$n}{$br};
			print OUT "\t".$$t{$chn}{branch};
		}
		print OUT "\n";
	}
	close OUT;		
	return 1;
}

sub parse_rooted_tree_with_branchlengths {

	my ($str, $t) = @_; 
	#print $str, "\n";
	$str =~ s/\s//g;
	my $nn=0;
	while ($str =~ m/(\w+):(\d+\.\d+)/g) {
		#print "leaf $1, branch $2\n";
		$$t{$nn}{"name"}=$1;
		$$t{$nn}{"branch"}=$2;
		$$t{$nn}{"leaf"}=1;		
		$nn++;
	}
	$$t{Nnodes} = 2*$nn - 1;
	#print keys %{$t}, "\n";
	my $i=0;
	#print "$nn leafs found\n";
	while ($nn < $$t{Nnodes}) {
		LOOP: for my $l (0 .. $nn-1) {
			next if ($$t{$l}{anc} ne "");
			for my $r ( 0 .. $nn-1) {
				next if ($r == $l);
				next if ($$t{$r}{anc} ne "");
				my $test = "(".$$t{$l}{"name"}.":".$$t{$l}{"branch"}.",".$$t{$r}{"name"}.":".$$t{$r}{"branch"}.")";
				#print "$test\n";
				
				for my $p (0 .. length($str) - length($test) -1) {
					my $match = substr($str,$p,length($test));
					if ($match eq $test) {
						#print "$test $match\n$str\n";
						$$t{$nn}{"name"} = $test;
						$$t{$nn}{left}=$l; $$t{$nn}{right}=$r;
						$$t{$l}{anc}=$nn; $$t{$r}{anc}=$nn;
						my $brstr = substr($str,$p+length($test),length($str)-length($test)-$p);
						$$t{$nn}{branch}=$1 if ($brstr =~ m/:(\d+\.\d+)/);	
						if ($brstr eq ";") {
							$$t{$nn}{root}=1;	
						}
						#print $$t{$nn}{branch}, "\n";
						$nn++;
						last LOOP;
					}	
				}
			}
		}
		$i++;
		if ($i>100) {
			print "couldn't parse\n"; return 0;
		} 		
	}
	return 1;
}


sub root {
	my ($t) = (@_);
	for my $n ( 0 .. $$t{Nnodes} - 1 ) {
		return $n if ($$t{$n}{root}==1);
	}
}

sub tree_string {
	my ($t, $n) = @_;
	
	return $$t{$n}{"name"} if ($$t{$n}{"leaf"}==1);
	
	return "(".tree_string($t,$$t{$n}{left}).",".tree_string($t,$$t{$n}{right}).")";
}
