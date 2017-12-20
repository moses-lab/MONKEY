#!/usr/bin/perl

#(c) 2006 by Alan Moses 
my $MIN_BRANCH=0.001;

die "usage: [NHX file]\nconverts nhx format to the monkey tree format \n" unless (@ARGV);
TREE: while (@ARGV) {
my $file = shift @ARGV;

my $trstr=""; 
open (FILE,"<", $file) or die "$file not found \n";
while (my $line=<FILE>) { chomp $line; $trstr .= $line; }
close FILE;
$trstr=~s/\s//g;
my %t; 
unless (parse_nhx2($trstr,\%t)) { print STDERR "$file not parsed\n"; next TREE; }
print STDERR "parsed: ".tree_string(\%t,get_root(\%t)), "\n";
my $trif=0; my $branchless=0;
for my $n (keys %{$t{node}}) {
	$trif++ if (scalar(@{$t{child}{$n}})>2);
	next if ($t{root}{$n}==1);
	if (($t{br}{$n} eq "") || ($t{br}{$n}<=0)) { $branchless++; $t{br}{$n}=$MIN_BRANCH; }
}
if ($trif>0) { print STDERR $trif. " multifucations (>2 leafs per node) found ...\n";	next TREE; }
if ($branchless>0) { 
	print STDERR $branchless. " nodes had no branches. (set to $MIN_BRANCH)\n"; 
}

print STDERR $file.".tree written \n" if (print_tree_file(\%t, $file.".tree"));
}

sub print_tree_file {
	my ($t, $fname) = @_;
	open (OUT, ">", $fname);
	print OUT "t ";
	print OUT $$t{Nnodes};
	print OUT "\n";
	
	for (my $n=$$t{Nnodes}-1;$n>=0;$n--) {
		next if ($$t{leaf}{$n}==1);
		if ($$t{root}{$n}==1) { print OUT "r $n"; }
		else { print OUT "n $n"; }
			
		for my $br ("1","0") { 
			my $chn=$$t{child}{$n}[$br];
			if ($$t{leaf}{$chn}==1) { print OUT "\tl $chn"; }
			else { print OUT "\tn $chn"; }
		}
		print OUT "\n";
	}
	print OUT "b\n";
	for (my $n=0;$n<$$t{Nnodes};$n++) {
		next if ($$t{leaf}{$n}==1);
		print OUT "$n ";			
		for my $br ("1","0") { 
			my $chn=$$t{child}{$n}[$br];
			print OUT "\t".$$t{br}{$chn};
		}
		print OUT "\n";
	}
	close OUT;		
	return 1;
}


sub parse_nhx2 {

	my ($str,$t) = @_; my $PMESS=0;
	return 1 unless (($str =~ m/\,/)&&($str =~ m/;/)); #trees must have more than one node.
	warn $str if ($PMESS); 
	$$t{Nnodes}=0;
	#find the leaves - the must have either a '(' or ',' in front and a ')' or ',' behind, and contain no '(', 
							#but essentially everythign else
	my $stuff = '[a-zA-Z0-9\.\=\_\:\s\&\$\[\]\-\{\}\/\%\!\'\+\?]';	
	my %NH;
	#my $stuff2 = '[a-zA-Z0-9\.\=\_\:\s\&\$\[\]\-\{\}\/\,]';					
	while ($str =~m/[\(\,]($stuff+)[\,\)]/g) { #find the leaves		
		my $ln = $$t{Nnodes};
		$$t{leaf}{$ln} = 1;
		create_node($1,$t);
		$NH{$$t{node}{$ln}}=$ln;
		pos($str) = pos($str)-1; # leaf match can overlap by 1
		print "$ln $1: ".$$t{name}{$ln}." ".$$t{br}{$ln}."\n" if ($PMESS);
	}
	#die;
	my $iter=1; my $max_ITER=1000; 
	while ($iter<$max_ITER) {
		
		my $nodesearch=0;
		NODE: for my $n (sort numerically keys %{$$t{node}}) { # find ancestors - nodes that contain this one
			next if (($$t{anc}{$n} ne "")||($$t{root}{$n}==1));
			my $nstr=$$t{node}{$n};
			#print "looking for ancestor for $n: \n".$$t{node}{$n}."\n...\n";
			$nodesearch++;
			$nstr =~ s/\W/\./g; #print "$nstr\n";
			next unless ($str =~ m/[\(\,]$nstr[\,\)]/);
			my $unresolved=1; my @ch = ($n); my $looptime=0;
			while ($unresolved) {
				$looptime++;
				last if ($looptime>10); #maximum number of daughters for each node
				if ($str =~ m/\($nstr\,(.+)/) {
					my $rest=$1;
					print "looking for $n rest=".substr($rest,0,10)."...\n" if ($PMESS);				
					my $nextnode = next_node_right($rest);
					my $nnstr = $nextnode;	$nnstr =~ s/\W/./g;
					print "next node is $nextnode\n" if ($PMESS);
					die "problem with $nextnode\n" unless ($nextnode);
					$nextnode .=$1 if ($str=~m/\($nstr\,$nnstr($stuff*)/);	
					
					$nnstr= $nextnode; $nnstr =~ s/\W/./g;				
					
					my $rn; #check the new sister
					if ($NH{$nextnode}) {
						$rn = $NH{$nextnode};
						print "old node $rn\n" if ($PMESS);
					}	
					else {
						$rn= $$t{Nnodes};
						#$nextnode = $1 if ($str =~m/$nextnode
						create_node($nextnode,$t);
						$NH{$nextnode}=$rn;
						print "new node $rn\n" if ($PMESS);
					}
					push @ch, $rn;
					#check the ancestor
					$nstr = $nstr."\,".$nnstr;
					if ($str =~ m/(\($nstr\)$stuff*)/) {
						$unresolved=0; my $nn; my $newn=$1;
						if ($NH{$newn}) {
							$nn = $NH{$newn};							
						}
						else {	
							$nn= $$t{Nnodes};
							create_node($newn,$t);
							$NH{$newn}=$nn;
						}	
						for my $m (@ch) { $$t{anc}{$m}=$nn; }
						print "$iter : resolved $nn (@ch) $1\n" if ($PMESS);
						if ($$t{node}{$nn}.";" eq $str) {
							print "root found\n" if ($PMESS);
							$$t{root}{$nn}=1;
						}
						
					}
					else { print "$nstr unresolved\n" if ($PMESS); }
				}
				elsif ($str =~ m/\,$nstr[\,\)](.+)/) { $unresolved=0; }
				#elsif ($$t{node}{$n}.";" eq $str) { $$t{root}{$n}=1; $unresolved=0; }
				else { print "node$n $nstr not found\n" if ($PMESS); }
			}			
			
		}
		#die;
		$iter++;
		last if ($nodesearch==0);
	}
	my $rn = get_root($t);
	if ($rn ne "") {
		my $alldef=1;
		for my $n (keys %{$$t{node}}) {
			next if ($n ==$rn);
			if ($$t{anc}{$n} eq "") { 
				$alldef=0;
				print "no anc for $n\n" if ($PMESS);
			}
		}
		if ($alldef) {
			find_children($t);
				#print STDERR "parsed ...\n";
			return 1;
		}		
	}			
	return 0;
}

sub next_node_right {
	my $stuff = '[a-zA-Z0-9\.\=\_\:\s\&\$\[\]\-\{\}\'\/]';
	my ($rstr) = @_;
	my $lb=0; my $rb=0; my $pos;
	while ($rstr =~ m/([\(\)\,])/g) {
		my $brac=$1; #my $ext=$2;
		if ($brac eq '(') { $lb++; }
		elsif ($brac eq ')') {
			$pos=pos($rstr);			
			return substr($rstr,0,$pos-1) if ($lb==0);						 
			$rb++;
			return substr($rstr,0,$pos) if ($lb==$rb);			
		}
		elsif ($brac eq ',') {
			$pos=pos($rstr);
			#warn "middle node in multif ..." if ($lb==0);			
			return substr($rstr,0,$pos-1) if ($lb==0);
		}
	}
	warn "node find failed $pos and $rstr";
	return "";			
}


sub tree_string {
	my ($t,$n)=@_;
	#print $n, "\n";
	if ($$t{leaf}{$n}) {
		#print $$t{name}{$n}."\n";
		return ($$t{br}{$n} ne "") ? $$t{name}{$n}.":".$$t{br}{$n} : $$t{name}{$n};
	}
	else {
		my $str="(";
		for my $ch (@{$$t{child}{$n}}) {
			$str.="," unless (length($str)==1);
			my $temp=tree_string($t,$ch);
			$str.=$temp;
			#print "$n $ch $str ".$$t{br}{$n}."\n";
		}
		return $str.")".$$t{name}{$n}.":".$$t{br}{$n} if ($$t{br}{$n} ne "");
		return $str.")".$$t{name}{$n}.";" if ($$t{root}{$n} ==1 );
		return $str.")".$$t{name}{$n};
	}			
	
}

sub find_children {
	my ($t)=@_;
	#print "findign childern...\n";
	for my $n (keys %{$$t{node}}) { next if ($$t{leaf}{$n}==1);
		for my $m (keys %{$$t{node}}) { 
			push @{$$t{child}{$n}}, $m if ($$t{anc}{$m} ==$n)
		}
		print "$n ".$$t{name}{$n}." ".$$t{anc}{$n}." has no children \n" unless (@{$$t{child}{$n}});
		#for my $c (@{$$t{child}{$n}}) { print $c. " "; }
		#print "\n";
	}
	return 1;		
}

sub get_root {
	my ($t)=@_;
	return "" unless (scalar (keys %{$$t{root}})>0);
	for my $n (keys %{$$t{root}}) {
		return $n if ($$t{root}{$n}==1) ;		
	}
	return "";	
}

sub create_node {
	my ($nstr,$t)=@_;
	my $stuff = '[a-zA-Z0-9\.\=\_\:\s\&\$\[\]\-\{\}\'\/]';
	my $tagstuff='[a-zA-Z0-9\.\=\_\:\s\&\$\-\{\}\/]';
	my $num = $$t{Nnodes};
	$$t{node}{$num} = $nstr;
	$$t{tag}{$num}=$1 if ($nstr=~m/\[($tagstuff+)\]\Z/);	
	$nstr =~ s/\[.+?\]//g;
	#print $nstr."\n";
	$$t{name}{$num} = $1 if ( $nstr =~ m/($stuff+)\:[\d\.]+\Z/);
	$$t{br}{$num} = $1 if ($nstr =~ m/\S+\:([\d\.]+)\Z/);	
	$$t{name}{$num} = $nstr if (($$t{name}{$num} eq "")&&($$t{leaf}{$num}==1));
	#if ($$t{br}{$num} eq "") {
	#	warn "no branch found for $nstr\n" if ($PMESS);
	#}	
	#print "$num ".
	$$t{Nnodes}++;
	return 1;
}

sub numerically { $a <=> $b }
