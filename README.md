MONKEY
----
Version: 2.0 | Last Updated: May 2006

**MONKEY** is a set of programs designed to search alignments of non-coding DNA sequence for matches to matrices representing the sequence specificity of transcription factors.

![Alt text](/_img/monkey.png?raw=true "AA Probabilities")


## Published Papers

AM Moses, DY Chiang, DA Pollard, VN Iyer, MB Eisen. [*MONKEY: identifying conserved transcription-factor binding sites in multiple alignments using a binding site-specific evolutionary model.*](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2004-5-12-r98) Genome biology, 2004.

AM Moses et al. [*Large-scale turnover of functional transcription factor binding sites in Drosophila.*](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0020130) PLoS computational biology, 2006.

## Contents

1. Introduction
2. Getting started with monkey
3. The matrix file
4. The alignment file
5. The tree file
6. The .occ output file
7. More options
8. rMonkey


### I introduction

This program is designed to search aligned non-coding DNA for matches to a known
motif matrix.   It runs on multiple alignments (2 or more sequences) of non-coding DNA
and provides a variety of information of about the putative motif instances that
it encounters.  

The basic operation of the program is to scan a single 'reference' genome (which can be
specified by the user) ignoring gaps.  For all regions of the 'reference' genome that may
be instances of the motif, the corresponding region of the alignment is examined using
a model of motif evolution.  The sequences from each genome corresponding to this
region are reported, along with the various scores and p-values that have been computed.  
Monkey uses a statistical model to report a pvalue that is an estimate of the probability
of the observed match in random genomic alignments.


### II Getting started with monkey

To install monkey on linux

type

```
tar -xvf monkey2.0.tar
```

this will create a directory called monkey2.0 .  Change to that directory, and
compile the program.

```
$ cd monkey2.0
$ make all
```

Monkey can now be run by typing ./Monkey

Monkey requires two command line arguments: the names of the matrix file and the fasta
alignment file. The formats for these will be discussed below. The simplest
run of monkey might look like  

```
$ ./Monkey gal4sites.matrix example/YBR018C_upstream.aln.fasta
```

which scans the alignment example/YBR018C_upstream.aln.fasta for sequences
matching the matrix in gal4sites.matrix.  The results should appear in a new
file called example/YBR018C_upstream.aln.occ.

Monkey can be run on windows under cygwin (http://www.cygwin.com/). Once the
make and g++ compiler have been installed, just follow the instructions for
linux.


### III The matrix file

Monkey must read in the frequency matrix representing the sequence specificty of the motif.
The format for this 'matrix' file is simply four columns of numbers, with each row
representing the frequencies of A, C, G and T for one position in the motif.
Frequencies should be greater than 0 and sum to one. An example matrix file is
provided with the distribution.


### IV The alignment file

Monkey reads in a 'multiple fasta alignment.'  The format for file is as a fasta file, with
gapped positions represented by '-' and nucleotides represented as A, C, G and
T. Other characters, such as 'a','c','g','t' or 'N' may be included, but
columns containing these will be ignored. Examples of these files are included with the
distribution and the Monkey distribution includes a perl script (clustal2fasta.pl) to convert
clustalw format files to multiple fasta. Typing

```
$ perl scripts/clustal2fasta.pl filename.aln
```

will create a multiple fasta file called filename.aln.fasta.


### V The tree file

If the fasta alignment file contains more than 2 sequences a tree file is required to
specify the tree topology and branchlengths. This can be estimated using one of many
availble phylogentic analysis packages (such as paml).  The tree file is structured as
follows

```
t 7
r 6     n 5     l 0
n 5     n 4     l 2
n 4     l 1     l 3
b
4       0.09787 0.13253
5       0.12665 0.23934
6       0.154675        0.154675
```

The first section of the file, beginning with t, specifies the tree topology. t is followed by
the total number of nodes in the tree. r stands for the root, n for interior nodes and l for
leaves. The line

```
r 6     n 5     l 0
```

means that the root is node 6 and it joins the internal node 5 with the leaf 0.  The numbers of
the leaves MUST correspond to the order of the sequences in the corresponding multi-fasta
alignment file, 0 being the first sequence, 1 the next and so on.  The numbers of the other nodes
must be self-consistent, but are arbitrary.

The second second section of the tree file starts with b and lists the branch lengths for the
interior nodes.  The line

```
4       0.09787 0.13253
```

means that the left branch of node 4 (to leaf 1) is of length 0.09787 substitutions per site, and
the right branch (to leaf 3) is of length 0.13253.  Phylogenetic trees can be computed by many
freely available programs (such as the PAML package.)


The distribution contains a perl script (dnd2tree.pl) to convert 'Newick' format files to the tree
format used for monkey.  This script requires that the species occur in left
to right order corresponding to the top to bottom order of the fasta file.  In
addition, branchlengths must be specified in 'X.X' format, (such that 0 must be
written as 0.0 or 0.000), the tree be strictly bifurcating (including the
root) and does not support internal node names.

Monkey2.0 includes a parsing script 'nhx2tree.pl' to convert newick files to
the format used for monkey.  This script should handle the full specification of
the extended newick format (http://www.genetics.wustl.edu/eddy/forester/NHX.html)
and will produce warnings if the tree is not strictly bifurcating, or does not
contain branchlengths (requirements for monkey and Rmonkey).  As with dnd2tree.pl,
the left to right order of species names in the tree should agree with the top
to bottom order in the multi-fasta file.  


If the alignment contains more than 2 sequences,  Monkey will by default
look for a file with the same name as the multiple fasta file, except with the extension '.tree'.
It is possible to specify a tree file using the -tree option. For example, typing

 ```
$ ./Monkey gal4sites.matrix example/YBR018C_upstream.aln.fasta -tree all.tree
```

uses all.tree instead of example/YBR018C_upstream.aln.tree

For pairwise alignments, no tree is required and Monkey will by default calculate the distance between
the two sequences.  If -tree is used it should be followed by a number - monkey will
take this to be the distance between the sequences (in substituions per site.)



### VI The .occ output file

Monkey will produce a file with the same name as the multiple fasta alignment file that
was given as input, except the extension will now be .occ instead of .fasta.  This
file is in tab delimited format.  Each row in this file represents a putative instance of
the motif. It is possible to specify a different filename using the -out option.  For example,

```
$ ./Monkey gal4sites.matrix example/YBR018C_upstream.aln.fasta -out file.occ
```

will write the output to file.occ.

The first 7 columns in the .occ file contain 1) the name of the fasta file and
matrix file, 2) a probabilistic estimate of which strand the binding site is on,
where 1 is the forward strand, 0 is the reverse, 3) the position in the
alignment, 4) the position in the ungapped version of the reference sequence,
5) the segment of reference sequence, 6) the single sequence liklihood ratio
score and 7) the p-value associated with that score.  

The next group of columns contain the sequences and single sequence liklihood
ratio scores for all of the other sequences in the alignment.  

Finally, the last three columns contain the p-value associated with the score
under the evolutionary model, the % identitiy of the subsequences in the
alignment, and the % ungapped.  If the percent ungapped is high, monkey will
heuristically modify the alignment and calculate a p-value.  If the percent
gaps is large, monkey will print out '-'s instead of the sequences and
p-values.


### VII More options

`-help`
To see all the options currently available type

```
$ ./Monkey -help
```

`-m`
Monkey allows the use of several scores based on different assumptions about
the evolution of the motif. The default is to use the Halpern-Bruno model.

`-b`
This specifies the background substitution model.  The Default is to use the
Jukes-Cantor (JC) model, but the HKY model is also available.  If the HKY
model is used (by adding -b HKY to the command line) the transition
transversion rate ratio, or kappa, should also be specified using the -kappa
option.  The default is kappa=2.0 .

`-list`
If the same tree is to be used for searches of many alignment files, there is
a considerable speed advantage to scanning them all at once.  Monkey allows
the user to provide a file containing a list of alignments.

```
$ ./Monkey gal4sites.matrix eg.list -tree all.tree -list
```

will scan all of the alignments listed in eg.list using the tree found in
all.tree, and write the output to a single file called eg.list.occ.

`-scr`
Prints out the evolutionary score as well as the associated p-value.

`-ref`
specifies which sequence in the alignment file is to be treated as the
reference.  This is done by adding -ref n, where n is the position in the file
of the desired reference sequence, starting with 0.

`-cut`
This option can be used to reduce the total amount of hits produced in a scan
of many alignments.  This is the minimum score that a binding sites must have
in the reference genome before it is considered for evolutionary analysis.
The default is to use a score of zero.

`-freq`
To specify a background frequency file use -freq followed by the filename
This file has the following form.

```
A 0.3
C 0.2
G 0.2
T 0.3
```

This means that average frequencies in the intergenic region for A, C, G and T are 30%, 20%
20% and 30% respectively.  If a file is not specified monkey will calculate these frequencies from the sequences given.


### VIII rMonkey


Rmonkey is another program to detect matches to specificity matrices is multiple sequence alignments.  


Differences between rMonkey and the original monkey.  

rMonkey uses a different heuristic to identify binding sites in alignments (see
Moses et al 2006).  Breifly, rather than adjusting the alignment based on the
number of gaps in the cloumns that correspond to the binding site, rmonkey
selects an reference match, and then choses the best overlapping matches to
the matrix in each other sequence.  This makes rmonkey more 'agressive' at
assigning orthology to binding sites, and it can often consider non-orthologous
sequences to be aligned binding sites,  Unlike the original monkey, in which the
'reference' species was held constant for a given run, rmonkey includes an option
(specified by adding -any to the command line) which will allow the 'reference'
sequence for matches to be in any species.

In addition to the p-value assusiated with the evolutionary generalization of the
likelihood ratio score (as in monkey), rmonkey computes p-values associated with
a conditional likelihood ratio, which assumes a match in one 'reference' species
has already been observed (see Moses et al, 2006).   These p-values can be take
a long time and use a lot of memory to calculate; they can be disabled using the
-M option.

A final difference between rmonkey and monkey is that cutoffs in rmonkey are in
p-values rather than likelihood scores, and p-values are returned for each sequence
by default.  Single sequence p-values are computed by a different method, and will
differ slightly (but be more accurate) than those produced by the original monkey.  


#### Inputs

The input files to rmonkey are the same as for the original monkey, but some of the options have slightly different behaviours in rmonkey.

`-scr`
produces scores for each sequence and the evolutionary log likelihood ratio instead
of p-values

`-cut`
specifies a single sequence p-value cutoff for a match to be included in
evolutionary analysis


#### Outputs

The rmonkey output file is very similar to that of the original monkey.  rMonkey
produces a tab delimted file with one row corresponding to each match to the matrix.  
The first 5 columns contain the names, strand,  aligned and ungapped positions in
the reference sequence, and the sequence for the match in the reference.  Note that
positions in the aligned sequences are now given as 'start'-'stop'.  The next
column gives the single species p-value for the match in the reference. The next
sets of columns give positions, sequence and single species p-values for each of
the orthologous matches.  The next column gives the p-value associated with the
likelihood ratio score for a conserved match.  This is the p-value computed in
the original monkey.  If a sequence contains too many gaps, rmonkey reports a 'g'
instead of the '-' character used in the original monkey.

Finally, the last three columns give the conditional likelihood ratio (T statistic,
see Moses et al. 2006), and associated p-values under either the HB null hypothesis
(a test for lack of conservation) or the HKY null hypothesis (a more conservative
test for conservation).


**When should you use rMonkey?**

-If you are trying to find all the binding sites that might be conserved (and
are willing to include some that are actually not conserved)

-If you are trying to be conservative in identifying non-conserved bindng sites.

**When should you not use rMonkey?**

-If you are doing an evolutionary analysis of the bases in the binding sites -
in rmonkey the alignment is modified in a way that depends on the matrix.

-If you are interested in a region that contains many partially overlapping
binding sites, and want to preserve the orthology relationships between them.



#### Copyright Statement

version 2.0
copyright 2004-2006 by Alan Moses. This software is provided "as is" without
warranty of any kind.  The author assumes no responsibility for the results it
produces or conclusions based thereupon.  It is distributed free of charge for
academic use only.  Permission to copy and use it is granted free of charge
provided that no fee is charged and this copyright notice is not removed.
