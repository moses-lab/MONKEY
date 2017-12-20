using namespace std;

#include <iostream.h>
#include <fstream.h>
#include <math.h>
//#include <time.h>
//#include <stdio.h>
#include <stdlib.h>
//#include <map>
#include <string>

const int MAX_ALIGN_SIZE=15; //maximum number of species
const int MAX_ALIGNMENTS = 10000;
//const int MAX_W=20; // biggest motif size
double INITIAL_BRANCH=0.1;
bool PRINT_TO_SCREEN=1;
const bool BOTH_STRANDS=1;

const double MOTIF_PRIOR=0.001;
double thr=0.0;

const double MIN_RATE=0.001;
const int MAX_HITS=10000;
const int EVOL_MODEL=2; // 0= JC bg rate,  1= JC 1/2 bg rate 2=HB selection
const int N_MODELS=5;
double FUDGE=1.0; //correction because the HB model is based on the synonymous rate
double N_BINS=200.0;

struct sequence { // sequence object - 2 copies of the sequence, a string, and as ints
	int length;
	string name;
	string bases;
	int * asints;

	sequence();
	bool integerize();
	void erase();
	bool print(int start, int stop);
	bool is_base(int pos);
	//double Pr(int start, int stop, double * f);
};

struct node {
	bool leaf; bool root;
	struct node * anc;
	struct node * left;
	struct node * right;
	double abr,lbr,rbr;
	//double * seq_post;
	int base; int number;

	node();
	bool print();
};

struct tree { //hold a list of nodes, a pointer to the root and a some substitution matrices
	struct node * nodelist;
	int Nnodes;
	struct node * root;
	double ** psubs[4][4];

	tree();
	bool read_from_file(string treefile);
	bool initialize_pairwise(double dis);
	bool print(struct sequence * seq);
	bool compute_PJC();
	bool compute_PHKY(double * f, double kap);
	string tree_string( struct node * n, struct sequence * seq);
};

struct alignment { // alignment - holds some sequences, positions of the gaps. tree
	int length;
	int Nseqs;
	int gapped;
	int bad;
	int leftend;
	int rightend;
	
	string name;
	struct sequence * seq;
	int * nextgap;
	int * nextbad;
	long int * asints;
	struct tree phyl;
	double treelength;
	
	alignment();
	bool read_bg_freq(string freqfile);
	bool read_fasta(string filename);
	bool print();
	float div();
	
	bool compute_nextgap();
	void integerize();
	
	double * bgfreq;
	double kappa;
	
	bool read_tree(string filename);
	
};


struct motif { // holds a pwm and vector of evol. rates
	int length;
	string name;
	double ** f;
	double pri;
	double * K;
	double br;
	motif * rc;
	double ***** psubs;
	int *** pssm;
	int * minscr;	
	double ** pvalue;
	//double ** score_to_LR;
	//struct tree * phyl;

	motif();
	//bool intialize_from_word(int * word, int start, int stop, double pseudocounts);
	bool initialize_neutral_motif(int w, double pseudocounts);
	bool read_from_file(string mfile);
	bool print();
	double halpern_bruno(int pos);
	struct motif make_rc();
	//bool optimize_by_EM( struct SequenceData * sd,bool update_prior, bool update_matrix);
	double Pr(int s[], int start);
};

double p_subtree2 (struct sequence * seq, struct node * n, struct tree * phy, int base, int pos);
double p_subtree_hb (struct sequence * seq, struct node * n, int base, int sp, int mp, struct motif * m);
double p_subtree_jc (struct sequence * seq, struct node * n,  int base, int pos, double k);

double PJC (int a, int b, double alpha);


int gap_count (int * seq, int start, int stop);
struct sequence remove_gaps (struct sequence * seq, int start, int size);
double pw_pid (struct sequence * seq1,struct sequence * seq2, int start1, int start2, int length);

int main(int argc, char * argv[]) {	
	if (argc<=1) { cout << "usage: [pwm file] [fasta alignment] (options) ... or try -help for more details" <<endl; return 0; }
	int read_alignment_list(struct alignment * ali, string file);
	void help();
	//cout << argc <<endl;
	int ref_num=0,mmod=2,bgmod=0, Nalign;
	
	string freqfile = "?";
	string mfile = "?";
	string treefile = "?";
	string fastafile = "?";
	string outfile = "?";
	
	bool HEURISTIC=1; bool PVALS=1 ; bool LLRS=0; bool LIST=0; bool EVALS=0; 
	double KAPPA = 2.0; bool DAN=0;
	
	struct alignment * ALI = new alignment [MAX_ALIGNMENTS];	
	struct motif * M = new motif [2]; 
	
	int i=1;
	while (i<argc) {
		//cout << argv[i] <<endl;
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'm') {
				if (argv[i+1][0] == 'm') mmod=0;
				else if (argv[i+1][0] == 'J') mmod=1;
				else if (argv[i+1][0] == 'H') mmod=2;
				else if (argv[i+1][0] == 's') mmod=3;
				else { cout << "didn't recognize "<< argv[i+1] <<endl; help(); return 0; }
				i++;
			}	
			else if (argv[i][1] == 'b') {
				if (argv[i+1][0] == 'J') bgmod=0;
				else if (argv[i+1][0] == 'H') bgmod=1;
				else { cout << "didn't recognize "<< argv[i+1]<<endl; help(); return 0; }			
				i++;
			}
			else if (argv[i][1] == 's') { LLRS=1; }
			else if (argv[i][1] == 'c') { thr = atof(argv[i+1]); i++; }
			else if (argv[i][1] == 'k') { KAPPA = atof(argv[i+1]); i++; }
			else if (argv[i][1] == 't') { treefile = argv[i+1];	i++; }
			else if (argv[i][1] == 'f') { freqfile = argv[i+1];	i++; }
			else if (argv[i][1] == 'r') { ref_num = atoi(argv[i+1]); i++; }
			else if (argv[i][1] == 'h') { help(); return 0; }
			else if (argv[i][1] == 'l') { LIST=1; }
			else if (argv[i][1] == 'u') { HEURISTIC=0; }
			else if (argv[i][1] == 'e') { EVALS=1; }
			else if (argv[i][1] == 'o') { outfile=argv[i+1]; i++; }
			else if (argv[i][1] == 'D') { cout << "hey dan"<< endl; DAN=1; i++; }
			else { cout << argv[i] << " is not a recognized option"<<endl; help(); return 0; }
			
		}
		else if (mfile[0] == '?') { 
			mfile = argv[i]; 
			if (!M[1].read_from_file(mfile)) { 
				cerr <<"matrix file not read\n"; 
				return 0; 
			}
		}
		else if (fastafile[0] == '?') { fastafile = argv[i]; }
		else { cout << " too many arguments: "<<argv[i]<<endl; help(); return 0; }
				
		i++;
	}
	
	if (LIST) {
		Nalign = read_alignment_list(ALI, fastafile);
		if (Nalign==0) return 0; 
		else if (PRINT_TO_SCREEN) cerr << Nalign << " alignments read"<<endl;
	}
	else {
		Nalign=1;	
		if (!ALI[0].read_fasta(fastafile)) { cerr << "error reading alignment "<< fastafile<<endl; return 0; }
	}	
		
	int a;		
	//return 0;
	//cout << ALI[0].name <<endl; return 0;
	for (a=0;a<Nalign;a++) {		//cout << ALI[a].name<<endl;
		if (ALI[a].Nseqs>2) {		 			
			string temp;
			if (treefile[0]=='?') temp=ALI[a].name + ".tree";
			else temp = treefile;
			if (!ALI[a].read_tree(temp)) { cerr << "tree file "<<temp<<" not found"<<endl; return 0; }		
		}
		else if (treefile[0]!='?') {
			float k= atof(treefile.c_str());
			//if (PRINT_TO_SCREEN) cout << "using k="<< k<<endl;
			ALI[a].phyl.initialize_pairwise(k);
		}
		if (freqfile[0]!='?') ALI[a].read_bg_freq(freqfile);
		//else ALI[a].print();
	}
	
	//if (ALI[0].Nseqs>2) ALI[0].phyl.print(ALI[0].seq);
	if (ALI[0].Nseqs>2) 
		cerr << "tree interpreted as: " <<ALI[0].phyl.tree_string(ALI[0].phyl.root,ALI[0].seq) << endl;
	M[0].initialize_neutral_motif(M[1].length,0.0); //make a "motif" for the background model

	int p,m,s,b;

	struct motif * RC= new motif [2];
	for (m=0;m<2;m++) {
		RC[m]=M[m].make_rc();
		M[m].rc = &RC[m];
		M[m].pssm = new int ** [N_MODELS];
		M[m].pvalue = new double * [N_MODELS];
		M[m].minscr = new int [N_MODELS];
		//M[m].score_to_LR = new double * [N_MODELS];
		RC[m].pssm = new int ** [N_MODELS];
		RC[m].pvalue = new double * [N_MODELS];
		RC[m].minscr = new int [N_MODELS];
		//RC[m].score_to_LR = new double * [N_MODELS];
	}
	
	for (p=0;p<M[1].length;p++) for (b=0;b<4;b++) M[0].f[p][b]=ALI[0].bgfreq[b];
			
	//double logLR (struct motif * mod, struct alignment * ali, int pos);
	double logLR2 (struct motif * m, struct sequence * seq, struct tree * phy, double * bgf, int pos, int e_mod);
	double skip_gap_log_LR (struct motif * M, struct motif * B, struct sequence * seq, int start);
	double keep_gap_log_LR (struct motif * M, struct motif * B, struct sequence * seq, int start);
	struct sequence remove_gaps (struct sequence * seq, int start, int size);
	double scr;//double rcscr;
	int ugp=0;
	bool compute_pHB (struct motif * m, struct tree * t, int bg, double * bgf, double kap);
		
	for (a=0;a<Nalign;a++) { 
		if (bgmod==0) ALI[a].phyl.compute_PJC();
		else if (bgmod==1) { 
			ALI[a].kappa=KAPPA;
			ALI[a].phyl.compute_PHKY(ALI[a].bgfreq, ALI[a].kappa);
		}			
	}
	if (mmod==2) {
		if (LIST) {
			if (treefile[0] != '?') {		
				if (!compute_pHB (&M[1],&ALI[0].phyl,bgmod,ALI[0].bgfreq, ALI[0].kappa)) return 0;
				compute_pHB (&RC[1], &ALI[0].phyl,bgmod,ALI[0].bgfreq, ALI[0].kappa);
			}		
			else { cerr << "run seperately"<<endl; return 0; }	
		}
		else { 
			if (!compute_pHB (&M[1],&ALI[0].phyl, bgmod,ALI[0].bgfreq, ALI[0].kappa)) return 0;
			compute_pHB (&RC[1], &ALI[0].phyl, bgmod,ALI[0].bgfreq, ALI[0].kappa);
		}		
	}	
	//return 1;
	void compute_pssm_pvals (struct motif * m, struct motif * bg, struct tree * t, int N, int e_mod, int ref);
	
	if (PVALS) {
		if (PRINT_TO_SCREEN) cerr << "working on model " << mmod << "...";
		compute_pssm_pvals(&M[1], &M[0], &ALI[0].phyl, ALI[0].Nseqs,mmod,ref_num);	
		if (PRINT_TO_SCREEN) cerr<<"done"<<endl;
		if (PRINT_TO_SCREEN) cerr << "working on model " << 4 << "...";
		compute_pssm_pvals(&M[1], &M[0], &ALI[0].phyl, ALI[0].Nseqs,4,ref_num);	
		if (PRINT_TO_SCREEN) cerr<<"done"<<endl;
	}	
	//return 1;	
	if (outfile[0] == '?') {
		if (LIST) outfile=fastafile+".occ";
		else outfile=ALI[0].name+".occ";
	}
	else if (PRINT_TO_SCREEN) cerr << "writing to "<< outfile<<endl;	
	//return 1;
	ofstream out(outfile.c_str()); //open the output file and write the headers
	
	out <<"NAME\tSTRAND\tPOS\tUG\t";
	out << ALI[0].seq[ref_num].name << "\tlog(LR)\t";
	if (EVALS) out << "E(LR>lr)\t";
	if (PVALS) out << "p(LR>lr)\t";
	for (s=0;s<ALI[0].Nseqs;s++) if (s!=ref_num) out << ALI[0].seq[s].name<<"\tlog(LR)\t"; 		
	
	if (EVALS) out << "E(s>S)\t";
	if ((LLRS)&&(PVALS)) out <<"S\tp(s>S)\tpID\tpUG\n";
	else if (LLRS) out <<"S\tpID\tpUG\n";
	else if (PVALS) out <<"p(s>S)\tpID\tpUG\n";
	else out <<"pID\tpUG\n";	
	
	//return 1;
	//int gap_count (int * seq, int start, int stop);
	double pid (struct sequence * seq, int ns, int start, int size);
	double pug (struct sequence * seq, int ns, int start, int size);
	struct sequence optimized_seq(struct sequence * ref, struct sequence * orth, int start, int length, int offset);
	double aligned_motif_pval(struct motif * m, long int * ali, int ns, int pos, int e_mod);
	double gapped_aligned_motif_pval(struct motif * m, struct motif * bg, struct alignment * ali, int ns, int pos, int e_mod, int ref); 
	long int col_to_int (struct sequence * seq, int p, int n);	
	long int gapped_complem_int (long int i, int l);
	void general_variant(long int mn, int ml, int * seq, int N);
	void generate_variant (long int mn, int ml, int * seq);
	
	int total_bp=0;
	
	for (a=0;a<Nalign;a++) {
		ugp=0;
		int * spos = new int [ALI[a].Nseqs];
		for (s=0;s<ALI[a].Nseqs;s++) spos[s]=0;
	for (p=0;p<ALI[a].length-M[1].length;p++) { 
		if (ALI[a].seq[ref_num].is_base(p)) {
			scr=skip_gap_log_LR(&M[1], &M[0], &ALI[a].seq[ref_num], p);
			
			if (scr>thr) {						
				double * sp_scr= new double [ALI[a].Nseqs];				
				
				sp_scr[ref_num]=scr;
				//sp_pos[ref_num]=p;
				int offset=gap_count(ALI[a].seq[ref_num].asints, p, p+M[1].length);
				out << ALI[a].name << " "<<M[1].name << "\t"; //"= "<<scr <<endl;
				
				double ns=1.0; double average=scr;
				//if (PRINT_TO_SCREEN) cout << p<<" "<<ALI[a].nextgap[p]<<endl;
				//if (PRINT_TO_SCREEN) for (s=0;s<ALI[a].Nseqs;s++) ALI[a].seq[s].print(p,p+M[1].length);				
				if (((ALI[a].nextgap[p]<=(p+M[1].length))||(ALI[a].nextbad[p]<=(p+M[1].length)))&&(HEURISTIC==1)) {
					struct sequence * bestalign = new struct sequence [ALI[a].Nseqs];
 					for (s=0;s<ALI[a].Nseqs;s++) {					
 						if (s == ref_num) {							
 							bestalign[s] = remove_gaps(&ALI[a].seq[s],p, M[1].length); 
 							double strand=M[1].Pr(bestalign[s].asints,0)/(M[1].Pr(bestalign[s].asints,0)+M[1].rc->Pr(bestalign[s].asints,0));						
 							out <<strand<<"\t"<< p <<"\t"<<ugp<<"\t"<<bestalign[s].bases<<"\t"<<sp_scr[ref_num]<<"\t"; 						
 						}
 						else {						
 							bestalign[s] = optimized_seq(&ALI[a].seq[ref_num],&ALI[a].seq[s],p, M[1].length, offset);
 							sp_scr[s]=skip_gap_log_LR(&M[1], &M[0], &bestalign[s],0);
 						} 						
 					} 						 					
  					bool usall=1; int pos;
  					for (s=0;s<ALI[a].Nseqs;s++) if (bestalign[s].length!=M[1].length) usall=0;
 					if (usall) { 
						// double aligned_motif_pval(struct motif * m, struct motif * bg, struct alignment * ali, int ns, int pos, int e_mod) 
 						long int * asints = new long int [M[1].length];
  						for (pos=0;pos<M[1].length;pos++) asints[pos] = col_to_int(bestalign, pos, ALI[a].Nseqs);
						if (EVALS) out << ((double)(ALI[a].length - M[1].length))*aligned_motif_pval(&M[1],asints,ALI[a].Nseqs,0,4)<<"\t";
 						if (PVALS) out << aligned_motif_pval(&M[1],asints,ALI[a].Nseqs,0,4)<<"\t";						
 						//cout << bestalign[ref_num].bases<<endl;
						for (s=0;s<ALI[a].Nseqs;s++) if (s!=ref_num) {
						  //cout <<p<<" "<< bestalign[s].bases<<endl;
						  	if (DAN) out << spos[s] <<"\t";
  							out <<bestalign[s].bases<<"\t"<<sp_scr[s]<<"\t";
  							average+=sp_scr[s];
  							ns+=1.0;
  						}
  						int emo; 						
						if (EVALS) out <<((double)(ALI[a].length - M[1].length))*aligned_motif_pval(&M[1],asints,ALI[a].Nseqs,0,mmod)<<"\t";
  						for (emo=0;emo<3;emo++) if (emo==mmod) {							
  							if (LLRS) out << logLR2(&M[1],bestalign,&ALI[a].phyl,M[0].f[0],0,emo) <<"\t";
  							if (PVALS) out <<aligned_motif_pval(&M[1],asints,ALI[a].Nseqs,0,emo)<<"\t";
  						}
						if (mmod==3) {
  							if (LLRS) out <<average/ns<<"\t";
 							if (PVALS) out <<aligned_motif_pval(&M[1],asints,ALI[a].Nseqs,0,3)<<"\t";
						}
  						out <<pid(bestalign,ALI[a].Nseqs,0,M[1].length)<<"\t";
  						out <<pug(ALI[a].seq,ALI[a].Nseqs,p,M[1].length);						
 						out <<endl;
						delete asints;
  					}
  					else {
 						if (EVALS) out<<((double)(ALI[a].length - M[1].length))*gapped_aligned_motif_pval(&M[1],&M[0],&ALI[a],ALI[a].Nseqs,p,4, ref_num)<<"\t";
						if (PVALS) out<<gapped_aligned_motif_pval(&M[1],&M[0],&ALI[a],ALI[a].Nseqs,p,4, ref_num)<<"\t";
						
 						for (s=0;s<ALI[a].Nseqs;s++) {
 							if (s!=ref_num) {
 								if (bestalign[s].length==M[1].length) out <<bestalign[s].bases<<"\t"<<sp_scr[s]<<"\t";
 								else out<< "-\t-\t";
 							}	
 						}	
 						out << "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t"<<endl;	
						total_bp--;
 					}	
					//for (s=0;s<ALI[a].Nseqs;s++) bestalign[s].erase();
 					//bestalign=0;
					//delete bestalign;				  
				}
				else {
 					int pos;
					sp_scr[ref_num] = keep_gap_log_LR(&M[1],&M[0],&ALI[a].seq[ref_num],p);
 					double strand=M[1].Pr(ALI[a].seq[ref_num].asints,p)/(M[1].Pr(ALI[a].seq[ref_num].asints,p)+M[1].rc->Pr(ALI[a].seq[ref_num].asints,p));
 					out <<strand<<"\t"<< p <<"\t"<<ugp<<"\t"; for (pos=0;pos<M[1].length;pos++) out << ALI[a].seq[ref_num].bases[p+pos];
 					out<<"\t"<<sp_scr[ref_num]<<"\t";
					if (EVALS) out<<((double)(ALI[a].length - M[1].length))*gapped_aligned_motif_pval(&M[1],&M[0],&ALI[a],ALI[a].Nseqs,p,4, ref_num)<<"\t";
					if (PVALS) out<<gapped_aligned_motif_pval(&M[1],&M[0],&ALI[a],ALI[a].Nseqs,p,4, ref_num)<<"\t";
					
 					for (s=0;s<ALI[a].Nseqs;s++) if (s!=ref_num) { 
						//sp_scr[s]=skip_gap_log_LR(&M[1],&M[0],&ALI[a].seq[s],p);
						//sp_scr[s]=log(M[1].Pr(ALI[a].seq[s].asints,p)/M[0].Pr(ALI[a].seq[s].asints,p));
						sp_scr[s]=keep_gap_log_LR(&M[1],&M[0],&ALI[a].seq[s],p);
						if (DAN) out << spos[s] <<"\t";
 						for (pos=0;pos<M[1].length;pos++) out << ALI[a].seq[s].bases[p+pos];
 						out << "\t"<<sp_scr[s] <<"\t";
 						average+=sp_scr[s];
 						ns+=1.0;
 					}	
 					int emo;
					if (EVALS) out <<((double)(ALI[a].length - M[1].length))*gapped_aligned_motif_pval(&M[1],&M[0],&ALI[a],ALI[a].Nseqs,p,mmod, ref_num)<<"\t";
					if (mmod==3) {
						if (LLRS) out <<average/ns<<"\t";
						if (PVALS) out <<gapped_aligned_motif_pval(&M[1],&M[0],&ALI[a],ALI[a].Nseqs,p,3, ref_num)<<"\t";
					}	
 					for (emo=0;emo<3;emo++) if (emo==mmod) {
 						if (LLRS) out << logLR2(&M[1],ALI[a].seq,&ALI[a].phyl,M[0].f[0],p,emo) <<"\t";
 						if (PVALS) out <<gapped_aligned_motif_pval(&M[1],&M[0],&ALI[a],ALI[a].Nseqs,p,emo, ref_num)<<"\t";
 					}
					out <<pid(ALI[a].seq,ALI[a].Nseqs,p,M[1].length)<<"\t";
 					out <<pug(ALI[a].seq,ALI[a].Nseqs,p,M[1].length);					
					out <<endl;						
				}
			}
			ugp++;
			total_bp++;
		}
		for (s=0;s<ALI[a].Nseqs;s++) if (ALI[a].seq[s].is_base(p)) spos[s]++;	
	}
	}	
	//cout <<"log(LR)\n";
	out <<total_bp<<" bp searched"<<endl;
}


void help () {
	//cout << "help contents here"<<endl;
	cout << "you are using monkey v1.0"<<endl<<"usage: [matrix file] [fasta alignment] (options)"<<endl;
	cout << "\toption\texplanation"<<endl;
	cout << "\t-m\tmotif model : mk,JC0.5,HB or simple (default is HB)"<<endl;
	cout << "\t-b\tbackground model : JC,HKY (default is JC)"<<endl;
	cout << "\t-kappa\ttransition transversion rate ratio for HKY model (default is 2.0)"<<endl;
	cout << "\t-tree\ttree file name (default is alignment.tree) or pairwise distance (default is to compute it)"<<endl;
	cout << "\t-list\tsecond argument is a text file that is a list of fasta alignments (instead of a fasta alignment)"<<endl;
	cout << "\t-help\tdisplay this mesage"<<endl;
	cout << "\t-scr\tprint scores (as well as pvalues)"<<endl;
	cout << "\t-freq\tbackground frequency file (default is to compute them from the sequences)"<<endl;
	cout << "\t-ref\treference sequence number (used in heuristics; default is 0, the first one in the file)"<<endl;
	cout << "\t-cut\theuristic score cutoff for reference genome sites to test (default is 0.0)"<<endl;
	cout << "\t-out\t.occ output filename (default is fasta alignment.occ)"<<endl;
	 
	return;
}


int read_alignment_list(struct alignment * ali, string filename) {
	ifstream data(filename.c_str());	
	if (!data) { cout << filename<<" not found" << endl; return 0; }
	int length=-1;
	while (!data.fail()) {
		string line;
		getline(data,line);
		//cout << line << endl;
		length++;
	}
	if (PRINT_TO_SCREEN) cout << length << " alignemnts found" << endl;
	//ali = new alignment [length];
	ifstream data2(filename.c_str());
	string fastafile;
	int a;
	int not_read=0;
	for (a=0;a<length;a++) {
		data2>>fastafile;
		//cout << a << " "<< fastafile<<endl;
		if (!ali[a].read_fasta(fastafile)) not_read++;
	}	
	if ((PRINT_TO_SCREEN)&&(not_read>0)) cout << not_read << " alignemnts not read" << endl;
	return length;
}


 double aligned_motif_pval(struct motif * m, long int * ali, int ns, int pos, int e_mod) {
 	long int complem_int (long int i, int l);
 	int p, scr,rcscr;
	scr=rcscr=0;
	for (p=0;p<m->length;p++) if (ali[p+pos] !=-1) {
		scr+=m->pssm[e_mod][p][ali[p+pos]];
		rcscr+=m->pssm[e_mod][m->length-p-1][complem_int(ali[p+pos],ns)];
	}
	
	return ( (scr>rcscr) ? m->pvalue[e_mod][scr] : m->pvalue[e_mod][rcscr] );	
 }
 
 double gapped_aligned_motif_pval(struct motif * m, struct motif * bg, struct alignment * ali, int ns, int pos, int e_mod, int ref) {
 	long int complem_int (long int i, int l);
	int * C= new int [4]; C[0]=3; C[1]=2; C[2]=1; C[3]=0;
	
	//long int gapped_to_ungapped_int(long int gi, int l);
 	int p, scr,rcscr;
	scr=rcscr=0;
	//int gapscr =0; int rcgapscr=0;
	for (p=0;p<m->length;p++) {
		//long int ugint = gapped_to_ungapped_int(ali->asints[p+pos],ns);
		if (ali->asints[p+pos] !=-1) {
			scr+=m->pssm[e_mod][p][ali->asints[p+pos]];
			rcscr+=m->pssm[e_mod][m->length-p-1][complem_int(ali->asints[p+pos],ns)];
			//scr+=m->pssm[e_mod][p][ali[p+pos]];
			//cout << ali[p+pos]<<" "<<gapped_complem_int(ali[p+pos],ns)<<endl;
			//rcscr+=m->pssm[e_mod][m->length-p-1][complem_int(ali[p+pos],ns)];
		}
		else {
			//gapped=1;
			int b,s; double pm=0.0; double pbg=0.0; double prc=0.0; double llr=0.0; double rcllr=0.0;
			//double pbgrc=0.0;
			for (b=0;b<4;b++) {
				
				//double bcond = p_subtree_jc(ali->seq,ali->phyl.root,b,p+pos,1.0) ;
				double bcond = p_subtree2(ali->seq,ali->phyl.root,&ali->phyl,b,p+pos);
				pbg+=bg->f[p][b]*bcond;
			//	pbgrc+=bg->f[p][C[b]]*bcond;
				if (e_mod==0) { pm += m->f[p][b]*bcond; prc += m->rc->f[p][b]*bcond; }
				else if (e_mod==1) { 
					double mcond = p_subtree_jc(ali->seq,ali->phyl.root,b,p+pos,0.5);
					pm += m->f[p][b]*mcond; prc += m->rc->f[p][b]*mcond;
				}	
				else if (e_mod==2) {
					pm += m->f[p][b]*p_subtree_hb(ali->seq,ali->phyl.root,b,p+pos,p,m);
					prc += m->rc->f[p][b]*p_subtree_hb(ali->seq,ali->phyl.root,b,p+pos,p,m->rc);
				}	
			}
			if (e_mod<3) {
				scr += (int) (N_BINS*log(pm/pbg));
				rcscr += (int) (N_BINS*log(prc/pbg));
			}
			else if (e_mod==3) {				
				for (s=0;s<ns;s++) if (ali->seq[s].is_base(p+pos)) {
					int b = ali->seq[s].asints[p+pos];
					llr+=log(m->f[p][b]/bg->f[p][b]);
					rcllr+=log(m->rc->f[p][b]/bg->f[p][C[b]]);
				}
				scr += (int) (N_BINS*llr/((int) ns));
				rcscr+=(int) (N_BINS*rcllr/((int) ns));
			}
			else if (e_mod==4) {					
				if (ali->seq[ref].is_base(p+pos)) {
					int b = ali->seq[ref].asints[p+pos];
					llr=log(m->f[p][b]/bg->f[p][b]);
					rcllr=log(m->rc->f[p][b]/bg->f[p][C[b]]);
					//rcllr=log(m->f[m->length-p-1][C[b]]/bg->f[p][C[b]]);
				}
				scr += (int) (N_BINS*llr);
				rcscr+=(int) (N_BINS*rcllr);
				//int b;
				//cout << p <<" "<< (b=ali->seq[ref].asints[p+pos])<<" " << (((int) (N_BINS*llr)) - m->minscr[e_mod]) <<" "<<(scr-m->minscr[e_mod])<< endl;
			
				//cout << m->f[p][b] << " "<<m->rc->f[m->length-p-1][C[b]]<<endl;
			}		
			
			scr -= m->minscr[e_mod];			
			rcscr -= m->minscr[e_mod];
		}	
		//cout<< pos <<" " <<e_mod<< " "<<scr<< " "<<m->pvalue[e_mod][scr]<< " "<<rcscr<< " "<<m->pvalue[e_mod][rcscr]<<endl;
	}			
	//cout<< pos <<" " <<e_mod<< " "<<m->minscr[e_mod]<<" "<<scr<< " "<<m->pvalue[e_mod][scr]<< " "<<rcscr<< " "<<m->pvalue[e_mod][rcscr]<<endl;
	delete C;
	return ( (scr>rcscr) ? m->pvalue[e_mod][scr] : m->pvalue[e_mod][rcscr] );	
 }
 
 
 //void compute_pssm_pvals
 
 void compute_pssm_pvals (struct motif * m, struct motif * bg, struct tree * t, int N, int e_mod, int ref) {
 	void generate_variant (long int mn, int ml, int * seq);
	//void general_variant(long int mn, int ml, int * seq, int N);
	long int Nalpha = (long int) pow(4.0,N);
	//long int Nalpha = (long int) pow(5.0,N);
	m->pssm[e_mod] = new int * [m->length];
	
	int p,b,s;
	long int i;
	int min=10000;	
	//int worst=10000;
	double ** bgdist = new double * [m->length];
	//double 
	int best=0;
	for (p=0;p<m->length;p++) {
		//double check=0.0;
		m->pssm[e_mod][p] = new int [Nalpha];
		bgdist[p] = new double [Nalpha];
		for (i=0;i<Nalpha;i++) {
			struct sequence * col = new struct sequence [N];
			int * variant = new int [N];
			generate_variant(i,N,variant);
			//general_variant(i,N,variant,5);
			double llr=0.0;
			//bool gapped=0;
			for (s=0;s<N;s++) {
				if (e_mod==3) llr += log(m->f[p][variant[s]]/bg->f[p][variant[s]]);
				//if (variant[s]==4) gapped=1;
				col[s].asints = new int [1];
				col[s].asints[0]=variant[s];				
				//if (e_mod==4) cout << variant[s];
			}
			if (e_mod==3) m->pssm[e_mod][p][i] = (int) (N_BINS*(llr/((double) N)));
			if (e_mod==4) m->pssm[e_mod][p][i] = (variant[ref] == 5) ? 0 : ((int) (N_BINS*log(m->f[p][variant[ref]]/bg->f[p][variant[ref]])));
			//else {
			//if (e_mod==4) { cout << " |"<<p<<"| "; for (b=0;b<4;b++) cout << m->f[p][b] << " ";}
			double pm=0.0; 
			double pbg=0.0;			
			for (b=0;b<4;b++) {
				double bcond = p_subtree2(col,t->root,t,b,0) ;
				pbg+=bg->f[p][b]*bcond;
				if (e_mod==0) pm += m->f[p][b]*bcond;
				if (e_mod==1) pm += m->f[p][b]*p_subtree_jc(col,t->root,b,0,0.5);
				if (e_mod==2) pm += m->f[p][b]*p_subtree_hb(col,t->root,b,0,p,m);
				//if (e_mod==4) cout <<p_subtree_hb(col,t->root,b,0,p,m) <<" ";				
			}
			bgdist[p][i] = pbg;
			//if (e_mod==4) cout << endl	;
			//check += pbg;			
			if (e_mod<3) m->pssm[e_mod][p][i] = ((int) (N_BINS*log(pm/pbg)));
			//if (e_mod==2) cout << i << " "<<pm<<" "<<pbg<<" "<<m->pssm[e_mod][p][i]<<endl;
			//}	
			if (m->pssm[e_mod][p][i]<min) min=m->pssm[e_mod][p][i];
		}
		//cout << check << "=? 1"<<endl;
	}
	m->minscr[e_mod]=min;
	for (p=0;p<m->length;p++) { 
		int posbest=0;
		for (i=0;i<Nalpha;i++) { 
			m->pssm[e_mod][p][i]-=min;			
			if (m->pssm[e_mod][p][i]>posbest) posbest=m->pssm[e_mod][p][i];
		}
		best+=posbest;
	}	
	//if (e_mod==4) for (p=0;p<m->length;p++) { for (i=0;i<Nalpha;i++) cout << m->pssm[e_mod][p][i]<< " "; cout <<endl; }
	double ** pdf = new double * [m->length];
	for (p=0;p<m->length;p++) pdf[p]= new double [best+1];	
	for (p=0;p<m->length;p++) for (i=0;i<(best+1);i++) pdf[p][i] = -1.0;
	m->pvalue[e_mod] = new double [best+1];
	
	double score_pdf (int scr, int p, long int Na,double ** pbg, int ** pssm, double ** pdf);
	
	m->pvalue[e_mod][best] = score_pdf(best,m->length-1,Nalpha,bgdist,m->pssm[e_mod],pdf);
	//m->pvalue[e_mod][best] = score_pdf(459,2,Nalpha,bgdist,m->pssm[e_mod],pdf);
	for (s=best-1;s>=0;s--)  {		
		m->pvalue[e_mod][s] = score_pdf(s,m->length-1,Nalpha,bgdist,m->pssm[e_mod],pdf)
								+ m->pvalue[e_mod][s+1] ;
		//if (e_mod==4) cout << s <<" "<<score_pdf(s,m->length-1,Nalpha,bgdist,m->pssm[e_mod],pdf)<< " "<<m->pvalue[e_mod][s]<<endl;	
		//cout << s << " "<< pdf[2][s]<<endl;
	}	
	delete pdf;
	delete bgdist;							
 }
 
 double score_pdf (int scr, int p, long int Na, double ** pbg, int ** pssm, double ** pdf) {
 	if (scr<0) return 0.0;
 	if (p==-1) return ((scr==0) ? 1.0 : 0.0);
	//cout <<p<<" "<<scr<<" "<<pdf[p][scr]<<endl;
	if (pdf[p][scr] != -1.0) return pdf[p][scr];
	long int i;
	pdf[p][scr]=0.0;
	for (i=0;i<Na;i++) pdf[p][scr]+=score_pdf(scr-pssm[p][i], p-1, Na, pbg, pssm, pdf) * pbg[p][i];
	//cout <<"p("<< scr<<" after "<<p <<")";
	//cout << "="<< pdf[p][scr]<<endl;
	//delete i;
	return pdf[p][scr];
 }
 
 long int gapped_to_ungapped_int(long int gi, int l) {
 	void general_variant (long int mn, int ml, int * seq, int N);
	long int seq_to_num (int * seq, int start, int ml);
 	int * s = new int [l];
	general_variant(gi,l,s,5);
	int p; bool gapped=0;
	for (p=0;p<l;p++) if (s[p]==4) gapped=1;
	if (gapped) return -1;
	else return seq_to_num(s,0,l);
 }
 
 long int col_to_int (struct sequence * seq, int p, int n) {
 	long int seq_to_num (int * seq, int start, int ml);
 	int s; int * asints = new int [n];	
	for (s=0;s<n;s++) asints[s] = seq[s].asints[p];
	//for (s=0;s<n;s++) cout << asints[s] << " "<<seq[s].bases[p]<< " "<<seq[s].asints[p]<<endl;
	long int i = seq_to_num(asints,0,n);
	delete asints;
	return i;
 }
 
 long int gapped_col_to_int (struct sequence * seq, int p, int n) {
 	long int gapped_seq_to_num(int * seq, int start, int ml);
	int s; int * asints = new int [n];	
	for (s=0;s<n;s++) asints[s] = seq[s].asints[p];
	long int i = gapped_seq_to_num(asints,0,n);
	delete asints;
	return i;
 }
 
 long int complem_int (long int i, int l) {
 	if (i<=-1) return -1;
 	long int seq_to_num (int * seq, int start, int ml);
	void generate_variant (long int mn, int ml, int * seq);
 	int * seq = new int [l];
	generate_variant(i,l,seq);
	//long int * rc = new long int [l];
	int p;
	int * com = new int [4]; com[0]=3; com[1]=2; com[2]=1; com[3]=0;
	//int com[4];  com[0]=3; com[1]=2; com[2]=1; com[3]=0; 
	for (p=0;p<l;p++) seq[p] = com[seq[p]];	
	long int ci = seq_to_num(seq,0,l);
	delete seq;
	delete com;
	return ci; 
 }
 
 long int gapped_complem_int (long int i, int l) {
 	long int gapped_seq_to_num (int * seq, int start, int ml);
	void general_variant (long int mn, int ml, int * seq, int N);
 	int * seq = new int [l];
	general_variant(i,l,seq,5);
	//cout << i << " "<<l<<endl;
	//long int * rc = new long int [l];
	int p;
	int * com = new int [4]; com[0]=3; com[1]=2; com[2]=1; com[3]=0;
	//int com[4];  com[0]=3; com[1]=2; com[2]=1; com[3]=0; 
	///for (p=0;p<l;p++) cout << seq[p] << " "<<endl;
	for (p=0;p<l;p++) seq[p] = ((seq[p]==-1)||(seq[p]==4)) ? 4 : com[seq[p]];	
	//for (p=0;p<l;p++) cout << seq[p] << " "<<endl;
	//long int ci = seq_to_num(seq,0,l);
	long int ci = gapped_seq_to_num(seq,0,l);
	delete seq;
	delete com;
	return ci; 
 }
 
void general_variant (long int mn, int ml, int * seq, int N) {
	int p,b;
	for (p=0;p<ml;p++) {
		b = mn % N;
		seq[p] = b;
		mn /= N;
	}
}
 
long int gapped_seq_to_num (int * seq, int start, int ml) {
	long int mn=0;
	int p;
	for (p=(ml-1);p>=0;p--) { 
		mn *= 5;
		mn += ((seq[start+p]==-1) ? 4 : seq[start+p]);
	}
	return mn;	
}

void generate_variant (long int mn, int ml, int * seq) {
	//int * seq = new int [ml];
	int p,b;
	for (p=0;p<ml;p++) {
		b = mn % 4;
		seq[p] = b;
		mn /= 4;
	}
	//print_seq(seq,ml);
	//return seq;
}


long int seq_to_num (int * seq, int start, int ml) {
	long int mn=0;
	int p;
	for (p=(ml-1);p>=0;p--) {
		mn*=4;
		mn+=seq[start+p];
	}
	return mn;
}


bool compute_pHB (struct motif * m, struct tree * t, int bg, double * bgf, double kap) {
	//cout << m->name<<endl;
	void matrix_exponential (double ** a, double ** b, int n);
	m->psubs = new double **** [m->length];
	int p,a,b,n;
	for (p=0;p<m->length;p++) {
		m->psubs[p] = new double *** [t->Nnodes];
		for (n=0;n<t->Nnodes;n++) m->psubs[p][n] = new double ** [t->Nnodes];
		
		for (n=0;n<t->Nnodes;n++) if (!t->nodelist[n].root) {
			double ** R = new double * [4];	
			m->psubs[p][t->nodelist[n].anc->number][n] = new double * [4];
			for (a=0;a<4;a++) {
				R[a] = new double [4];
				m->psubs[p][t->nodelist[n].anc->number][n][a] = new double [4];
			}
				
			if (bg==0) {
				for (a=0;a<4;a++) {
					R[a][a]=0.0;
					for (b=0;b<4;b++) if (a != b) {					
						R[a][b] = (t->nodelist[n].abr/3.0);
						if (m->f[p][a] != m->f[p][b]) 
							R[a][b] *= (log(m->f[p][b]/m->f[p][a])/(1.0 - (m->f[p][a]/m->f[p][b])));
												
						R[a][a] -= R[a][b];	
					}	
				}
			}	
			else if (bg==1) {
				double ** Q = new double * [4];
				for (a=0;a<4;a++) {
					double tot=0.0;  Q[a] = new double [4];
					for (b=0;b<4;b++) if (a != b) {	
						Q[a][b] = bgf[b];
						if (((a==0)&&(b==2)) || ((a==2)&&(b==0)) || ((a==1)&&(b==3)) || ((a==3)&&(b==1)))
							Q[a][b] *= kap;
						tot+=Q[a][b];	
					}				
					for (b=0;b<4;b++) if (a!=b) Q[a][b] *= (t->nodelist[n].abr/tot);
				}
				for (a=0;a<4;a++) {
					double idrate=0.0;	
					for (b=0;b<4;b++) if (a!=b) {
						double temp = (m->f[p][b]*Q[b][a])/(m->f[p][a]*Q[a][b]);
						//double temp = 0.5;
						if (m->f[p][a] != m->f[p][b]) R[a][b] = Q[a][b] * (log(temp)/(1.0 - (1.0/temp)));
						else R[a][b] = Q[a][b];
						idrate -= R[a][b];
					}
					R[a][a] = idrate;						
				}	
			}						 
			matrix_exponential(R, m->psubs[p][t->nodelist[n].anc->number][n], 4);
//  			void print_matrix (double ** m, int n);
//  			for (a=0;a<4;a++) cout << m->f[p][a] << " ";
//  			cout << " "<<p<<" "<<t->nodelist[n].anc->number <<"->"<<n<<" k="<<t->nodelist[n].abr<<endl;
//  			cout << "rates:"<<endl; print_matrix(R,4);
//  			cout << "probabilities:"<<endl; print_matrix(m->psubs[p][t->nodelist[n].anc->number][n], 4);
// 			
		}		
	}
// 	for (a=0;a<4;a++) for (b=0;b<4;b++) psubs[a][b] = new double * [Nnodes];
// 	for (n=0;n<Nnodes;n++) for (a=0;a<4;a++) for (b=0;b<4;b++) psubs[a][b][n] = new double [Nnodes];
// 	for (n=0;n<Nnodes;n++) if (!nodelist[n].root) for (a=0;a<4;a++) for (b=0;b<4;b++) {
// 		if (a != b) psubs[a][b][nodelist[n].anc->number][n] = 0.25 - 0.25*exp(-(4.0/3.0)*nodelist[n].abr);
// 		else psubs[a][b][nodelist[n].anc->number][n] = 0.25 + 0.75*exp(-(4.0/3.0)*nodelist[n].abr);
// 	}
	
	return 1;
}

void print_matrix (double ** m, int n) {
	int a,b;
	for (a=0;a<n;a++) { 
		for (b=0;b<n;b++) cout << m[a][b]<< " ";
		cout << endl;
	}
}		

void matrix_exponential (double ** a, double ** b, int n) { //b= exp(a) a is n x n
	double temp[n][n];
	double temp2[n][n];
	int i,j,k;
	for (i=0;i<n;i++) for (j=0;j<n;j++) {
		b[i][j] = a[i][j]; 
		if (i==j) b[i][j]+=1.0;
		temp[i][j]=0.0;
		for (k=0;k<n;k++) temp[i][j]+=a[i][k]*a[k][j];
	}
	double fact=1.0;	
	double MAX_ITER=200;
	int iter;
	for (iter=2;iter<MAX_ITER;iter++) {
		fact /= ((double) iter);
		for (i=0;i<n;i++) for (j=0;j<n;j++) {
			b[i][j]+= temp[i][j]*fact;
			temp2[i][j]=0.0;
			for (k=0;k<n;k++) temp2[i][j]+=a[i][k]*temp[k][j];
		}
		for (i=0;i<n;i++) for (j=0;j<n;j++) temp[i][j] = temp2[i][j];
	}

}

struct sequence optimized_seq(struct sequence * ref, struct sequence * orth, int start, int length, int offset) {
	int s;
	//int gap_count (int * seq, int start, int stop);
	//struct sequence remove_gaps (struct sequence * seq, int start, int size);
	//double pw_pid (struct sequence * seq1,struct sequence * seq2, int start1, int start2, int length);
	int MAX_GAPS=5;
	struct sequence bestseq;
	int minus,plus,ngaps;	
	ngaps = gap_count(orth->asints,start,start+length);
	
	struct sequence refseq=remove_gaps(ref,start,length);	
	if ((ngaps>MAX_GAPS)||(offset>MAX_GAPS)) { 
	  bestseq.length=0; 
	  return bestseq; 
	}
	
	minus=plus=offset+ngaps;
	int begin=start;
	if (minus==0) bestseq = remove_gaps(orth, start, length); 	
	else {
		//cout << ngaps << " gaps and "<<offset<< " gaps in reference" <<endl;
		int ugbp=0;
		int end = ((start+plus+length)<=orth->length) ? (start+plus) : (orth->length-length) ;	
		while ((begin>0)&&(ugbp<minus)) {
			begin--; 
			if (orth->is_base(begin)) ugbp++;
		}		
			//find the best ungapped local alignment									
		//cout << "looking from " << begin << " to " << end <<" to match "<<refseq.bases<<endl;	
		double bestscr=0.0;
		int bestq,q;
		for (q=begin;q<=end;q++) if (orth->is_base(q)) {
			struct sequence temp = remove_gaps(orth,q,length);
			if (temp.length!=length) return temp;		
			double scr=pw_pid(&refseq,&temp,0,0,length);
			if (scr>bestscr) { bestq=q; bestscr=scr; }
			else if ((scr == bestscr)&& (abs(bestq-start)>abs(q-start))) { bestq=q; bestscr=scr; }
			//cout << q << " " << scr << " " <<temp.bases<<endl;	
		}
			//store the best ungapped wmer
		bestseq = remove_gaps(orth, bestq, length); 
	}
	//delete refseq;
	refseq.erase();
	return bestseq;	
}

double pw_pid (struct sequence * seq1,struct sequence * seq2, int start1, int start2, int length) {
	int p; int id=0;
	for (p=0;p<length;p++) 
		if ((seq1->is_base(p+start1)) && (seq1->asints[p+start1]==seq2->asints[p+start2])) id++;
	//cout << id << " " <<endl;	
	return ((double) id / (double) length);
}

int gap_count (int * seq, int start, int stop) {
	int pos;
	int n=0;
	for (pos=start;pos<stop;pos++) if ((seq[pos]==4)||(seq[pos]==-1)) n++;
	return n;
}

double pid (struct sequence * seq, int ns, int start, int size) {
	int s,p; float id=0;
	for (p=start;p<(size+start);p++) {
		bool same=1;
		for (s=1;s<ns;s++) if (seq[s].asints[p] != seq[0].asints[p]) same=0;
		if (same) id+=1.0;
	}
	return id/((double) size);
}

double pug (struct sequence * seq, int ns, int start, int size) {
	int s,p; float ug=0;
	for (p=start;p<(size+start);p++) {
		bool no_gap=1;
		for (s=0;s<ns;s++) if ((seq[s].asints[p] == 4) || (seq[s].asints[p] == -1)) no_gap=0;
		if (no_gap) ug+=1.0;
	}
	return ug/((double) size);
}


struct sequence remove_gaps (struct sequence * seq, int start, int size) {
	struct sequence ug;
	ug.name = seq->name + "_sub";
	int pos=start;
	while (pos<seq->length) {
		if (seq->is_base(pos)) ug.bases += seq->bases[pos];
		if(ug.bases.size()>=size) pos=seq->length;
		//cout << ug.bases <<endl;
		pos++;
	}
	ug.integerize();
	//ug.print(0,size-1);
	return ug;
}		
	
double skip_gap_log_LR (struct motif * M, struct motif * B, struct sequence * seq, int start) {
	int pos=0; int skipped=0; 
	int * C= new int [4]; C[0]=3; C[1]=2; C[2]=1; C[3]=0;
	double llr=0.0; double rcllr=0.0;
	while (pos<M->length) {
		int i=start+pos+skipped;
		if (i >= seq->length) return -1000.0;		
		else if (seq->is_base(i)) {
			llr += log( M->f[pos][seq->asints[i]] / B->f[pos][seq->asints[i]] );
			rcllr += log( M->rc->f[pos][seq->asints[i]] / B->f[pos][C[seq->asints[i]]] );
			pos++;
		}			
		else skipped++;
	}
	delete C;
	if (rcllr>llr) return rcllr;
	else return llr;	
}

double keep_gap_log_LR (struct motif * M, struct motif * B, struct sequence * seq, int start) {
	int pos;
	double llr=0.0; double rcllr=0.0;
	for (pos=0;pos<M->length;pos++) {
		int i=start+pos;
		if (i >= seq->length) return -1000.0;		
		else if (seq->is_base(i)) {
			llr += log( M->f[pos][seq->asints[i]] / B->f[pos][seq->asints[i]] );
			rcllr += log( M->rc->f[pos][seq->asints[i]] / B->f[pos][seq->asints[i]] );			
		}			
		//else skipped++;
	}
	
	if (rcllr>llr) return rcllr;
	else return llr;	
}

double div_to_subs_JC (double pid) { // Jukes-cantor parameter as a function of Pid
	if (pid>0.25) {
		double k=-0.75*log((4.0*pid-1.0)/3.0);
		return ((k>MIN_RATE) ? k : MIN_RATE);
	}
	else return -1.0;
}

double p_subtree2 (struct sequence * seq, struct node * n, struct tree * phy, int base, int pos) {
	if (n->leaf) {
		if (seq[n->number].is_base(pos)) return ( (seq[n->number].asints[pos] == base) ? 1.0 : 0.0 );
		else return 1.0;
	}	
	int lb,rb; double Prl=0.0; double Prr=0.0;
	for (lb=0;lb<4;lb++) Prl+= ( p_subtree2(seq,n->left,phy,lb,pos) * phy->psubs[base][lb][n->number][n->left->number] );	
	for (rb=0;rb<4;rb++) Prr+= ( p_subtree2(seq,n->right,phy,rb,pos) * phy->psubs[base][rb][n->number][n->right->number] );
	return Prr*Prl;
}	

double p_subtree_hb (struct sequence * seq, struct node * n, int base, int sp, int mp, struct motif * m) {
	//if (n->leaf) return ( (seq[n->number].asints[sp] == base) ? 1.0 : 0.0 );
	if (n->leaf) { 
		if (seq[n->number].is_base(sp)) return ( (seq[n->number].asints[sp] == base) ? 1.0 : 0.0 );
		else return 1.0;
	}	
	int lb,rb; double Prl=0.0; double Prr=0.0;
	//for (lb=0;lb<4;lb++) Prl+= ( p_subtree_hb(seq,n->left,lb,pos,f,bgf) * PHB(base,lb,n->lbr,f,bgf) );
	//for (rb=0;rb<4;rb++) Prr+= ( p_subtree_hb(seq,n->right,rb,pos,f,bgf) * PHB(base,rb,n->rbr,f,bgf) );
	for (lb=0;lb<4;lb++) Prl+= ( p_subtree_hb(seq,n->left,lb,sp,mp,m) * m->psubs[mp][n->number][n->left->number][base][lb] );
	for (rb=0;rb<4;rb++) Prr+= ( p_subtree_hb(seq,n->right,rb,sp,mp,m) * m->psubs[mp][n->number][n->right->number][base][rb] );
	return Prr*Prl;
}

double p_subtree_jc (struct sequence * seq, struct node * n,  int base, int pos, double k) {
	//if (n->leaf) return ( (seq[n->number].asints[pos] == base) ? 1.0 : 0.0 );
	if (n->leaf) { 
		if (seq[n->number].is_base(pos)) return ( (seq[n->number].asints[pos] == base) ? 1.0 : 0.0 );
		else return 1.0;
	}	
	int lb,rb; double Prl=0.0; double Prr=0.0;
	for (lb=0;lb<4;lb++) Prl+= ( p_subtree_jc(seq,n->left,lb,pos,k) * PJC(base,lb,k*n->lbr) );
	for (rb=0;rb<4;rb++) Prr+= ( p_subtree_jc(seq,n->right,rb,pos,k) * PJC(base,rb,k*n->rbr) );
	return Prr*Prl;
}

double logLR2 (struct motif * m, struct sequence * seq, struct tree * phy, double * bgf, int pos, int e_mod) {
	double pm,prc,pbg,lhbg,lhm,lhrc;
	int p,b;	
	//double bcond[4]; double mcond[4]; double rcond[4];
	double * bcond = new double [4];
	double * mcond = new double [4];
	double * rcond = new double [4];
	
	pm=prc=pbg=1.0;
	for (p=0;p<m->length;p++) {
		for (b=0;b<4;b++) bcond[b]= p_subtree2(seq,phy->root,phy,b,p+pos);
		if (e_mod==0) { for (b=0;b<4;b++) mcond[b]=rcond[b]=bcond[b]; }
		else if (e_mod==1) { for (b=0;b<4;b++) mcond[b]=rcond[b]=p_subtree_jc(seq,phy->root,b,p+pos,0.5); }
		else if (e_mod==2) { for (b=0;b<4;b++) mcond[b]=p_subtree_hb(seq,phy->root,b,p+pos,p,m); }
		lhbg=lhm=lhrc=0.0;
		for (b=0;b<4;b++) {
			lhbg+= bgf[b] * bcond[b];
			lhm+= m->f[p][b] * mcond[b];
			
			//if (e_mod!=0) 
			//cout << m->f[p][b]<<" x "<<  mcond[b]<<"\t";
		}
		pm*=lhm;
		pbg*=lhbg;
		//if (e_mod!=0) { 
		//cout <<"mod"<<e_mod <<" pos"<<p+pos<<" "; for (b=0;b<4;b++) cout <<bcond[b] << " "; cout <<endl;
		//}
		if (BOTH_STRANDS) { 
			if (e_mod==2) for (b=0;b<4;b++) rcond[b]=p_subtree_hb(seq,phy->root,b,p+pos,p,m->rc);
			for (b=0;b<4;b++) lhrc += m->rc->f[p][b] * rcond[b];
			prc *= lhrc;
		}	
	}
	//if (PRINT_TO_SCREEN) {
	//	cout << e_mod<<": pm="<<log(pm) << " (prc="<<log(prc)<<") pbg="<< log(pbg) <<"LR= " ;
	//	if ((BOTH_STRANDS)&&(prc>pm)) cout << log(prc)-log(pbg);
	//	else cout << log(pm)-log(pbg);
	//	cout <<endl;
	//}	
	
	delete rcond; delete mcond; delete bcond;
	if ((BOTH_STRANDS)&&(prc>pm)) return (log(prc)-log(pbg));	
	return (log(pm)-log(pbg));
}

double PJC (int a, int b, double alpha) { //jukes-cantor probability of a to b in alpha//
	//cout << a<<"->"<<b<<endl;
	if (a!=b) return (0.25-0.25*exp((-4.0/3.0)*alpha));
	else return (0.25+0.75*exp((-4.0/3.0)*alpha));

}

sequence::sequence() {
	length=0;
}

//double sequence::Pr(int start, int stop, double * f) {
//	int p; double pr=1.0;
//	for (p=start;p<stop;p++) if ((asints[p]!=4)&&(asints[p]!=-1)) pr *= f[p];
//	return pr;
//}
void sequence::erase() {
  bases.erase();
  delete asints;
}
bool sequence::print(int start, int stop) {
	if (start>=stop) return 0;
	if (stop>length) return 0;
	int p,b; char B[5]; B[0]='A'; B[1]='C'; B[2]='G'; B[3]='T'; B[4]='-';
	cout << name <<"("<<start<<"-"<<stop-1<<")\t";
	for (p=start;p<stop;p++) {
		if (asints[p]==-1) cout <<"N";
		else cout <<B[asints[p]];
	}
	cout << endl;
	return 1;
}
bool sequence::integerize() {
	length=bases.size();
	int p,b;
	//char B[4];
	asints = new int [length];
	//B[0]='A'; B[1]='C'; B[2]='G'; B[3]='T';
	//for (p=0;p<length;p++) for (b=0;b<4;b++) if (B[b] == bases[p]) asints[p]=b;
	//for (p=0;p<length;p++) cout <<seq[p]<<endl;
	for (p=0;p<length;p++) asints[p] = (
				(bases[p]=='A') ? 0 : (
					(bases[p]=='C') ? 1 : (
						(bases[p]=='G') ? 2 : (
							(bases[p]=='T') ? 3 : (
								(bases[p]=='-') ? 4 : -1
									)
								)
							)
						)
					);

	return 1;
}

bool sequence::is_base(int pos) {
	if ((asints[pos]==4)||(asints[pos]==-1)) return 0;
	else return 1;
}	
	

node::node() {
	leaf=0;
	root=0;
	//seq_post = new double [4];
	//int b;
	//for (b=0;b<4;b++) seq_post[b]=0.0;
	base=4;
	lbr=abr=rbr=INITIAL_BRANCH;

}

bool node::print() {
	char B[4]; B[0]='A'; B[1]='C'; B[2]='G'; B[3]='T';
	int b;
	if (root) cout << "\troot"<<endl;
	else cout <<"\tn"<<anc->number<<" "<<endl;
	cout <<"\t|"<<endl;
	if (leaf) cout << "\tleaf" <<endl;
	else cout <<"n"<< left->number << "\t-\tn" <<right->number<<endl;
	return 1;
}

tree::tree() {
	Nnodes=0;
}


bool tree::print(struct sequence * seq) {
	cout << tree_string(root, seq) << endl;
	//int n=0;
// 	for (n=0;n<Nnodes;n++) {
// 	//while (!done) {
// 		cout << "node " <<n<<endl;
// 		nodelist[n].print();
// 	}
	return 1;
}

string tree::tree_string( struct node * n, struct sequence * seq) {
	if (n->leaf) return seq[n->number].name;
	else return "("+tree_string(n->left,seq) +","+ tree_string(n->right,seq)+")";
}

bool tree::compute_PJC() {
	int a,b,n;
	for (a=0;a<4;a++) for (b=0;b<4;b++) psubs[a][b] = new double * [Nnodes];
	for (n=0;n<Nnodes;n++) for (a=0;a<4;a++) for (b=0;b<4;b++) psubs[a][b][n] = new double [Nnodes];
	for (n=0;n<Nnodes;n++) if (!nodelist[n].root) for (a=0;a<4;a++) for (b=0;b<4;b++) {
		if (a != b) psubs[a][b][nodelist[n].anc->number][n] = 0.25 - 0.25*exp(-(4.0/3.0)*nodelist[n].abr);
		else psubs[a][b][nodelist[n].anc->number][n] = 0.25 + 0.75*exp(-(4.0/3.0)*nodelist[n].abr);
	}
	//for (n=0;n<Nnodes;n++) if (!nodelist[n].root) for (a=0;a<4;a++) for (b=0;b<4;b++) {
		//cout << a<<"->"<<b<< " "<<nodelist[n].anc->number<<"<"<<nodelist[n].abr<<">"<< n<<":";
	//	cout<<psubs[a][b][nodelist[n].anc->number][n]<<endl;
	//}	
	return 1;
}
bool tree::compute_PHKY(double * f, double kap) {
	int a,b,n;
	void matrix_exponential (double ** a, double ** b, int n);
	for (a=0;a<4;a++) for (b=0;b<4;b++) psubs[a][b] = new double * [Nnodes];
	for (n=0;n<Nnodes;n++) for (a=0;a<4;a++) for (b=0;b<4;b++) psubs[a][b][n] = new double [Nnodes];
	for (n=0;n<Nnodes;n++) if (!nodelist[n].root) {
		double ** R = new double * [4];
		double ** P = new double * [4];		
		for (a=0;a<4;a++) {
			double tot=0.0; double idrate=0.0;
			R[a] = new double [4];
			P[a] = new double [4];
			for (b=0;b<4;b++) if (a != b) {
				R[a][b] = f[b];
				//R[a][b]=1.0;
				if (((a==0)&&(b==2)) || ((a==2)&&(b==0)) || ((a==1)&&(b==3)) || ((a==3)&&(b==1)))
					R[a][b] *= kap;
				tot+=R[a][b];
			}	
			for (b=0;b<4;b++) if (a!=b) {
				R[a][b] *= (nodelist[n].abr/tot);
				idrate -= R[a][b];
			}	
			R[a][a] = idrate;			
		}
		//for (a=0;a<4;a++) { for (b=0;b<4;b++) cout << R[a][b]<< " "; cout <<endl; }
		//cout << endl;
		matrix_exponential(R, P, 4);
		for (a=0;a<4;a++) for (b=0;b<4;b++) psubs[a][b][nodelist[n].anc->number][n] = P[a][b];
		//for (a=0;a<4;a++) { for (b=0;b<4;b++) cout << psubs[a][b][nodelist[n].anc->number][n]<< " "; cout <<endl; }
		//cout << endl;
		delete R;
		delete P;
	}	
		
	return 1;	
}

bool tree::initialize_pairwise(double dis) {
	int n;
	Nnodes=3;
	nodelist=new struct node [3];
	nodelist[2].number=2;
	root= &nodelist[2];
	nodelist[2].root=1;
	nodelist[2].left = &nodelist[0];
	nodelist[2].right= &nodelist[1];
	nodelist[2].lbr=nodelist[2].rbr=0.5*dis;
	for (n=0;n<2;n++) {
		nodelist[n].number=n;
		nodelist[n].leaf=1;
		nodelist[n].anc= &nodelist[2];
		nodelist[n].abr=0.5*dis;
	}
	return 1;
}
bool tree::read_from_file(string treefile) {
	ifstream treedata(treefile.c_str());
	if (!treedata) return 0;
	char comm;
	int thisnum,num;
	treedata>>comm;
	if (comm == 't') {
		treedata>>Nnodes;
		//cout << Nnodes <<endl;
		nodelist=new struct node [Nnodes];
	}
	int i=0;
	while (i < Nnodes) {
		treedata>>comm; //print();
		//if (!comm) return 0;
		if ((comm == 'r')||(comm == 'n')) {
			treedata>> thisnum;
			//cout <<thisnum<<endl;
			//cout << i << endl;
			nodelist[thisnum].number=thisnum;
			if (comm=='r') {
				root = &nodelist[thisnum];
				nodelist[thisnum].root=1;
			}
			treedata>>comm;
			if (comm == 'l') {
				treedata>>num;
				nodelist[thisnum].left = &nodelist[num];
				nodelist[num].leaf=1;
				nodelist[num].number=num;
				nodelist[num].anc = &nodelist[thisnum];
				//nodelist[num].lbr=nodelist[num].rbr=-1.0;
				i++;
			}
			else if (comm == 'n') {
				treedata>>num;
				nodelist[thisnum].left = &nodelist[num];
				nodelist[num].anc= &nodelist[thisnum];
			}
			treedata>>comm;
			if (comm == 'l') {
				treedata>>num;
				nodelist[thisnum].right= &nodelist[num];
				nodelist[num].leaf=1;
				nodelist[num].number=num;
				nodelist[num].anc= &nodelist[thisnum];
				i++;
			}
			else if (comm == 'n') {
				treedata>>num;
				nodelist[thisnum].right= &nodelist[num];
				nodelist[num].anc= &nodelist[thisnum];
			}
			i++;
		}
		//else {cout << i ; print();}
	//print();
	}
	treedata>>comm;
	if (comm == 'b') {
		for (i=0;i<Nnodes;i++) if (!nodelist[i].leaf) {
			treedata>>num; if (num!=i) return 0;
			treedata>>nodelist[num].lbr;
			treedata>>nodelist[num].rbr;
			nodelist[num].left->abr=nodelist[num].lbr;
			nodelist[num].right->abr=nodelist[num].rbr;
		}
	}
	//cout << "read "<<endl;
	return 1;

}

alignment::alignment() {
	length=0;
	Nseqs=0;
	gapped=0;
	bad=0;
	seq = new struct sequence [MAX_ALIGN_SIZE];
	bgfreq = new double [4];
}
bool alignment::read_bg_freq(string freqfile) {
	char B[4]; B[0]='A'; B[1]='C'; B[2]='G'; B[3]='T';
	ifstream data(freqfile.c_str());
	if (!data) return 0;
	//else if (PRINT_TO_SCREEN) cerr << "reading "<<freqfile<<" ";
	int b;
	for (b=0;b<4;b++) {
		char base;
		data>>base;
		if (B[b] == base) { 
		  data>>bgfreq[b];  
		  //if (PRINT_TO_SCREEN) cerr<<B[b]<< ": "<<bgfreq[b]<<" "; 
		}
		else { cerr << "coudn't parse " << freqfile <<endl; return 0; }
	}
	// if (PRINT_TO_SCREEN) cerr<<endl;
	char comm;
	data>>comm;
	if (comm == 'f')  data >> FUDGE; 
	else if (comm == 'k') data >> kappa;
	//	if (PRINT_TO_SCREEN) cerr << "ki/ks=" << 1.0/FUDGE<<endl; 
	
	return 1;

}

bool alignment::print() {
	int b; char B[4]; B[0]='A'; B[1]='C'; B[2]='G'; B[3]='T';
	cout << name << " " <<length<<" bp (" <<(length-gapped-bad) << " ungapped)" <<endl;
	for (b=0;b<4;b++) cout << B[b] <<": "<<bgfreq[b]<<" ";
	cout <<endl;
	return 1;
}

void alignment::integerize() {
	int p, b,s;
	asints = new long int [length];
	long int seq_to_num (int * seq, int start, int ml);
	//long int gapped_seq_to_num (int * seq, int start, int ml);
	for (p=0;p<length;p++) {
		bool gapped=0;
		for (s=0;s<Nseqs;s++) if (!seq[s].is_base(p)) gapped=1;
		if (gapped) { asints[p] =-1; }
		else { 
			int * col = new int [Nseqs];
			for (s=0;s<Nseqs;s++) col[s] = seq[s].asints[p];
		//asints[p] = gapped_seq_to_num(col,0,Nseqs);			
			asints[p] = seq_to_num(col,0,Nseqs);
		}
		//else asints[p]=-1;
		//cout << p << " " << asints[p] << " ";
		//for (s=0;s<Nseqs;s++) cout << " "<<seq[s].asints[p];
		//cout << endl;
	}
}

bool alignment::compute_nextgap() {
	char B[4]; B[0]='A'; B[1]='C'; B[2]='G'; B[3]='T';
	int p,b,s;
	int lastgap=length;
	int lastbad=length;
	leftend=0; rightend=-1;
	//double total=0.0;
	nextgap= new int [length];
	nextbad= new int [length];
	double totalbases=0.0;
	gapped=0;
	bad=0;
	for (b=0;b<4;b++) bgfreq[b]=0.0;
	for (p=length-1;p>=0;p--) {
		for (s=0;s<Nseqs;s++) {
			if ((seq[s].bases[p] != 'A') && (seq[s].bases[p] != 'C')
				&& (seq[s].bases[p] != 'G') && (seq[s].bases[p] != 'T')) {
				if (seq[s].bases[p] == '-') lastgap=p; //seq[s].asints[p]=-1;
				else lastbad=p;
			}
			else for (b=0;b<4;b++) if (seq[s].bases[p] == B[b]) {
				bgfreq[b]+=1.0;
				totalbases+=1.0;
			}

		}
		nextgap[p]=lastgap;
		nextbad[p]=lastbad;
		if (p==lastgap) gapped ++;
		else leftend=p;
		if (p==lastbad) bad ++;
		if (rightend==-1) if (gapped != (length-p)) rightend=p;
	}
	for (b=0;b<4;b++) bgfreq[b] /= totalbases;
	//for (p=0;p<length;p++) {
		//for (s=0;s<Nseqs;s++) cout << seq[s].bases[p];
		// cout << p <<" "<< nextgap[p] <<" "<<endl;
	//}
	return 1;
}

float alignment::div() {
	float id=0.0;
	//float ungapped=0.0;
	int p,s1,s2;
	for (p=0;p<length;p++) {
		int thisposdiff=0;
		bool thisgapped=0;
		for (s1=0;s1<Nseqs;s1++) if (seq[s1].asints[p]==4) thisgapped=1;
		if (!thisgapped) {
			//for (s1=0;s1<Nseqs;s1++)
			for (s2=1;s2<Nseqs;s2++)
				if (seq[0].asints[p] != seq[s2].asints[p]) thisposdiff=1;
			//ungapped+=1.0;
			if (!thisposdiff) id+=1.0;
		}
	}
	id /= ((float) (length-gapped-bad));
	return (1.0-id);
}

bool alignment::read_tree(string filename) {
	int p,n,b;
	//if (PRINT_TO_SCREEN) cerr << "reading tree from " <<filename<<"...";
	//phyl= new struct tree [length];
	if (!phyl.read_from_file(filename)) return 0;
	treelength=0.0;
	for (n=0;n<phyl.Nnodes;n++) if (!phyl.nodelist[n].root)
		treelength += phyl.nodelist[n].abr;
	//cout << global_tree.root->lbr<<endl;
	//double eqf[4]; eqf[0]=eqf[1]=eqf[2]=eqf[3]=0.25;
	//double loglh=0.0;
	//for (p=0;p<length;p++) if (nextgap[p]!=p) {
	//for (p=0;p<1;p++) {
		//if (!phyl[p].read_from_file(filename)) return 0;
		//phyl[p].print();
		//ALI[0].phyl[p].print();
		//for (n=0;n<phyl[p].Nnodes;n++) {
		//	if (phyl[p].nodelist[n].leaf) {
		//		phyl[p].nodelist[n].base=seq[n].asints[p];
				//for (b=0;b<4;b++) ALI[0].phyl[p].nodelist[n].seq_post[b]=0.0;
 		//		phyl[p].nodelist[n].seq_post[phyl[p].nodelist[n].base]=1.0;
		//	}
		//	else if ((global_tree.root->lbr==INITIAL_BRANCH)&&(global_tree.root->rbr==INITIAL_BRANCH)) {
		//			for (b=0;b<4;b++) phyl[p].nodelist[n].seq_post[b]=bgf[b];
		//	}
		//	else if (EVOL_MODEL==1) compute_post(phyl[p].nodelist + n, bgf);
		//	else if (EVOL_MODEL==0) compute_post_JC(phyl[p].nodelist + n, bgf);
		//}
	//	loglh+=log(compute_post_JC(phyl[p].root,bgf));
		//compute_post(phyl[p].root->left);
		//compute_post(phyl[p].root->right);
		//loglh+=compute_post(phyl[p].root, bgf);
	//}
	//cout << "initial bg log liklihood " << loglh <<endl;
	//if (EVOL_MODEL==0) {
	//if ((global_tree.root->lbr==INITIAL_BRANCH)&&(global_tree.root->rbr==INITIAL_BRANCH))
	//	optimize_branch_lengths(EPSILON,MAX_ITER);
	if ((phyl.root->lbr!=INITIAL_BRANCH)&&(phyl.root->rbr!=INITIAL_BRANCH))
	  //	if (PRINT_TO_SCREEN) cerr << "branch lengths found" <<endl;
	  //else if (PRINT_TO_SCREEN) cerr << "branch lengths not found" <<endl;
	//}
	//loglh=0.0;
	//for (p=0;p<length;p++) if (nextgap[p]!=p)
	//	loglh+=log(compute_post_JC(phyl[p].root,bgf));
	//cout << "maximum bg liklihood "<<loglh<<endl;
	return 1;
}

bool alignment::read_fasta(string filename) {
	ifstream data(filename.c_str());
	if (!data) return 0;
	//else if (PRINT_TO_SCREEN) cerr << "reading "<<filename<<"...";
	//string seqdata;
	Nseqs=-1;
	bool foundmax=0;
	int p; bool ext=1;
	for (p=filename.size();p>=0;p--) { //get rid of the last extension
		if (ext) { if (filename[p]=='.') ext=0; }
		else name=filename[p] + name;
		//cout << filename[p] <<" "<<ext<<" "<<name<<endl;
	}
	if (ext) name=filename; // if there was no extension
	//name[filename.size()]='\0';
	//name=filename;
	while ((!data.fail())&&(!foundmax)) {
		string line;
		getline(data,line);
		int linelength=line.size();
		//cout << linelength << endl;
		if ( line[0] == '>' ) {
			Nseqs++;
			if (Nseqs<MAX_ALIGN_SIZE) { for (p=1;p<linelength;p++) seq[Nseqs].name += line[p]; }
			else foundmax=1;
			//if (Nseqs>0) if (seq[Nseqs].bases.size()>length)
			//		length=seq[Nseqs].bases.size();
				//seq[Nseqs].name = new char [linelength];
				//for (p=1;p<linelength;p++) {
				//	seq[Nseqs].name+=line[p];
					//cout << p << " " << line[p]<< endl;
				//}
				//seq[Nseqs].name[linelength-1]='\0';
				//seq[Nseqs].name=line.c_str();
				//cout << Nseqs<< " "<<seq[Nseqs].name << endl;
			//}
		}
		else if ((linelength>0)&&(!foundmax)) {
			seq[Nseqs].bases += line;
			//cout << line << endl;
			//for (p=0;p<linelength;p++) if (line[p] == '\n') cout << p <<endl;
		}
	}
	int s;
	if (Nseqs<MAX_ALIGN_SIZE) Nseqs++;
	for (s=0;s<Nseqs;s++) {
		seq[s].integerize();
		if (length < seq[s].length) length=seq[s].length;
	}
	compute_nextgap();
	integerize();
	//P_BG_COMPUTED=0;
	//pbg=new double [length];
	kappa=2.0;
	double div_to_subs_JC (double pid);
	if (Nseqs==2) { // initialise the pairwise tree

		treelength=div_to_subs_JC(1.0-div()); //JC distance
		phyl.initialize_pairwise(treelength);
		//double d=div();
		//double r=(1.0-d)/d;
		//double k=-0.75*log((1.0-(r/3.0))/(1.0-r));
		//phyl.initialize_pairwise(k);
		//int b;
		//phyl= new struct tree [length];
		//for (p=0;p<length;p++) if (nextgap[p]!=p) {
		//	phyl[p].Nnodes=3;
		//	phyl[p].nodelist= new struct node [3];
		//	phyl[p].nodelist[2].number=2;
		//	phyl[p].root= &phyl[p].nodelist[2];
		//	phyl[p].nodelist[2].root=1;
		//	phyl[p].nodelist[2].left = &phyl[p].nodelist[0];
		//	phyl[p].nodelist[2].right= &phyl[p].nodelist[1];
		//	for (b=0;b<4;b++) phyl[p].nodelist[2].seq_post[b]=bgfreq[b];
		//	for (s=0;s<Nseqs;s++) {
		//		phyl[p].nodelist[s].number=s;
		//		phyl[p].nodelist[s].leaf=1;
		//		phyl[p].nodelist[s].anc= &phyl[p].nodelist[2];
		//		phyl[p].nodelist[s].base=seq[s].asints[p];
		//		phyl[p].nodelist[s].seq_post[phyl[p].nodelist[s].base]=1.0;
		//	}
		//}
		//phyl.Nnodes=3;
		//phyl.nodelist=new struct node [3];
		//phyl.nodelist[2].number=2;
		//phyl.root= &phyl.nodelist[2];
		//phyl.nodelist[2].root=1;
		//phyl.nodelist[2].left = &phyl.nodelist[0];
		//phyl.nodelist[2].right= &phyl.nodelist[1];
		//for (s=0;s<Nseqs;s++) {
		//	phyl.nodelist[s].number=s;
		//	phyl.nodelist[s].leaf=1;
		//	phyl.nodelist[s].anc= &phyl.nodelist[2];
		//}
	}
	//if (PRINT_TO_SCREEN) { cerr << "done" <<" %id="<<1.0-div()<<" k="<<treelength<<endl; }
	if (Nseqs>1) return 1;
	else return 0;
}
motif::motif() {
	//name = new char [MAX_NAME_LENGTH];
	//f = new double * [4];
	pri=MOTIF_PRIOR;
}
double motif::Pr(int s[],int start) {
	int pos; double p=1.0;
	for (pos=0;pos<length;pos++) {
		int b=s[start+pos];		
		if ((b!=-1)&&(b!=4)) p *= f[pos][b];
	}	
	return p;
}

struct motif motif::make_rc() {
	struct motif rc;
	int b, pos;
	int revcom[4]; revcom[0]=3; revcom[1]=2; revcom[2]=1; revcom[3]=0;
	rc.length=length;
	rc.pri=pri;
	rc.K = new double [length];
	rc.f = new double * [length];
	for (pos=0;pos<length;pos++) rc.f[pos] = new double [5];
	for (pos=0;pos<length;pos++) {
		for (b=0;b<4;b++) rc.f[pos][b]=f[length-pos-1][revcom[b]];
		rc.K[pos]=K[length-pos-1];
		rc.f[pos][4]=f[length-pos-1][4];
	}
	int a;	
	return rc;
}
bool motif::initialize_neutral_motif(int w, double pseudocounts) {
	int b, pos;
	length=w;
	//cout << "initialising matrix of size "<< length <<endl;
	f = new double * [length];
	K = new double [length];
	for (pos=0;pos<length;pos++) f[pos] = new double [5];
	for (pos=0;pos<length;pos++) {
		for (b=0;b<4;b++) f[pos][b]=0.25*pseudocounts;
		K[pos]=1.0;
		f[pos][4]=1.0;
	}
	return 1;
}


bool motif::print() {
	char B[4]; B[0]='A'; B[1]='C'; B[2]='G'; B[3]='T';
	int b,pos;
	cout << "pos\t";
	for (b=0;b<4;b++) cout << B[b] <<"\t";
	cout << "Kp"<<"\t"<<"Khb"<<endl;
	for (pos=0;pos<length;pos++) {
		cout << pos << "\t";
		for (b=0;b<4;b++) cout << 0.01*((int) (100.0*f[pos][b])) <<"\t";
		cout << K[pos]<< "\t"<<halpern_bruno(pos)<<endl;

	}
	return 1;
}
double motif::halpern_bruno(int pos) {
	double hbrate=0.0;
	int b1,b2;
	for (b1=0;b1<4;b1++) for (b2=0;b2<4;b2++) if (b2!=b1) {
		if (f[pos][b1]==f[pos][b2]) hbrate+=f[pos][b1];
		else hbrate+=(f[pos][b1]*log(f[pos][b2]/f[pos][b1])/(1-(f[pos][b1]/f[pos][b2])));
	}
	return (hbrate/3.0);
}
bool motif::read_from_file(string filename) {
	ifstream data(filename.c_str());
	if (!data) { cout << filename<<" not found" << endl; return 0; }
	length=-1;
	while (!data.fail()) {
		string line;
		getline(data,line);
		//cout << line << endl;
		length++;
	}
	if (PRINT_TO_SCREEN) cout << "motif of length " << length << " found" << endl;
	int b,pos;
	f = new double * [length];
	for (pos=0;pos<length;pos++) f[pos] = new double [5];
	//int pos=0;
	ifstream datatoo(filename.c_str());
	for (pos=0;pos<length;pos++) {		
		for (b=0;b<4;b++) datatoo>>f[pos][b];			
		f[pos][4]=1.0;
	}

	K=new double [length];
	for (pos=0;pos<length;pos++) K[pos]=1.0;
	name=filename;
	return 1;
}
