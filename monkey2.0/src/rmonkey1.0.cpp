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
int MAX_GAPS=5;
const double MOTIF_PRIOR=0.001;
double thr=0.01;

const double MIN_RATE=0.001;
const int MAX_HITS=10000;
const int EVOL_MODEL=2; // 0= JC bg rate,  1= JC 1/2 bg rate 2=HB selection
const int N_MODELS=100;
double FUDGE=1.0; //correction because the HB model is based on the synonymous rate
double N_BINS=200.0;
double MODEL_ONE_RATE=0.5;
bool HEURISTIC=1;


struct sequence { // sequence object - 2 copies of the sequence, a string, and as ints
	int length;
	string name;
	string bases;
	int * asints;

	sequence();
	//~sequence();
	bool integerize();
	void erase();
	bool print(int start, int stop);
	bool is_base(int pos);
	int next_base(int pos);
	int ungapped_pos(int pos);
	//double Pr(int start, int stop, double * f);
};

struct node {
	bool leaf; bool root;
	struct node * anc;
	struct node * left;
	struct node * right;
	double abr,lbr,rbr;
	//double * seq_post;
	int base; int number; int state;

	node();
	bool print();
	struct node * sister();
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
	bool copy_topology(struct tree * t);
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


struct hit {
	int * apos;
	int * aend;
	int * ugpos;
	int refs; int leftend; int rightend;
	struct alignment ali;
	double sspval; double evolpval;
	double strand;

	hit();
	bool print();

};

double p_subtree2 (struct sequence * seq, struct node * n, struct tree * phy, int base, int pos);
double p_subtree_hb (struct sequence * seq, struct node * n, int base, int sp, int mp, struct motif * m);
double p_subtree_jc (struct sequence * seq, struct node * n,  int base, int pos, double k);

double PJC (int a, int b, double alpha);


int gap_count (int * seq, int start, int stop);

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

	bool PVALS=1 ; bool LLRS=0; bool LIST=0; bool EVALS=0;
	bool ANY_SP=0;  int C_MEM=0;
	double KAPPA = 2.0; bool DAN=1; bool BLOCK_MODE=0; int WIDTH = 10;

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
			//else if (argv[i][1] == 'e') { EVALS=1; }
			else if (argv[i][1] == 'q') { PRINT_TO_SCREEN=0; }
			else if (argv[i][1] == 'o') { outfile=argv[i+1]; i++; }
			else if (argv[i][1] == 'D') { cout << "hey dan"<< endl; DAN=1; }
			else if (argv[i][1] == 'R') { MODEL_ONE_RATE = atof(argv[i+1]); i++; }
			//else if (argv[i][1] == 'B') { BLOCK_MODE = 1; }
			else if (argv[i][1] == 'a') { ANY_SP=1; }
			else if (argv[i][1] == 'M') { C_MEM=atoi(argv[i+1]); i++; }
			//else if (argv[i][1] == 'W') { WIDTH=atoi(argv[i+1]); i++; }
			else { cout << argv[i] << " is not a recognized option"<<endl; help(); return 0; }

		}
		else if (mfile[0] == '?') { mfile = argv[i]; }
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

	//if (BLOCK_MODE) {
	//	if (mmod!=1) { cerr << " use -m JC for conserved blocks"<<endl; return 0; }
	//	thr=-1.0;
	//	WIDTH = atoi(mfile.c_str());
	//	M[1].initialize_neutral_motif(WIDTH,0.0);
	//	if (PRINT_TO_SCREEN) cerr << "testing blocks of size "<<WIDTH <<" at "<<MODEL_ONE_RATE<<endl;
	//}
	//else 
	if (!M[1].read_from_file(mfile)) {
		cerr <<"matrix file not read\n";
		return 0;
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
	if ((ALI[0].Nseqs>2)&&(PRINT_TO_SCREEN))
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
	//if (BLOCK_MODE) for (p=0;p<M[1].length;p++) for (b=0;b<4;b++) M[1].f[p][b]=ALI[0].bgfreq[b];
	for (p=0;p<M[1].length;p++) for (b=0;b<4;b++) M[0].f[p][b]=ALI[0].bgfreq[b];

	//double logLR (struct motif * mod, struct alignment * ali, int pos);
	int remove_gaps (struct sequence * seq, struct sequence * ug, int start, int size);
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
		if (C_MEM<2) {
		if (PRINT_TO_SCREEN) cerr << "working on model " << 5 << "..";
		int r; 
		for (r=0;r<ALI[0].Nseqs;r++) {
			if ((r==ref_num) || (ANY_SP)) {
				compute_pssm_pvals(&M[1], &M[0], &ALI[0].phyl,ALI[0].Nseqs,5+r,r);		
				if (PRINT_TO_SCREEN) cerr <<".";
			}	
		}	
		if (PRINT_TO_SCREEN) cerr<<"done"<<endl;
		}
		if (C_MEM==0) {
			int r; 
		if (PRINT_TO_SCREEN) cerr << "working on model "<<5+ALI[0].Nseqs << "..";
		 
		for (r=0;r<ALI[0].Nseqs;r++) {
			if ((r==ref_num) || (ANY_SP)) {
				compute_pssm_pvals(&M[1], &M[0],&ALI[0].phyl,ALI[0].Nseqs,5+ALI[0].Nseqs+r,r);		
				if (PRINT_TO_SCREEN) cerr <<".";
			}	
		}
		//compute_pssm_pvals(&M[1], &M[0], &ALI[0].phyl,ALI[0].Nseqs,6,ref_num);
		if (PRINT_TO_SCREEN) cerr<<"done"<<endl;
		}


	}
	//return 1;
	if (outfile[0] == '?') {
		if (LIST) outfile=fastafile+".occ";
		else outfile=ALI[0].name+".occ";
	}
	else if (PRINT_TO_SCREEN) cerr << "writing to "<< outfile<<endl;
	//return 1;
	ofstream out(outfile.c_str()); //open the output file and write the headers

	out <<"NAME\tSTRAND\tPOS(0-based)\tUG(1-based)\t";
	out << ALI[0].seq[ref_num].name;
	if (LLRS) out << "\tlog(LR)\t";
	else out << "\tp(LR>lr)\t";
	//if (EVALS) out << "E(LR>lr)\t";
	//if (PVALS) out << "p(LR>lr)\t";
	for (s=0;s<ALI[0].Nseqs;s++) if (s!=ref_num) {
	//	if (DAN) out << "POS\t";
		out <<"POS\t"<< ALI[0].seq[s].name<<"\t";
		if (LLRS) out << "log(LR)\t";
		else out << "p(LR>lr)\t";
	}
	if (LLRS) out << "S^";
	else out << "p(S^>s)";
	out<<"\tT\tp(T<t|HB)\tp(T>t|HKY)\t";//<<endl;
	//return 1;
	//int gap_count (int * seq, int start, int stop);
	double pid (struct sequence * seq, int ns, int start, int size);
	double pug (struct sequence * seq, int ns, int start, int size);
	int optimized_seq(struct sequence * ref, struct sequence * orth, struct sequence * news, int start, int length, int offset);
	double aligned_motif_pval(struct motif * m, long int * ali, int ns, int pos, int e_mod);
	double gapped_aligned_motif_pval(struct motif * m, struct motif * bg, struct alignment * ali, int ns, int pos, int e_mod, int ref);
	double ss_pval(struct motif * m, int * seq, int pos);
	long int col_to_int (struct sequence * seq, int p, int n);
	long int gapped_complem_int (long int i, int l);
	void general_variant(long int mn, int ml, int * seq, int N);
	void generate_variant (long int mn, int ml, int * seq);
	void generate_two_state_tree (int n, struct tree * st);
	
	double logLR2 (struct motif * m, struct sequence * seq, struct tree * phy, double * bgf, int pos, int e_mod);
	double skip_gap_log_LR (struct motif * M, struct motif * B, struct sequence * seq, int start);
	double skip_gap_ss_pval (struct motif * M,  struct sequence * seq, int start);
	double keep_gap_log_LR (struct motif * M, struct motif * B, struct sequence * seq, int start);
	

	out<<endl;
	//return 1;
	
	int total_bp=0;
	int h;
	for (a=0;a<Nalign;a++) {
		double * ssp = new double [ALI[a].length];
		int * bestsp = new int [ALI[a].length];
		struct hit * hitlist;
		int nabove=0;
		if (ANY_SP) {
			//
			for (p=0;p<ALI[a].length-M[1].length;p++) {
				double bp=1.0; int bs;
				for (s=0;s<ALI[a].Nseqs;s++) {
					double pval = skip_gap_ss_pval(&M[1], &ALI[a].seq[s], p);
					if (pval<bp) { bs=s; bp=pval; }
				}
				ssp[p] = bp; bestsp[p] = bs;
				if (ssp[p]<thr) nabove++;
				//if (ssp[p]<thr) { cout << p << " "<<ssp[p]<<" "<<bestsp[p]<<endl; }
			}
			hitlist = new hit [nabove];
			 for (h=0;h<nabove;h++) {
				hitlist[h].apos = new int [ALI[a].Nseqs];
				//0hitlist[h].refs=-1;
			}
		}
		else {
			for (p=0;p<ALI[a].length-M[1].length;p++) {
				//if (ALI[a].seq[ref_num].is_base(p)) {
					if ((p!=0) && (p % 10000==0)&&(PRINT_TO_SCREEN==1)) cerr << ".";
					ssp[p] = skip_gap_ss_pval(&M[1],  &ALI[a].seq[ref_num], p);
					bestsp[p] = ref_num;
					if (ssp[p]<thr) nabove++;
				//if (ssp[p]<thr)
				 //cout << p <<" "<<ssp[p]<< endl;
			//if (ALI[a].seq[ref_num].is_base(p)) ugp++;
			}
			if (ALI[a].length>10000) cerr << endl;
		//cerr << nabove << " ss hits in " << ALI[a].name <<endl;
			hitlist = new hit [nabove];
			for (h=0;h<nabove;h++) {
				hitlist[h].apos = new int [ALI[a].Nseqs];
				//hitlist[h].refs=ref_num;
			}
		}
		bool find_hits ( struct alignment * a, struct motif * m, struct hit * hlist ,
						double * sspvals, int * sssp, double cut, int * hnum, int l, int r);
		int * HN = new int [1]; HN[0]=0;
		find_hits(&ALI[a],&M[1],hitlist,ssp,bestsp,thr,HN,0,ALI[a].length-M[1].length);


		delete [] bestsp; // = new int [ALI[a].length];


		//bool find_top_hit ( double * ssplist, struct hit * tophit, int left, int right, double cut);

		//find_top_hit(ssp,&hitlist[0],0,ALI[a].length-M[1].length,thr);
		//bool make_hit_alignment ( struct alignment * a, struct motif * m, struct hit * h );
		//make_hit_alignment(&ALI[a],&M[1],&hitlist[0]);
		//cout << " top hit:" <<endl; hitlist[0].print();
		 int ugp=0; int s;
		double ss_pval(struct motif * m, int * seq, int pos);
		 double aligned_motif_pval(struct motif * m, long int * ali, int ns, int pos, int e_mod);
		 double TO_statistic (struct motif * m, struct motif * bg, struct tree * t, struct alignment * ali, int pos, int ref, bool strand);
		double Tstat_pval(struct motif * m, long int * ali, int ns, int pos, bool strand, int e_mod);
		double test_leaf_loss(int s, struct alignment * ali, struct tree * t, int pos, struct motif * m, double * bgf, bool strand);
		double test_model_change_below (struct node * n, struct alignment * ali, struct tree * t, int pos,	struct motif * m, double * bgf, bool strand, int mod_above);
		double score_given_state_tree (struct motif * M, struct motif * B, struct sequence * seq, struct tree * phy, struct tree * st1, struct tree * st2, int pos);
		double state_tree_score (struct tree * st1, struct tree * st2, double br, double dr);
		double TO_statistic_given_state_tree (struct motif * m, struct motif * bg, struct alignment * ali, struct tree * t, struct tree * st1, struct tree * st2, int pos, int ref);
		//return 1;

		for (p=0;p<ALI[a].length-M[1].length;p++) {
			if (ALI[a].seq[ref_num].is_base(p)) ugp++;
			//if (ssp[p]<thr)
			for (h=0;h<HN[0];h++) if (hitlist[h].apos[ref_num] == p) {
				bool good=1;
				out << ALI[a].name<<" "<<M[1].name<<"\t"<<hitlist[h].strand
				<<"\t"<<hitlist[h].apos[ref_num]<<"-"<<hitlist[h].aend[ref_num]
				<<"\t"<<ugp<<"\t"<<hitlist[h].ali.seq[ref_num].bases<<"\t";
				if (hitlist[h].ali.seq[ref_num].length < M[1].length) {
					out << "g\t"; good=0;
				}
				else if (LLRS) { out << skip_gap_log_LR (&M[1],&M[0], &hitlist[h].ali.seq[ref_num],0)<<"\t"; }			
				else { out <<ss_pval(&M[1],hitlist[h].ali.seq[ref_num].asints,0)<<"\t"; }
				for (s=0;s<hitlist[h].ali.Nseqs;s++) if (s!=ref_num) {
					out <<hitlist[h].apos[s]<<"-"<<hitlist[h].aend[s]<<"\t"
					<<hitlist[h].ali.seq[s].bases<<"\t";
					if (hitlist[h].ali.seq[s].length < M[1].length) {
						out << "g\t"; good=0;
					}
					else if (LLRS) { out << skip_gap_log_LR (&M[1],&M[0], &hitlist[h].ali.seq[s],0)<<"\t"; }			
					else { out<<ss_pval(&M[1],hitlist[h].ali.seq[s].asints,0)<<"\t"; }
				}
				//	<<ss_pval(&M[1],hitlist[h].ali.seq[s].asints,0)<<"\t";
				//bool ungapped=1;

//					int pos; for (pos=0;pos<M[1].length;pos++) {
//								cout << hitlist[h].ali.asints[pos] << endl;
//
//							}
//							hitlist[h].print(); return 0;

				if (good) {
					if (LLRS) out << logLR2 (&M[1], hitlist[h].ali.seq, &ALI[0].phyl, ALI[0].bgfreq,0,mmod);
					else out <<aligned_motif_pval(&M[1],hitlist[h].ali.asints,hitlist[h].ali.Nseqs,0,2);
				}	
				else { out << "g"; }
				//cout << hitlist[h].apos[ref_num]<<" "<<TO_statistic(&M[1], &M[0], &ALI[a].phyl, &hitlist[h].ali,0,ref_num);
				//cout << " ";
				//cout <<aligned_motif_pval(&M[1],hitlist[h].ali.asints,hitlist[h].ali.Nseqs,0,2);
				//cout << endl;
				//out << "\t"<< hitlist[h].refs;// " was the reference"<<endl;
				if (good) {

				
				if (hitlist[h].strand<0.5) {
					//if (hitlist[h].refs==ref_num) 
					out << "\t"<<TO_statistic(&M[1], &M[0], &ALI[a].phyl, &hitlist[h].ali,0,hitlist[h].refs,0);

					out <<"\t";
					if (C_MEM<2) {
						out<<Tstat_pval(&M[1],hitlist[h].ali.asints,hitlist[h].ali.Nseqs,0,0,5+hitlist[h].refs);
					}
					else { out <<"nc"; }
					out<< "\t";
					if (C_MEM==0) { 
						out<<Tstat_pval(&M[1],hitlist[h].ali.asints,hitlist[h].ali.Nseqs,0,0,5+ALI[0].Nseqs+hitlist[h].refs);
					}
					else { out << "nc"; }
					//else out <<"\tnr\tnr";
		//			for (s=0;s<hitlist[h].ali.Nseqs;s++) {
		//				out << "\t"<< test_leaf_loss(s,&hitlist[h].ali,&ALI[a].phyl,0,&M[1],M[0].f[p],0);
					//	out << "\t" << test_model_change_below (ALI[a].phyl.nodelist[s].anc, &hitlist[h].ali, &ALI[a].phyl,0,&M[1],M[0].f[0],0,1);
					//}
				}
				else {
					//if (hitlist[h].refs==ref_num) 
					out << "\t"<<TO_statistic(&M[1], &M[0], &ALI[a].phyl, &hitlist[h].ali,0,hitlist[h].refs,1);

					out << "\t";
					if (C_MEM<2) {				
						out<<Tstat_pval(&M[1],hitlist[h].ali.asints,hitlist[h].ali.Nseqs,0,1,5+hitlist[h].refs);
					}
					else { out << "nc"; }
					out<< "\t";
					if	(C_MEM==0) {
						 out << Tstat_pval(&M[1],hitlist[h].ali.asints,hitlist[h].ali.Nseqs,0,1,5+hitlist[h].refs+ALI[0].Nseqs);
					}
					else { 	 out<<"nc"; }
					//else out <<"\tnr\tnr";
		//			for (s=0;s<hitlist[h].ali.Nseqs;s++) {
		//				out << "\t"<< test_leaf_loss(s,&hitlist[h].ali,&ALI[a].phyl,0,&M[1],M[0].f[p],1);
					//	out << "\t" << test_model_change_below (ALI[a].phyl.nodelist[s].anc, &hitlist[h].ali, &ALI[a].phyl,0,&M[1],M[0].f[0],1,1);
					//}
				
				}
				}
				else { out <<"\tg\tg\tg"; }
				out << endl;
			}
		}
		delete [] hitlist; delete [] ssp; // = new double [ALI[a].length];
	}
}

bool find_hits ( struct alignment * a, struct motif * m, struct hit * hlist ,
				double * sspvals, int * sssp, double cut, int * hnum, int l, int r) {

	bool find_top_hit ( double * ssplist, int * sssp, struct hit * tophit, int left, int right, double cut, struct alignment * ali, int w);
	bool make_hit_alignment ( struct alignment * a, struct motif * m, struct hit * h , int rightbound, int leftbound);
	//cout << "searching " << l << " to " << r <<endl;
	if (find_top_hit(sspvals,sssp,&hlist[hnum[0]],l,r,cut,a,m->length)) {
		//cout << "hit "<<hnum[0]<<" "<<l<< "-"<<r<<": "<<hlist[hnum[0]].leftend<< " in "<<hlist[hnum[0]].refs<<endl ;
		make_hit_alignment(a,m,&hlist[hnum[0]],l,r);
		//hlist[hnum[0]].print();
		int hitleft=hlist[hnum[0]].leftend;
		int hitright=hlist[hnum[0]].rightend;
		hnum[0]++;
		
		find_hits(a,m,hlist,sspvals,sssp,cut,hnum,l,hitleft-m->length);
		//find_hits(a,m,hlist,sspvals,sssp,cut,hnum,hitright,r-m->length);
		find_hits(a,m,hlist,sspvals,sssp,cut,hnum,hitright,r);
	}
	//else { cout << "none found"<<endl; }
	return 1;
}

bool make_hit_alignment ( struct alignment * a, struct motif * m, struct hit * h , int leftbound, int rightbound) {
	int s,p;
	int remove_gaps (struct sequence * seq, struct sequence * ug, int start, int size);
	int n_skipped_gaps (struct sequence * seq, int start, int w);
	h->ali.seq = new struct sequence [ a->Nseqs ];
	h->ali.Nseqs=a->Nseqs;	int ngaps_in_ref=remove_gaps (&a->seq[h->refs],&h->ali.seq[h->refs],h->apos[h->refs],m->length);
	h->aend = new int [ a->Nseqs ];
	h->aend[h->refs] = h->apos[h->refs] + m->length + ngaps_in_ref;
	//cout << h->ali.seq[h->refs].bases << endl;
	h->rightend=h->aend[h->refs]; h->leftend=h->apos[h->refs];

	h->strand = m->Pr(h->ali.seq[h->refs].asints,0)/(m->Pr(h->ali.seq[h->refs].asints,0)+m->rc->Pr(h->ali.seq[h->refs].asints,0));
	//cout << " strand = "<< strand << " " <<endl;
	//a->seq[h->refs].print(h->apos[h->refs],h->apos[h->refs]+m->length);
	double skip_gap_ss_pval_given_strand (struct motif * M, struct sequence * seq, int start, bool strand);
	
	for (s=0;s<a->Nseqs;s++) {//if (s!=h->refs) {

		int left = h->apos[h->refs]; 
		int right = (HEURISTIC) ? h->aend[h->refs]-1 : left+1;
		if (left <leftbound) cerr << "warning: "<<left <<" < "<<leftbound<<endl;
		int ugb=0; int bestp=left; double bestpval=1.0;
		//cout << s<<" "<<left <<endl;
		if (HEURISTIC) {
		
		while ((left>leftbound)&&(ugb<(m->length-1)))	{
			if (a->seq[s].is_base(left)) ugb++;	left--;
		}
		}
		//if (s==h->refs) cout << left <<" "<<leftbound ;
		//if (s==h->refs) cout << " "<<h->apos[h->refs] << endl;

		for (p=left;p<right;p++) {
			double ssp = (h->strand >= 0.5) ?
				skip_gap_ss_pval_given_strand (m, &a->seq[s], p, 1)
								:
				skip_gap_ss_pval_given_strand (m, &a->seq[s], p, 0);
			//if (s==h->refs) cout << p << " " << ssp <<endl;
			if (ssp<=bestpval) { //use the rightmost 
				int n = n_skipped_gaps(&a->seq[s],p,m->length);
				if ((p + n)<=rightbound) { bestp=p; bestpval=ssp; }
			}
		}
		h->apos[s] = h->aend[s] = bestp;
		h->aend[s] += m->length;
		int nsk = n_skipped_gaps(&a->seq[s],bestp,m->length);
		if ((nsk<50)&&(HEURISTIC)) h->aend[s] += remove_gaps (&a->seq[s],&h->ali.seq[s],bestp,m->length);
		else remove_gaps (&a->seq[s],&h->ali.seq[s],bestp,m->length);
		
		if (h->aend[s] > h->rightend) h->rightend=h->aend[s];
		if (h->apos[s] < h->leftend) h->leftend=h->apos[s];
		//if (h->refs==s) cout << "* ";
		//cout << h->apos[s]<<"-"<<h->aend[s]<< " "<< h->ali.seq[s].bases << " " << h->ali.seq[s].name << " "<<leftbound<< " "<<rightbound<<endl;

	}
	//cout <<h->sspval<<endl;
	h->ali.length=m->length;
	h->ali.integerize();

	//h->ali.print();
	return 1;
}

bool find_top_hit ( double * ssplist, int * sssp, struct hit * tophit, int left, int right, double cut, struct alignment * ali, int w) {
	int n_skipped_gaps (struct sequence * seq, int start, int w);
	int p; double bestpval=1.0; int bestp=-1;
//cout <<	"searching " << left << " to "<<right<<endl;
	for (p=left;p<right;p++) if ((ssplist[p]<bestpval) &&
	(ali->seq[sssp[p]].is_base(p))) { //don't start with a gap
		if ((n_skipped_gaps(&ali->seq[sssp[p]],p,w)+p)<right) {
			bestp=p; bestpval=ssplist[p];
			//cout << p << " : "<< ssplist[p] << " " <<sssp[p]<<endl;
		}
		
	}
	if (bestpval>=cut) return 0;
	//if (bestp==-1) return 0;
	tophit->refs = sssp[bestp];

	tophit->apos[tophit->refs]=bestp; tophit->sspval=bestpval;
	//cout << left<<"-"<<right<<" "<<tophit->apos[tophit->refs]<<" "<<bestpval<<" "<<cut<<endl;
	tophit->leftend=bestp; tophit->rightend=bestp;
	return 1;

}

void help () {
	//cout << "help contents here"<<endl;
	cout << "you are using Rmonkey v1.0"<<endl<<"usage: [matrix file] [fasta alignment] (options)"<<endl;
	cout << "\toption\texplanation"<<endl;
	//cout << "\t-m\tmotif model : mk,JC,HB or simple (default is HB)"<<endl;
	//cout << "\t-R\trate for JC motif model (default is 0.5)"<<endl;
	cout << "\t-b\tbackground model : JC,HKY (default is JC)"<<endl;
	cout << "\t-kappa\ttransition transversion rate ratio for HKY model (default is 2.0)"<<endl;
	cout << "\t-tree\ttree file name (default is alignment.tree) or pairwise distance (default is to compute it)"<<endl;
	cout << "\t-list\tsecond argument is a text file that is a list of fasta alignments (instead of a fasta alignment)"<<endl;
	cout << "\t-help\tdisplay this mesage"<<endl;
	cout << "\t-scr\tprint log LR scores (instead of pvalues)"<<endl;
	cout << "\t-freq\tbackground frequency file (default is to compute them from the sequences)"<<endl;
	cout << "\t-ref\treference sequence number (used in heuristics; default is 0, the first one in the file)"<<endl;
	cout << "\t-any\tallow match in any sequence to be the reference"<<endl;	
	cout << "\t-cut\tpvalue cutoff (default is 0.01)"<<endl;
	cout << "\t-out\t.occ output filename (default is fasta alignment.occ)"<<endl;
	cout << "\t-M\t 0,1,2 (save memory by not calculating some p-values 0=default, calculate all)"<<endl;
	cout << "\t-quiet\tprint less messages"<<endl;
	//cout << "\t-n\tmin number of sequences that have a match

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
		//cout << p << " "<< ali[p+pos]<<" "<<e_mod << " "<<m->pssm[e_mod][p][ali[p+pos]]<< " "<<m->pssm[e_mod][m->length-p-1][complem_int(ali[p+pos],ns)]<<endl;
	}

	return ( (scr>rcscr) ? m->pvalue[e_mod][scr] : m->pvalue[e_mod][rcscr] );
 }

double Tstat_pval(struct motif * m, long int * ali, int ns, int pos, bool strand,int e_mod) {
 	long int complem_int (long int i, int l);
 	int p, scr,rcscr;
	scr=rcscr=0;
	for (p=0;p<m->length;p++) if (ali[p+pos] !=-1) {
		if (strand==1) scr+=m->pssm[e_mod][p][ali[p+pos]];
		else scr+=m->pssm[e_mod][m->length-p-1][complem_int(ali[p+pos],ns)];
		//cout << p << " "<< ali[p+pos]<<" "<<m->pssm[5][p][ali[p+pos]]<< " "<<m->pssm[5][m->length-p-1][complem_int(ali[p+pos],ns)]<<endl;
	}
	return m->pvalue[e_mod][scr];
	//cout << scr << " " <<m->pvalue[5][scr]<< " "<<rcscr<<" "<<m->pvalue[5][rcscr]<<endl;
	//return ( (scr>rcscr) ? m->pvalue[5][scr] : m->pvalue[5][rcscr] );
 }

// double ss_pval_given_strand(struct motif * m, int * seq, int pos, bool strand) {
 //	int * C= new int [4]; C[0]=3; C[1]=2; C[2]=1; C[3]=0;
//	int p, scr,rcscr;
//	scr=rcscr=0;
//	for (p=0;p<m->length;p++) if (seq[p+pos] !=-1) {
//		scr += m->pssm[4][p][seq[p+pos]];
//		rcscr+=m->pssm[4][m->length-p-1][C[seq[p+pos]]];
//	}
//	delete [] C;
//	//cout << scr <<" "<< rcscr<< " "<<m->pvalue[4][scr]<<" "<<m->pvalue[4][rcscr]<<endl;
//	return ( (scr>rcscr) ? m->pvalue[4][scr] : m->pvalue[4][rcscr] );
// }

bool strand_w_min_pval (struct motif * m, int * seq, int pos) {
	int * C= new int [4]; C[0]=3; C[1]=2; C[2]=1; C[3]=0;
	int p, scr,rcscr;
	scr=rcscr=0;
	for (p=0;p<m->length;p++) if (seq[p+pos] !=-1) {
		scr += m->pssm[4][p][seq[p+pos]];
		rcscr+=m->pssm[4][m->length-p-1][C[seq[p+pos]]];
	}
	delete [] C;
	//cout << scr <<" "<< rcscr<< " "<<m->pvalue[4][scr]<<" "<<m->pvalue[4][rcscr]<<endl;
	return ( (scr>rcscr) ? 1 : 0 );
}

 double ss_pval(struct motif * m, int * seq, int pos) {
 	int * C= new int [4]; C[0]=3; C[1]=2; C[2]=1; C[3]=0;
	int p, scr,rcscr;
	scr=rcscr=0;
	for (p=0;p<m->length;p++) if (seq[p+pos] !=-1) {
		scr += m->pssm[4][p][seq[p+pos]];
		rcscr+=m->pssm[4][m->length-p-1][C[seq[p+pos]]];
	}
	delete [] C;
	//cout << scr <<" "<< rcscr<< " "<<m->pvalue[4][scr]<<" "<<m->pvalue[4][rcscr]<<endl;
	return ( (scr>rcscr) ? m->pvalue[4][scr] : m->pvalue[4][rcscr] );
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
		 	if (e_mod!=4) {
				scr+=m->pssm[e_mod][p][ali->asints[p+pos]];
				rcscr+=m->pssm[e_mod][m->length-p-1][complem_int(ali->asints[p+pos],ns)];
			}
			else {
				scr+=m->pssm[e_mod][p][ali->seq[ref].asints[p+pos]];
				rcscr+=m->pssm[e_mod][m->length-p-1][C[ali->seq[ref].asints[p+pos]]];
			}
			//scr+=m->pssm[e_mod][p][ali[p+pos]];
			//cout << ali[p+pos]<<" "<<gapped_complem_int(ali[p+pos],ns)<<endl;
			//rcscr+=m->pssm[e_mod][m->length-p-1][complem_int(ali[p+pos],ns)];
			//if (e_mod==4) cout << p << " "<< e_mod << " "<<m->pssm[e_mod][p][ali->asints[p+pos]]<< " "<<m->pssm[e_mod][m->length-p-1][complem_int(ali->asints[p+pos],ns)]<<endl;
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
					double mcond = p_subtree_jc(ali->seq,ali->phyl.root,b,p+pos,MODEL_ONE_RATE);
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
				//else cout << "this position if gapped"<<endl;
				scr += (int) (N_BINS*llr);
				rcscr+=(int) (N_BINS*rcllr);
				//int b;
				//cout << p <<" "<< ali->seq[ref].asints[p+pos]<<" " <<llr<<" "<< (((int) (N_BINS*llr)) - m->minscr[e_mod]) <<" "<<scr<< endl;

				//cout << m->f[p][b] << " "<<m->rc->f[m->length-p-1][C[b]]<<endl;
			}

			scr -= m->minscr[e_mod];
			rcscr -= m->minscr[e_mod];
		}
		//cout<< pos <<" " <<e_mod<< " "<<scr<< " "<<m->pvalue[e_mod][scr]<< " "<<rcscr<< " "<<m->pvalue[e_mod][rcscr]<<endl;
	}
	//cout<< pos <<" " <<e_mod<< " "<<m->minscr[e_mod]<<" "<<scr<< " "<<m->pvalue[e_mod][scr]<< " "<<rcscr<< " "<<m->pvalue[e_mod][rcscr]<<endl;
	delete [] C;
	return ( (scr>rcscr) ? m->pvalue[e_mod][scr] : m->pvalue[e_mod][rcscr] );
 }



double TO_statistic (struct motif * m, struct motif * bg, struct tree * t, struct alignment * ali,
					int pos, int ref, bool strand) {
	//int * C= new int [4]; C[0]=3; C[1]=2; C[2]=1; C[3]=0;
	double llh=0.0; double llhrc=0.0; int p; int b;
	for (p=0;p<m->length;p++) {
		 double pm=0.0; double pbg=0.0; double prc=0.0;
		for (b=0;b<4;b++) {
			if (strand) pm += m->f[p][b]*p_subtree_hb(ali->seq,t->root,b,pos+p,p,m);
			else pm += m->rc->f[p][b]*p_subtree_hb(ali->seq,t->root,b,pos+p,p,m->rc);
			pbg += bg->f[p][b]*p_subtree2(ali->seq,t->root,t,b,pos+p);
		}
		llh+=log(pm/pbg); //llhrc += log(prc/pbg);
		//cout << p<<" "<<pm << " "<<prc<<" "<<pbg<<endl;
	}
	//cout << llh<<" " <<llhrc<<endl;
	//return 1.0;

	double sscorr=0.0;
	double marginalized_px_hb (struct node * n, int p, int base, struct motif * m, struct node * x, int * seq);
	double marginalized_px_bg (struct node * n, int base, struct tree * t, struct node * x, int * seq);
	double p_above (int base, struct node * n, struct alignment * ali, struct tree * t, int pos, double * bgf);
	double  p_above_hb (struct sequence * seq, struct node * n, int base, int sp, int mp, struct motif * m);
	//struct node * sis = t->nodelist[ref].sister(); int c;
	for (p=0;p<m->length;p++) {
		int * col = new int [ali->Nseqs]; int s;
		for (s=0;s<ali->Nseqs;s++) col[s] = ali->seq[s].asints[p+pos];
		double prefgbg = 0.0;
		//int refb = ali->seq[ref].asints[p+pos];
		//cout <<p<<" " << refb << " " <<endl;

		for (b=0;b<4;b++) //{
			prefgbg += bg->f[p][b] * marginalized_px_bg (t->root, b,t, &t->nodelist[ref], col);
			//double sispr=0.0;
			//= p_subtree2(ali->seq, sis, t,
		//prefgbg += p_above(refb,&t->nodelist[ref],ali,t,pos+p,bg->f[p]);
		//}
		//cout << p << " "<<prefgbg<<endl;
		double prefgm=0.0;
		if (strand) 	{
			for (b=0;b<4;b++)
				prefgm += m->f[p][b] * marginalized_px_hb (t->root, p,b,m, &t->nodelist[ref], col);
				//prefgm  += p_above_hb (ali->seq, &t->nodelist[ref],refb, pos+p, p, m) ;
			//}
		}
		else {
			for (b=0;b<4;b++)
				prefgm += m->rc->f[p][b] * marginalized_px_hb (t->root, p,b,m->rc, &t->nodelist[ref], col);
				//prefgm += p_above_hb (ali->seq, &t->nodelist[ref],refb,pos+p,p,m->rc);
		}
		//cout << prefgm << " "<< prefgbg<<" "<<p<<endl;
		sscorr += log (prefgm/prefgbg);
	}
	//cout << llh <<" "<<llhrc <<" "<<sscorr<<endl;
	//if (llh>llhrc) return (llh - sscorr);
	//else 
	return (llh - sscorr);
	//double ssllr= log(m->Pr(ali->seq[ref].asints,0)/bg->Pr(ali->seq[ref].asints,0));
	//cout << " "<<sscorr <<" "<<ssllr<<endl;
	//return (llh - sscorr);
//
////
}


double p_above (int base, struct node * n, struct alignment * ali, struct tree * t,
					int pos, double * bgf) {
	if (n->root==1) return bgf[base];
	int b,c;
	double pr=0.0;
	struct node * sis = n->sister();
	for (b=0;b<4;b++)  {
		double ps=0.0;
		for (c=0;c<4;c++) ps+= (p_subtree2(ali->seq,sis,t,c,pos)
								* t->psubs[b][c][n->anc->number][sis->number]);
		pr += (ps * p_above(b,n->anc,ali,t,pos,bgf)
				* t->psubs[b][base][n->anc->number][n->number]);
	}
	//cout << base<<" "<<pr<<" "<<n->number<<" "<<sis->number<<endl;
	return pr;
}



//void internal_node_posteriors (struct node * n, double * post, struct alignment * ali, int pos,double * bgf) {
//	int b; double tot=0.0; double lh=0.0;
//	for (b=0;b<4;b++) post[b] = p_above(b,n,ali,pos,bgf)*p_subtree2(ali->seq, n, &ali->phyl, b, pos);
//	for (b=0;b<4;b++) { tot+=post[b]; lh += post[b]; }
//	for (b=0;b<4;b++) post[b] /= tot;
//	//cout <<lh;
//}


double  p_above_hb (struct sequence * seq, struct node * n, int base, int sp, int mp, struct motif * m) {
	if (n->root==1) return m->f[mp][base];
	int b,c; double pr=0.0;
	struct node * sis = n->sister();
	for (b=0;b<4;b++) {
		double ps = 0.0;
		for (c=0;c<4;c++) ps+= (p_subtree_hb(seq,sis,c,sp,mp,m) *
								m->psubs[mp][n->anc->number][sis->number][b][c] );
		pr+= (ps * p_above_hb(seq,n->anc,b,sp,mp,m)
				* m->psubs[mp][n->anc->number][n->number][b][base]);
	}
	return pr;
}

void internal_node_post_HB (struct node * n, double * post, struct alignment * ali, int pos, struct motif * m, int mp) {
	int b; double tot=0.0; double lh=0.0;
	for (b=0;b<4;b++) post[b] = p_above_hb(ali->seq,n,b,pos,mp,m) * p_subtree_hb(ali->seq,n,b,pos,mp,m);
	for (b=0;b<4;b++) { tot+=post[b]; lh += post[b]; }
	for (b=0;b<4;b++) post[b] /= tot;
	//cout <<lh;
}

double marginalized_px_hb (struct node * n, int p, int base, struct motif * m, struct node * x, int * seq) {
	bool is_parent (struct node * n, struct node * m);
	if (x->anc->number == n->number)
		return m->psubs[p][n->number][x->number][base][seq[x->number]];
	int b; double pr=0.0;
	if (is_parent(x,n->left)) for (b=0;b<4;b++)
		pr+= m->psubs[p][n->number][n->left->number][base][b] * marginalized_px_hb(n->left,p,b,m,x,seq);
	else if (is_parent(x,n->right)) for (b=0;b<4;b++)
		pr+= m->psubs[p][n->number][n->right->number][base][b]* marginalized_px_hb(n->right,p,b,m,x,seq);
	return pr;
}


double marginalized_px_bg (struct node * n, int base, struct tree * t, struct node * x, int * seq) {
	bool is_parent (struct node * n, struct node * m);
	if (x->anc->number == n->number)
		return t->psubs[base][seq[x->number]][n->number][x->number];
	int b; double pr=0.0;
	if (is_parent(x,n->right)) for (b=0;b<4;b++)
		pr+= t->psubs[base][b][n->number][n->right->number]*marginalized_px_bg(n->right,b,t,x,seq);
	else if (is_parent(x,n->left)) for (b=0;b<4;b++)
		pr+= t->psubs[base][b][n->number][n->left->number]*marginalized_px_bg(n->left,b,t,x,seq);
	return pr;
}

bool is_parent (struct node * n, struct node * m) { //is m a parent of n
	if (m->leaf) return 0;
	if (m->root) return 1;
	if (n->root) return 0;
	if (n->anc->number == m->number) return 1;
	if (is_parent(n,m->left)) return 1;
	if (is_parent(n,m->right)) return 1;
	return 0;
}


 void compute_pssm_pvals (struct motif * m, struct motif * bg, struct tree * t, int N, int e_mod, int ref) {
 	void generate_variant (long int mn, int ml, int * seq);
	//void general_variant(long int mn, int ml, int * seq, int N);
	double marginalized_px_hb (struct node * n, int p, int base, struct motif * m, struct node * x, int * seq);
	double marginalized_px_bg (struct node * n, int base, struct tree * t, struct node * x, int * seq);

	long int Nalpha = (long int) pow(4.0,N);
	if (e_mod==4) Nalpha=4;
	//long int Nalpha = (long int) pow(5.0,N);
	m->pssm[e_mod] = new int * [m->length];

	int p,b,s;
	long int i;
	int min=10000;
	//int worst=10000;
	double ** bgdist = new double * [m->length];
	double * bgpr = new double [Nalpha];
	//double
	int best=0;
	//cout <<endl;
	if (e_mod!=4) {
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
				col[s].length=1;
				col[s].asints[0]=variant[s];
				//if ((p==0)&&(e_mod==4)) cout << variant[s];
			}
			if (e_mod==3) m->pssm[e_mod][p][i] = (int) (N_BINS*(llr/((double) N)));
			//if (e_mod==4) m->pssm[e_mod][p][i] = (variant[ref] == 5) ? 0 : ((int) (N_BINS*log(m->f[p][variant[ref]]/bg->f[p][variant[ref]])));
			//else {
			//if (e_mod==2) { cout << " |"<<p<<"| "; for (b=0;b<4;b++) cout << m->f[p][b] << " ";}
			double pm=0.0;
			double pbg=0.0;
			double sscorr=0.0; // used only in e_mod==5
			//for (b=0;b<4;b++) {

				//calculate bg and bgdist
			if (p==0) {
				for (b=0;b<4;b++) pbg += bg->f[p][b]*p_subtree2(col,t->root,t,b,0);
				//if (e_mod==4) cout << "b"<<pbg<<" p"<<pow(bg->f[p][b], N)<<endl;
			}

			if (p==0) bgdist[p][i] = bgpr[i] = pbg;
			else bgdist[p][i] = pbg = bgpr[i];

				//calculate the motif probabilities
			if (e_mod==0) for (b=0;b<4;b++) pm += m->f[p][b]*p_subtree2(col,t->root,t,b,0);
			else if (e_mod==1) for (b=0;b<4;b++) pm += m->f[p][b]*p_subtree_jc(col,t->root,b,0,MODEL_ONE_RATE);
			else if (e_mod==2) for (b=0;b<4;b++) pm += m->f[p][b]*p_subtree_hb(col,t->root,b,0,p,m);
			else if (e_mod>=5) { //conditioned on reference models
				double prefgbg=0.0; double prefgm=0.0;
				for (b=0;b<4;b++) {
					pm += m->f[p][b]*p_subtree_hb(col,t->root,b,0,p,m);
					prefgbg += bg->f[p][b] * marginalized_px_bg (t->root, b,t, &t->nodelist[ref], variant);
					prefgm += m->f[p][b] * marginalized_px_hb (t->root, p,b,m, &t->nodelist[ref], variant);
				}
				sscorr+=log(prefgm/prefgbg);
				if ((e_mod>=5)&&(e_mod<(5+N))) bgdist[p][i]=pm; //use HB as bgdist
			}
			//for (b=0;b<4;b++) if (e_mod==2) cout <<t->nodelist[b].abr<< " "<<p_subtree_hb(col,t->root,b,0,p,m) <<" ";
			//cout <<endl;
			//}
			//if (e_mod==2) cout << p << " "<<pm <<" "<< pbg << endl	;
			//check += pbg;
			if (e_mod<3) m->pssm[e_mod][p][i] = ((int) (N_BINS*log(pm/pbg)));
			if ((e_mod>=5)&&(e_mod<(5+N))) m->pssm[e_mod][p][i] = ((int) (-1.0*N_BINS*(log(pm/pbg) - sscorr)));
			//reverse the sign for <=T
			if (e_mod>=(5+N)) m->pssm[e_mod][p][i] = ((int) (N_BINS*(log(pm/pbg) - sscorr)));
				
			//if ((p==5)&&(e_mod==2))
			//	cout <<" "<< i <<" "<< p<<" "<<variant[ref]<< " "<<pm<<" "<<bgdist[p][i]<<" "<<m->pssm[e_mod][p][i]<<endl;

			if (m->pssm[e_mod][p][i]<min) min=m->pssm[e_mod][p][i];
		}
		//cout << check << "=? 1"<<endl;
	}
	}
	else {
		for (p=0;p<m->length;p++) {
			m->pssm[e_mod][p] = new int [Nalpha];
			bgdist[p] = new double [Nalpha];
			for (i=0;i<Nalpha;i++) {
				bgdist[p][i] = bg->f[p][i];
				m->pssm[e_mod][p][i] = ((int) (N_BINS*log(m->f[p][i]/bg->f[p][i])));
				//cout << m->pssm[e_mod][p][i] << " ";
				if (m->pssm[e_mod][p][i]<min) min=m->pssm[e_mod][p][i];
			}
			//cout << endl;
		}

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
	//if (e_mod==2) for (p=0;p<m->length;p++) { for (i=0;i<Nalpha;i++) cout << m->pssm[e_mod][p][i]<< " "; cout <<endl; }
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
		//if ((e_mod==9)&&(score_pdf(s,m->length-1,Nalpha,bgdist,m->pssm[e_mod],pdf)!=0.0)) cout << s <<" "<<score_pdf(s,m->length-1,Nalpha,bgdist,m->pssm[e_mod],pdf)<< " "<<m->pvalue[e_mod][s]<<endl;
		//cout << s << " "<< pdf[2][s]<<endl;
	}
	delete [] pdf;
	delete [] bgpr;
	delete [] bgdist;
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
	delete [] asints;
	return i;
 }

 long int gapped_col_to_int (struct sequence * seq, int p, int n) {
 	long int gapped_seq_to_num(int * seq, int start, int ml);
	int s; int * asints = new int [n];
	for (s=0;s<n;s++) asints[s] = seq[s].asints[p];
	long int i = gapped_seq_to_num(asints,0,n);
	delete [] asints;
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
	delete [] seq;
	delete [] com;
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
	delete [] seq;
	delete [] com;
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
				double tot=0.0;
				for (a=0;a<4;a++) {
					  Q[a] = new double [4];
					for (b=0;b<4;b++) if (a != b) {
						Q[a][b] = bgf[b];
						if (((a==0)&&(b==2)) || ((a==2)&&(b==0)) || ((a==1)&&(b==3)) || ((a==3)&&(b==1)))
							Q[a][b] *= kap;
						tot+= (bgf[a]*Q[a][b]);
					}
				}
				for (a=0;a<4;a++) for (b=0;b<4;b++)
					if (a!=b) Q[a][b] *= (t->nodelist[n].abr/tot);
				for (a=0;a<4;a++) {
					double idrate=0.0;
					for (b=0;b<4;b++) if (a!=b) {
						double temp = (m->f[p][b]*Q[b][a])/(m->f[p][a]*Q[a][b]);
						//double temp = 0.5;
						if (fabs(1.0-(1.0/temp))>1e-100) R[a][b] = Q[a][b] * (log(temp)/(1.0 - (1.0/temp)));
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

int optimized_seq(struct sequence * ref, struct sequence * orth, struct sequence * news,int start, int length, int offset) {

	//					****						//
	//this is substantially different from monkey1.0//
	//mix and match with caution					//
	//					****						//

	int s;
	//int  (int * seq, int start, int stop);
	int remove_gaps (struct sequence * seq, struct sequence * ug, int start, int size);
	//struct sequence remove_gaps (struct sequence * seq, int start, int size);
	//double pw_pid (struct sequence * seq1,struct sequence * seq2, int start1, int start2, int length);

	//struct sequence bestseq;
	int minus,plus,ngaps;
	ngaps = gap_count(orth->asints,start,start+length);

	if ((ngaps>MAX_GAPS)||(offset>MAX_GAPS)) {
	  news->length=0;
	  return -1;
	}

	struct sequence refseq;
	remove_gaps(ref,&refseq,0,length);
	//cout <<start<< " "<< ref->bases << endl;
	minus=plus=offset+ngaps;
	int begin=start;
	if (minus==0) { remove_gaps(orth, news, start, length); refseq.erase(); return start; }
	else {
		//cout << ngaps << " gaps and "<<offset<< " gaps in reference" <<endl;
		int ugbp=0;
		int end = ((start+plus+length)<=orth->length) ? (start+plus) : (orth->length-length) ;
		while ((begin>0)&&(ugbp<minus)) {
			begin--;
			if (orth->is_base(begin)) ugbp++;
		}
			//find the best ungapped local alignment
		//cout << ref->name << endl;
		//cout << "looking from " << begin << " to " << end <<" to match "<<refseq.bases<< " in "<<orth->print(begin,end+length);
		double bestscr=0.0;
		int bestq=-1; int q;
		for (q=begin;q<=end;q++) if (orth->is_base(q)) {
			struct sequence temp;
			remove_gaps(orth,&temp,q,length);
			if (temp.length==length) {
				//remove_gaps(orth,news,q,length);
				//return;
				double scr=pw_pid(&refseq,&temp,0,0,length);
				if (scr>bestscr) { bestq=q; bestscr=scr; }
				else if ((scr == bestscr)&&(abs(bestq-start)>abs(q-start))) { bestq=q; bestscr=scr; }
			}
			//cout << q << " " << scr << " " <<temp.bases<<endl;
		}
		refseq.erase();
		if (bestq==-1) return -1;
		else {	//store the best ungapped wmer
			remove_gaps(orth, news,bestq, length);
			//cout <<news->bases<< " chosen"<<endl;
			return bestq;
		}

		//cout << news->bases << " is the sequence"<<endl;
	}
	//delete refseq;
	//cout << "erasing ";

	//cout << " done"<<endl;


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


int remove_gaps (struct sequence * seq, struct sequence * ug, int start, int size) {
	//struct sequence ug;
	ug->name = seq->name + "_sub"; ug->bases.erase();
	int pos=start; int skipped=0;
	//cout << pos << " "<<ug->name<<endl;
	while (pos<seq->length) {
		if (seq->is_base(pos)) ug->bases += seq->bases[pos];
		else skipped++;
		if ((skipped > 50)||(ug->bases.size()>=size)) pos=seq->length;
		//cout << pos <<" "<<ug->bases <<endl;
		pos++;
	}
	//cout << " out of the loop"<<endl;

	ug->integerize();
	//ug->print(0,ug->length);
	return skipped;
}

double skip_gap_ss_pval_given_strand (struct motif * M, struct sequence * seq, int start, bool strand) {
	int pos=0; int skipped=0;
	int * C= new int [4]; C[0]=3; C[1]=2; C[2]=1; C[3]=0;
	int scr=0; int rcscr=0;
	while (pos<M->length) {
		int i=start+pos+skipped;
		if ((i >= seq->length)||(skipped>50)) return 100.0;
		else if (seq->is_base(i)) {
			if (strand) scr += M->pssm[4][pos][seq->asints[i]];
			else scr+=M->pssm[4][M->length-pos-1][C[seq->asints[i]]];
			pos++;
		}
		else skipped++;
	}
	delete [] C;
	//cout << scr << " "<<rcscr<< " "<<M->pvalue[4][scr]<<" "<<M->pvalue[4][rcscr]<<endl;
	//if (rcscr>scr) return M->pvalue[4][rcscr];
	//else
	return M->pvalue[4][scr];
}

int n_skipped_gaps (struct sequence * seq, int start, int w) {
	int pos=0; int skipped=0;
	while (pos<w) {
		int i=start+pos+skipped;
		if ((i >= seq->length)||(skipped>50)) return skipped;
		if (seq->is_base(i)) pos++;
		else skipped++;
	}
	return skipped;
}

double skip_gap_ss_pval (struct motif * M,  struct sequence * seq, int start) {
	int pos=0; int skipped=0;
	int * C= new int [4]; C[0]=3; C[1]=2; C[2]=1; C[3]=0;
	int scr=0; int rcscr=0;
	while (pos<M->length) {
		int i=start+pos+skipped;
		if ((i >= seq->length)||(skipped>50)) return 100.0;
		else if (seq->is_base(i)) {
			scr += M->pssm[4][pos][seq->asints[i]];
			rcscr+=M->pssm[4][M->length-pos-1][C[seq->asints[i]]];
			pos++;
		}
		else skipped++;
	}
	delete [] C;
	//cout << scr << " "<<rcscr<< " "<<M->pvalue[4][scr]<<" "<<M->pvalue[4][rcscr]<<endl;
	if (rcscr>scr) return M->pvalue[4][rcscr];
	else return M->pvalue[4][scr];
}

double skip_gap_log_LR (struct motif * M, struct motif * B, struct sequence * seq, int start) {
	int pos=0; int skipped=0;
	int * C= new int [4]; C[0]=3; C[1]=2; C[2]=1; C[3]=0;
	double llr=0.0; double rcllr=0.0;
	while (pos<M->length) {
		int i=start+pos+skipped;
		if ((i >= seq->length)||(skipped>50)) return -1000.0;
		else if (seq->is_base(i)) {
			llr += log( M->f[pos][seq->asints[i]] / B->f[pos][seq->asints[i]] );
			rcllr += log( M->rc->f[pos][seq->asints[i]] / B->f[pos][C[seq->asints[i]]] );
			pos++;
		}
		else skipped++;
	}
	delete [] C;
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

	int n;
	pm=prc=pbg=1.0;

	if (e_mod==3) {
		int s,ns=0; for (n=0;n<phy->Nnodes;n++) if (phy->nodelist[n].leaf) ns++;
		lhbg=lhm=lhrc=1.0;
		for (p=0;p<m->length;p++) for (s=0;s<ns;s++) if (seq[s].is_base(p+pos)) {
			pm *= m->f[p][seq[s].asints[p+pos]];
			pbg *= bgf[seq[s].asints[p+pos]];
			if (BOTH_STRANDS) prc *=  m->rc->f[p][seq[s].asints[p+pos]];
		}
		if ((BOTH_STRANDS)&&(prc>pm)) return ((log(prc/pbg))/((double) ns));
		return ((log(pm/pbg))/((double) ns));
	}
	else {
	double * bcond = new double [4];
	double * mcond = new double [4];
	double * rcond = new double [4];
	for (p=0;p<m->length;p++) {
		for (b=0;b<4;b++) bcond[b]= p_subtree2(seq,phy->root,phy,b,p+pos);
		if (e_mod==0) { for (b=0;b<4;b++) mcond[b]=rcond[b]=bcond[b]; }
		else if (e_mod==1) { for (b=0;b<4;b++) mcond[b]=rcond[b]=p_subtree_jc(seq,phy->root,b,p+pos,MODEL_ONE_RATE); }
		else if (e_mod==2) { for (b=0;b<4;b++) mcond[b]=p_subtree_hb(seq,phy->root,b,p+pos,p,m); }
		lhbg=lhm=lhrc=0.0;
		for (b=0;b<4;b++) {
			lhbg+= bgf[b] * bcond[b];
			lhm+= m->f[p][b] * mcond[b];
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
	delete [] rcond; delete [] mcond; delete [] bcond;
	if ((BOTH_STRANDS)&&(prc>pm)) return (log(prc)-log(pbg));
	return (log(pm)-log(pbg));
	}
	//if (PRINT_TO_SCREEN) {
	//	cout << e_mod<<": pm="<<log(pm) << " (prc="<<log(prc)<<") pbg="<< log(pbg) <<"LR= " ;
	//	if ((BOTH_STRANDS)&&(prc>pm)) cout << log(prc)-log(pbg);
	//	else cout << log(pm)-log(pbg);
	//	cout <<endl;
	//}


}

double PJC (int a, int b, double alpha) { //jukes-cantor probability of a to b in alpha//
	//cout << a<<"->"<<b<<endl;
	if (a!=b) return (0.25-0.25*exp((-4.0/3.0)*alpha));
	else return (0.25+0.75*exp((-4.0/3.0)*alpha));

}

hit::hit() {
	sspval=evolpval=1.0;
	ali.Nseqs=0;
}

bool hit::print() {
	cout << apos[refs] <<"("<<leftend<<"-"<<rightend << ") "<<sspval<<endl;
	if (ali.Nseqs>0) {
		int s;
		for (s=0;s<ali.Nseqs;s++)
			cout << apos[s] << " "<<ali.seq[s].bases<< " "<<ali.seq[s].name<<endl;
	}
	return 1;
}

sequence::sequence() {
	length=0;
}

//sequence::~sequence() {
	//cout << "deleting "<<name<<" "<<length<<"bp ...";
	//name.erase();
  	//bases.erase();

  	//delete [] asints;
  	//cout<<"done"<<endl;
	//erase();
	//if (length>0) delete [] asints;
//}
//double sequence::Pr(int start, int stop, double * f) {
//	int p; double pr=1.0;
//	for (p=start;p<stop;p++) if ((asints[p]!=4)&&(asints[p]!=-1)) pr *= f[p];
//	return pr;
//}
void sequence::erase() {
	name.erase();
  bases.erase();
 // cerr << " about to erase...";
  delete [] asints;
 // cerr << " done"<<endl;
  return;
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
	//cout << length << endl;
	asints = new int [length];
	//cout << "declared " <<endl;
	//B[0]='A'; B[1]='C'; B[2]='G'; B[3]='T';
	//for (p=0;p<length;p++) for (b=0;b<4;b++) if (B[b] == bases[p]) asints[p]=b;
	//for (p=0;p<length;p++) cout <<bases[p]<<endl;
	//cout << bases <<endl;
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
	if ((pos<0)||(pos>=length)) return 0;
	if ((asints[pos]==4)||(asints[pos]==-1)) return 0;
	else return 1;
}
int sequence::next_base(int pos) {
	while (pos<length) {
		if (is_base(pos)) return pos;
		pos++;
	}
	return -1;
}
int sequence::ungapped_pos(int pos) {
	int p; int ugp=pos;
	for (p=0;p<pos;p++) if (asints[p]==4) ugp--;
	//for (p=0;p<pos;p++) cout << asints[p] <<" ";
	return ugp;
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

struct node * node::sister() {
	return (anc->left->number == number) ? anc->right : anc->left;
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


bool tree::copy_topology(struct tree * t) {
	int i;
	Nnodes=t->Nnodes;
	nodelist = new struct node [Nnodes];
	for (i=0;i<Nnodes;i++) {
		nodelist[i].number = t->nodelist[i].number;
		if (t->nodelist[i].leaf) { nodelist[i].leaf=1; nodelist[i].abr = t->nodelist[i].abr; }
		else {
			nodelist[i].left = &nodelist[t->nodelist[i].left->number];
			nodelist[i].lbr = t->nodelist[i].lbr;
			nodelist[i].right = &nodelist[t->nodelist[i].right->number];
			nodelist[i].rbr = t->nodelist[i].rbr;
			nodelist[i].leaf=0;
			nodelist[i].abr = t->nodelist[i].abr;
		}
		if (t->nodelist[i].root) { root = &nodelist[i]; nodelist[i].root=1; }
		else { nodelist[i].anc = &nodelist[t->nodelist[i].anc->number]; nodelist[i].root=0;}
	}
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
		double tot=0.0;
		for (a=0;a<4;a++) {

			R[a] = new double [4];
			P[a] = new double [4];
			for (b=0;b<4;b++) if (a != b) {
				R[a][b] = f[b];
				//R[a][b]=1.0;
				if (((a==0)&&(b==2)) || ((a==2)&&(b==0)) || ((a==1)&&(b==3)) || ((a==3)&&(b==1)))
					R[a][b] *= kap;
				tot+=(f[a]*R[a][b]);
			}
		}
		for (a=0;a<4;a++) {
			double idrate=0.0;
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
		delete [] R;
		delete [] P;
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
	//cout << " reading tree file"<<endl;
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
