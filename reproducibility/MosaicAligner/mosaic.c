#include "mosaic.h"

main(int argc, char *argv[]) {

	int i;
	long seed = -setseed();
	char fname[MAXNAME];
	struct data *my_data;
	struct pars *my_pars;
	FILE *ifp;

	if (argc<2) {
		printf("\n\nInput file name:");
		scanf("%s", &fname);
		ifp = fopen(fname, "r");
	}
	else ifp = fopen(argv[1], "r");
	if (!ifp) nrerror("Cannot open input file");

	my_data = read_fasta(ifp, 1, 0);  /*The 1 indicates it is amino acid data*/
	printf("\n\nRead %i sequences\n\n", my_data->nseq);
	if (DEBUG) print_sequences(my_data, stdout);

	my_pars = (struct pars *) malloc((size_t) sizeof(struct pars));
	get_pars(my_pars, my_data);

	for (i=1;i<=1/*my_data->nseq*/;i++) {
		align(my_data, my_pars, i);
	}

	exit(0);

}



void align(struct data *my_data, struct pars *my_pars, int target) {

	struct palign **align_set;
	struct kalign *kwise;

	align_set = (struct palign **) malloc((size_t) my_data->nseq*sizeof(struct palign));
	
	/*Do pairwise alignments*/
	/*
	for (i=1;i<=my_data->nseq;i++) if (i!= target) {
		align_set[i] = pair_align_viterbi(my_data, my_pars, target, i);
	}
	*/

	/*Do kwise alignment*/
	kwise = kalign(my_data, my_pars, target);


	
}



/*Pairwise alignment with Viterbi - makes target (s1) as mosaic of s1*/

struct palign * pair_align_viterbi(struct data *my_data, struct pars *my_pars, int s1, int s2) {


	int i, j, k, l1, l2, **trace, *sq1, *sq2, mxi, len, lx, max_pc_i;
	double **m, val[5], max_pc;
	struct palign *my_align;

	my_align = (struct palign *) malloc((size_t) sizeof(struct palign));

	l1 = my_data->seqs[s1]->length;
	l2 = my_data->seqs[s2]->length;
	sq1 = my_data->seqs[s1]->seq;
	sq2 = my_data->seqs[s2]->seq;

	my_align->col_max = dvector(1, l1);
	my_align->col_pos = ivector(1, l1);

	/*Initialise matrices: Note putative mosaic makes columns for clarity*/
	m =  dmatrix(0, l2, 0, l1);
	trace = imatrix(0, l2, 0, l1);

	/*Note that j always refers to columns, i to rows*/
	m[0][0]=0.0;
	trace[0][0]=1;
	for (j=1;j<=l1;j++) {m[0][j]= (double) j*my_pars->del; trace[0][j]=3;}
	for (i=1;i<=l2;i++) {m[i][0]= 0; trace[i][0]=2;}

	/*Identifiers for best previous column*/
	max_pc = 0.0;
	max_pc_i = 0;

	/*In trace - 1=match, 2=move up a row(del), 3=move left a column(insert), -x=recombine to position x in previous column (NB can be zero)*/
	for (j=1;j<=l1;j++) {/*Do alignment by columns = j = putative mosaic*/
		for (i=1;i<=l2;i++) {
			val[1]=m[i-1][j-1]+my_pars->s[sq1[j]][sq2[i]];			/*Match*/
			val[2]=m[i-1][j]+my_pars->del;							/*Delete*/
			val[3]=m[i][j-1]+my_pars->del;							/*Insert*/
			val[4]=max_pc+my_pars->rho+my_pars->s[sq1[j]][sq2[i]];	/*Recombination*/

			/*Select best move*/
			for (k=2,val[0]=val[1],mxi=1;k<=4;k++) if (val[k]>val[0]) {val[0]=val[k]; mxi=k;}

			m[i][j]=val[mxi];
			if (mxi<4) trace[i][j]=mxi;
			else trace[i][j]=-max_pc_i;
		}

		for (i=1, max_pc_i=0, max_pc=m[0][j]; i<=l2;i++) if (m[i][j]>max_pc) {max_pc=m[i][j]; max_pc_i=i;}
	}

	if (DEBUG) {
		printf("\n\nAlignment matrix\n\n");
		printf("          ");
		for (j=1;j<=l1;j++) printf("%5c",num2nuc(sq1[j], 1));
		printf("\n     ");
		for (j=0;j<=l1;j++) printf("%5.0lf",m[0][j]);
		for (i=1;i<=l2;i++) {
			printf("\n%5c",num2nuc(sq2[i], 1));
			for (j=0;j<=l1;j++) printf("%5.0lf",m[i][j]);
		}
		printf("\n\n");

		printf("Traceback matrix\n\n");
		printf("          ");
		for (j=1;j<=l1;j++) printf("%5c",num2nuc(sq1[j], 1));
		printf("\n     ");
		for (j=0;j<=l1;j++) printf("%5i",trace[0][j]);
		for (i=1;i<=l2;i++) {
			printf("\n%5c",num2nuc(sq2[i], 1));
			for (j=0;j<=l1;j++) printf("%5i",trace[i][j]);
		}
	}

/*Find position in last colum to start alignment nd calculate length of alignment*/
	for(i=1, mxi=0, val[0]=m[0][l1];i<=l2;i++) if (m[i][l1]>val[0]) {mxi=i; val[0]=m[i][l1];}
	i=mxi;
	j=l1;
	len=0;

	while (j>0) {
		if (trace[i][j]==1) {i-=1; j-=1;len++;}			/*Match*/
		else if (trace[i][j]==2) {i-=1; len++;}			/*Delete*/
		else if (trace[i][j]==3) {j-=1; len++;}			/*Insert*/
		else {j-=1;i=-trace[i][j]; len++;}				/*Recombination*/
	}
	my_align->global_trace = imatrix(1, len, 1, 3);
	my_align->global_val = dvector(1, len);

	i=mxi;j=l1;lx=len;
	my_align->length=len;
	my_align->s1=s1;
	my_align->s2=s2;

	while (j>0) {
		my_align->global_val[len]=m[i][j]; 
		my_align->global_trace[len][1]=trace[i][j];
		my_align->global_trace[len][2]=i;
		my_align->global_trace[len][3]=j;

		if (trace[i][j]==1) {i-=1; j-=1;}
		else if (trace[i][j]==2) {i-=1;}
		else if (trace[i][j]==3) {j-=1;}
		else {i=-trace[i][j]; j-=1;}
		len--;
	}

	if (DEBUG) {
		printf("\n\nScores and traces of best alignment\n\n  Tr   i   j  Value\n\n");
		for (i=1;i<=lx;i++) printf("%4i%4i%4i%7.0lf\n",my_align->global_trace[i][1], my_align->global_trace[i][2],my_align->global_trace[i][3], my_align->global_val[i]);

		print_pair_align(my_align, sq1, sq2);
	}

	free_dmatrix(m, 0, l2, 0, l1);
	free_imatrix(trace, 0, l2, 0, l1);

	return my_align;
}


/*
Set parameters: Need to update to use doubles for costs, but still OK as is for simple Viterbi
*/

void get_pars(struct pars *par, struct data *my_data) {

	int i, j;
	FILE *ifp;

	if (my_data->type==1) {
		par->s = dmatrix(1, 25, 1, 25);
		for (i=1;i<=25;i++) for (j=1;j<=25;j++) {
			if (i<21 && j<21) par->s[i][j]=bl62[i-1][j-1];
			else par->s[i][j]=0;
		}
	}
	else {
		par->s = dmatrix(1, 4, 1, 4);
		for (i=1;i<=4;i++) for (j=1;j<=4;j++) par->s[i][j]=-1;
		for (i=1;i<=4;i++) par->s[i][j]=4;
	}

	ifp = fopen("params.txt", "r");
	if(!ifp) {

		par->del=-3;
		par->rho=-20;

		printf("\n\nSetting values to defaults: rho=%.1lf, delta=%.1lf\nBLOSUM matrix for AA, DIV for NTs\n\n", par->rho, par->del);

	}
	else {
		nrerror("Reading user-defined inputs - not yet implemented");
	}
}



struct kalign * kalign(struct data *my_data, struct pars *my_pars, int target) {

	int i, j, k, l, **tnew, ***t, **max_pc_id, *sq1, *sq2, mxi, len;
	double **mnew, ***m, max_pc, val[5], max_c;
	struct kalign *kwise;

	kwise = (struct kalign *) malloc((size_t) sizeof(struct kalign));


	/*Intialise arrays of alignment matrices*/
	m = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	t = (int ***) malloc((size_t) (my_data->nseq+1)*sizeof(int **));

	for (k=1;k<=my_data->nseq;k++) if (k != target) {
		mnew = dmatrix(0, my_data->seqs[k]->length, 0, my_data->seqs[target]->length);
		tnew = imatrix(0, my_data->seqs[k]->length, 0, my_data->seqs[target]->length);
		for (j=0;j<=my_data->seqs[target]->length;j++) {mnew[0][j]=(double) j*my_pars->del; tnew[0][j]=1;}
		for (i=1;i<=my_data->seqs[k]->length;i++) {mnew[i][0]=0; tnew[i][0]=2;}
		tnew[0][0]=1;

		m[k] = mnew;
		t[k] = tnew;
	}

	max_pc = max_c = 0.0;
	max_pc_id = imatrix(1,2,0,my_data->seqs[target]->length);
	max_pc_id[1][0]=max_pc_id[2][0]=0; /*Sequence and row of previous column's best value*/
	sq1 = my_data->seqs[target]->seq;

	for (j=1;j<=my_data->seqs[target]->length;j++) {

		/*Stuff to identify next best row for column*/
		max_pc = max_c;
		max_c = (double) j*my_pars->del;
		if (target==1) max_pc_id[1][j]=2;
		else max_pc_id[1][j]=1;
		max_pc_id[2][j]=0;

		for (k=1;k<=my_data->nseq;k++) if (k != target) {
			sq2 = my_data->seqs[k]->seq;
			for (i=1;i<=my_data->seqs[k]->length;i++) {
				val[1]=m[k][i-1][j-1]+my_pars->s[sq1[j]][sq2[i]];		/*Match*/
				val[2]=m[k][i-1][j]+my_pars->del;						/*Delete*/
				val[3]=m[k][i][j-1]+my_pars->del;						/*Insert*/
				val[4]=max_pc+my_pars->rho+my_pars->s[sq1[j]][sq2[i]];	/*Recombine*/

				for (l=2, val[0]=val[1], mxi=1; l<=4; l++) if (val[l]>val[0]) {val[0]=val[l]; mxi=l;}
				m[k][i][j]=val[mxi];
				t[k][i][j]=mxi;

				if (m[k][i][j]>max_c) {max_c=m[k][i][j]; max_pc_id[1][j]=k; max_pc_id[2][j]=i;}
			}
		}
	}

	/*NB max_pc_id will have identified place to start alignment*/
	/*Following first calculates length of alignment and then stores it in kwise*/
	j = my_data->seqs[target]->length;
	i = max_pc_id[2][my_data->seqs[target]->length];
	k = max_pc_id[1][my_data->seqs[target]->length];
	len=0;

	while (j>0) {
		len++;
		if (t[k][i][j]==1) {i--; j--;}
		else if (t[k][i][j]==2) {i--;}
		else if (t[k][i][j]==3) {j--;}
		else {j--; i=max_pc_id[2][j]; k=max_pc_id[1][j];}
	}

	kwise->target = target;
	kwise->length = len;
	kwise->trace = imatrix(1, 4, 1, len);
	kwise->score = dvector(1, len);

	j = my_data->seqs[target]->length;
	i = max_pc_id[2][my_data->seqs[target]->length];
	k = max_pc_id[1][my_data->seqs[target]->length];
	l = len;

	while (j>0) {
		kwise->trace[1][l]=t[k][i][j];
		kwise->trace[2][l]=k;
		kwise->trace[3][l]=i;
		kwise->trace[4][l]=j;
		l--;
		if (t[k][i][j]==1) {i--; j--;}
		else if (t[k][i][j]==2) {i--;}
		else if (t[k][i][j]==3) {j--;}
		else {j--; i=max_pc_id[2][j]; k=max_pc_id[1][j];}
	}

	if (DEBUG) {
	}

	print_kalign(kwise, my_data, stdout);


	/*Free up memory at end of routine*/
	for (i=1;i<=my_data->nseq;i++) if (i != target) {
		free_dmatrix(m[i], 0, my_data->seqs[i]->length, 0, my_data->seqs[target]->length);
		free_imatrix(t[i], 0, my_data->seqs[i]->length, 0, my_data->seqs[target]->length);
	}
	free(m);
	free(t);
	free_imatrix(max_pc_id, 1, 2, 0, my_data->seqs[target]->length);

	return kwise;

}



