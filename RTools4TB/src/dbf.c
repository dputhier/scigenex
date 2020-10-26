/******************************************************************************
 C implementation of the DBF algorithm

 COMPILATION : R CMD SHLIB -o dbf dbf.c
 gcc -O3 -lm -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -o dbf dbf.c -pg -g

 Authors :
   J. TEXTORIS ; F. LOPEZ : 05-30-2007
   D Puthier : 10-18-2007 (!)
   A. Bergon : 11-26-2008
   S. Granjeaud : 02-20-2009
         o Better documentation (english part)
         o Misc performance tips (replace_root...)
		 o Verbose and clock outputs
		 o Command line with options (see main)
		 o Optional richer outputs (see fprint_selected)
   A.Bergon : 04-14-2009
         o add new distance computation
                 o spm : arithmetic mean of spearman & pearson distances [ (s+p)/2 ]
                 o spgm : geometric mean of spearman & pearson distances [ sqrt(s*p) ]
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.


                   _|_|_|_|_|    _|_|      _|_|_|    _|_|_|
                       _|      _|    _|  _|        _|
                       _|      _|_|_|_|  _|  _|_|  _|
                       _|      _|    _|  _|    _|  _|
                       _|      _|    _|    _|_|_|    _|_|_|


******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <search.h>

#define debugLevel 0

/* Structure that keeps track of value (distance) and its index or rank
*/
typedef struct {
	double val;
	int i;
} RANG_VAL ;

/*==========================================================================*/
/* R helper functions
   A matrix in C is coded here as mat[][]; mat[] is a vector of pointers to each
   line of data. Thus, data are arranged line-by-line.
   A matrix in R is coded as a contiguous area of data, arranged column-by-column.
*/

// Read a matrix in R format
double **read_double_matrix(int nl, int nc, double *m) {
	int i, L, C;
	double **mat = (double **)malloc(nl * sizeof(double *));

	for (i = 0; i < nl; i++) *(mat + i) = (double *)malloc(nc * sizeof(double));
	for (L = 0; L < nl; L++)
		for (C = 0; C < nc; C++)
			*(*(mat + L) + C) = m[L + C * nl];
	return mat;
}

// Write a matrix in R format
void write_double_matrix(int nl, int nc, double **mat, double *m) {
	int L, C;

	for (L = 0; L < nl; L++)
		for (C = 0; C < nc; C++)
			*(m + C * nl + L) = mat[L][C];
}

// Read a matrix in R format
char ***read_char_matrix(int nl, int nc, char **m) {
	int i, L, C;
	char ***mat = (char ***)malloc(nl * sizeof(char **));

	for (i = 0; i < nl; i++) *(mat + i) = (char **)malloc(nc * sizeof(char*));
	for (L = 0; L < nl; L++)
		for (C = 0; C < nc; C++)
			*(*(mat + L) + C) = strdup(m[L + C * nl]);
	return mat;
}

// Write a matrix in R format
void write_char_matrix(int nl, int nc, char ***mat, char **m) {
	int L, C;

	for (L = 0; L < nl; L++)
		for (C = 0; C < nc; C++)
			*(m + C * nl + L) = mat[L][C];
}

/*==========================================================================*/
/* Heap sort and K lowest retaining
   A heap is a binary tree; value of each child is less than the value of
   the parent; but values in a branch could be greater than some values in
   another branch.
   Functions to manage (double,index) pairs are extensively commented in the
   following lines. Functions to manage double are then defined.
*/

/* Add a value,index pair to the heap
   The heap is not full; the new value is added at the end of the heap; because
   the bottom of the heap is the lowest value in the branche, the new value is
   compared to all values in the branches and swaped with any value lower than
   itself; this process stops when the new value is at its right place in the
   heap.
   liste : pointer to the heap
   nb : the current heap size, ie the position where to insert the new value
   val, i : the (value,i) to be inserted in the heap
*/
void add_value(RANG_VAL *liste, int nb, double val, int i) {
	int f, tmp_i;
	double tmp_val;

	// Add at the end
	liste[nb].val = val;
	liste[nb].i = i;
	// Loop to find the right position in the heap (from bottom to top)
	while (nb != 0) {
		// Get the parent position in the heap
		f = (nb % 2) ? ((nb - 1) / 2) : (nb / 2 - 1) ;
		// Check if current value is at its right position
		if (val <= liste[f].val) break;
		// Swap the current value and its parent if needed
		tmp_val = liste[f].val;
		liste[f].val = val;
		liste[nb].val = tmp_val;
		tmp_i = liste[f].i;
		liste[f].i = i;
		liste[nb].i = tmp_i;
		nb = f; // Redo
	}
}

/* Remove root (ie the greatest value) and rebuild the heap
   In order to rebuild the heap, one of the lowest value in the heap is put at
   the root; then, the processing loop put it back at its right position by
   swaping this value with its child as long as needed. The last value of the
   heap is one of the lowest value in the heap.
   liste : pointer to the heap
   nb : the current heap size
   rv : pointer to retrieve the root
*/
void remove_root(RANG_VAL *liste, int nb, RANG_VAL *rv) {
	int f = 0, tmp_i, max_f;
	double tmp_val;

	// Retrieve the root if a pointer is given
	if (rv != NULL) {
		rv->val = liste[0].val;
		rv->i = liste[0].i;
	}
	// Decrease heap size
	nb--;
	// Put the last value at the root
	liste[0].val = liste[nb].val;
	liste[0].i = liste[nb].i;
	// Loop to find the right position in the heap (from top to bottom)
	while (f < nb) {
		// Get the first child position in the heap
		max_f = 2 * f + 1;
		// Stop if no child
		if (max_f >= nb) break;
		// Get the position of the child with the greatest value
		if ((max_f + 1) < nb)
			if (liste[max_f].val < liste[max_f + 1].val) max_f++;
		// Check if current value is at its right position
		if (liste[f].val >= liste[max_f].val) break;
		// Swap the current value and its parent if needed
		tmp_val = liste[f].val;
		liste[f].val = liste[max_f].val;
		liste[max_f].val = tmp_val;
		tmp_i = liste[f].i;
		liste[f].i = liste[max_f].i;
		liste[max_f].i = tmp_i;
		f = max_f; // Redo
	}
}

/* Replace root (ie the greatest value) and rebuild the heap
   In order to rebuild the heap, the new value is put at the root of the heap;
   then, the processing loop put it at its right position by swaping this value
   with its child as long as needed.
   liste : pointer to the heap
   nb : the current heap size
   val, i : the (value,i) to be inserted in the heap
*/
void replace_root(RANG_VAL *liste, int nb, double val, int i) {
	int f = 0, tmp_i, max_f;
	double tmp_val;

	// Put the new value at the root
	liste[0].val = val;
	liste[0].i = i;
	// Loop to find the right position in the heap (from top to bottom)
	while (f < nb) {
		// Get the first child position in the heap
		max_f = 2 * f + 1;
		// Stop if no child
		if (max_f >= nb) break;
		// Get the position of the child with the greatest value
		if ((max_f + 1) < nb)
			if (liste[max_f].val < liste[max_f + 1].val) max_f++;
		// Check if current value is at its right position
		if (liste[f].val >= liste[max_f].val) break;
		// Swap the current value and its parent if needed
		tmp_val = liste[f].val;
		liste[f].val = liste[max_f].val;
		liste[max_f].val = tmp_val;
		tmp_i = liste[f].i;
		liste[f].i = liste[max_f].i;
		liste[max_f].i = tmp_i;
		f = max_f; // Redo
	}
}

// Same functions to manage double

void add_double(double *liste, int nb, double val) {
	int f;
	double tmp_val;

	liste[nb] = val;
	while (nb != 0) {
		if (nb % 2)
			f = (nb - 1) / 2;
		else
			f = nb / 2 - 1;
		if (val > liste[f]) {
			tmp_val = liste[f];
			liste[f] = val;
			liste[nb] = tmp_val;
			nb = f;
		}
		else
			nb = 0;
	}
}

void remove_root_double(double *liste, int nb, double *ret) {
	int f = 0, max_f;
	double tmp_val;

	if (ret != NULL) *ret = liste[0];
	nb--;
	liste[0] = liste[nb];
	while (f < nb) {
		max_f = 2 * f + 1;
		if ((max_f + 1) < nb)
			if (liste[max_f] < liste[max_f + 1]) max_f++;
		if (max_f < nb)
			if (liste[f] < liste[max_f]) {
				tmp_val = liste[f];
				liste[f] = liste[max_f];
				liste[max_f] = tmp_val;
				f = max_f;
			}
			else
				f = nb;
		else
			f = nb;
	}
}

void replace_double(double *liste, int nb, double val) {
	int f = 0, max_f;
	double tmp_val;

	// Put the new value at the root
	liste[0] = val;
	// Loop to find the right position in the heap (from top to bottom)
	while (f < nb) {
		// Get the first child position in the heap
		max_f = 2 * f + 1;
		// Stop if no child
		if (max_f >= nb) break;
		// Get the position of the child with the greatest value
		if ((max_f + 1) < nb)
			if (liste[max_f] < liste[max_f + 1]) max_f++;
		// Check if current value is at its right position
		if (liste[f] >= liste[max_f]) break;
		// Swap the current value and its parent if needed
		tmp_val = liste[f];
		liste[f] = liste[max_f];
		liste[max_f] = tmp_val;
		f = max_f; // Redo
	}
}

/*==========================================================================*/
/* Rank transformation
*/

// Compares values for the RANG_VAL structure (used by qsort)
static int compRangVal(const void *v1 ,const void *v2) {
	RANG_VAL *vp1 = (RANG_VAL*)v1;
	RANG_VAL *vp2 = (RANG_VAL*)v2;
	if (vp1->val > vp2->val)
		return 1;
 	else if (vp1->val < vp2->val)
		return -1;
	return 0;
}

// Rank transformation (ranks start at 1)
double *rank_transform(double *data, int nb) {
	double vc;
	int i, n, s, i0, k;
	double avgRank;
	double *m = (double *)calloc(nb, sizeof(double));
	RANG_VAL *vp = (RANG_VAL *)calloc(nb, sizeof(RANG_VAL));

	// Copy value and rank
	for (i = 0; i < nb; i++) {
		vp[i].i = i + 1;
		vp[i].val = data[i];
	}
	// Sort by value along with rank
	qsort(vp, nb, sizeof(*vp), compRangVal);
	// Compute rank and copy to new array
	i = 0;
	i0 = i, n = 1, s = i + 1;
	vc = vp[i].val;
	i++;
	while (i < nb) {
		// Account for duplicate values
		if (vp[i].val == vc) {
			s += i + 1;
			n++;
		}
		else {
			// Set average rank
			avgRank = (double)s / (double)n;
			for (k = i0; k < i; k++)
				m[vp[k].i-1] = avgRank;
			i0 = i, n = 1, s = i + 1;
			vc=vp[i].val;
		}
		i++;
	}
	avgRank = (double)s / (double)n;
	for (k = i0; k < i; k++)
		m[vp[k].i - 1] = avgRank;
	// Release memory and return rank values
	free(vp);
	return m;
}

// Gene name comparison for extract heart
int compare(const void *pa, const void *pb) {
	return strcmp((char *)pa, (char *)pb);
}

// Log time spent
void logClock(clock_t clock1, clock_t clock0, int verbose) {
	if (verbose>1)
		fprintf(stderr,"Time = %7.2f\n",(clock1-clock0)/(double)CLOCKS_PER_SEC);
}

/*==========================================================================*/
/* Extra output */

void compute_core_genes( double *obs_dist, int data_nr, double cut_off,
	RANG_VAL **rv_triee, int K, int *nb, int **indices )
/* Compute the number of genes at selected cut-off
   Take into account that observed distance could appear more than once,
   so the number of kept genes must be exactly computed
*/
{
	int i, k, *p;
	// Compute the number of genes at the selected cut-off
	for (k = i = 0; i < data_nr; i++)
		//if (obs_dist[i] < cut_off) k2++ else break;
		if (rv_triee[i][K - 1].val < cut_off) k++;
	// Return values
	*nb = k;
	*indices = p = (int *)malloc(k * sizeof(**indices));
	// Fill the array of gene indices (NB sorted by construction)
	for (i = 0; i < data_nr; i++)
		if (rv_triee[i][K - 1].val < cut_off) *p++ = i;
}

void compute_extra_genes( double *obs_dist, int data_nr, double cut_off,
	RANG_VAL **rv_triee, int K, int L, int *nb, int **indices )
/* Compute the number of genes at selected cut-off
   Take into account that observed distance could appear more than once,
   so the number of kept genes must be exactly computed
*/
{
	int i, k, *p;
	// Compute the number of genes at the selected cut-off
	for (k = i = 0; i < data_nr; i++)
		//if (obs_dist[i] < cut_off) k2++ else break;
		if (rv_triee[i][K - L - 1].val < cut_off &&
			rv_triee[i][K - 1].val >= cut_off) k++;
	// Return values
	*nb = k;
	*indices = p = (int *)malloc(k * sizeof(**indices));
	// Fill the array of gene indices (NB sorted by construction)
	for (i = 0; i < data_nr; i++)
		if (rv_triee[i][K - L - 1].val < cut_off &&
			rv_triee[i][K - 1].val >= cut_off) *p++ = i;
}


/*==========================================================================*/
/* Extra output
   All outputs are written in a single file. Sections are separated by the
   following FORMAT_OPTION
*/

#define FORMAT_OPTION ">>%s\n"

void fprint_genes( FILE *fd, char *option, int nb, int *indices, char **names )
/* Print the list of gene names aka probe ids along with corresponding indices
*/
{
	int i;
	fprintf( fd, FORMAT_OPTION, option );
	for ( i=0 ; i<nb ; i++ )
		fprintf( fd, "%d\t%s\n", indices[i], names[indices[i]] );
}

void fprint_value( FILE *fd, char *option, int nb, int *indices, char **names,
	int data_nc, double **data )
/* Print the data for the selected list of genes
   data : raw, rank or gaussian like
   genes : list of selected genes
*/
{
	int i, j;
	fprintf( fd, FORMAT_OPTION, option );
	for ( i=0 ; i<nb ; i++ ) {
		fprintf( fd, "%s\t", names[indices[i]] );
		for ( j=0 ; j<data_nc ; j++ )
			fprintf( fd, j?"\t%f":"%f", data[indices[i]][j] );
		fprintf( fd, "\n" );
	}
}

void fprint_dists( FILE *fd, char *option, int data_nr, double *obs_dist,
	int nb_rdn, double **rand_dist )
/* Prints all the distances
*/
{
	int i, j;
	fprintf( fd, FORMAT_OPTION, option );
	for ( i=0 ; i<data_nr ; i++ ) {
		fprintf( fd, "%f", obs_dist[i] );
		for ( j=0 ; j<nb_rdn ; j++ )
			fprintf( fd, "\t%f", rand_dist[j][i] );
		fprintf( fd, "\n" );
	}
}

void fprint_knn_dists( FILE *fd, char *option, int data_nr, char **names,
	int K, RANG_VAL **rv_triee)
/* Prints the sorted distances per gene
*/
{
	int i, j;
	fprintf( fd, FORMAT_OPTION, option );
	for ( i=0 ; i<data_nr ; i++ ) {
		fprintf( fd, "%s", names[i] );
		for ( j=0 ; j<K ; j+=10 )
			fprintf( fd, "\t%f", rv_triee[j][i].val );
		fprintf( fd, "\n" );
	}
}

void fprint_stat_dists( FILE *fd, char *option )
/* Print the deciles of the distance array
   dists : array of distances
   To be coded
*/
{
	fprintf( fd, FORMAT_OPTION, option );
}

void fprint_thresholds( FILE *fd, char *option,
	double cut_off, int kept_genes, double fdr )
/* Print cut-off, kept genes and fdr
*/
{
	fprintf( fd, FORMAT_OPTION, option );
	fprintf( fd, "%f\n", cut_off );
	fprintf( fd, "%d\n", kept_genes );
	fprintf( fd, "%f\n", fdr );
}

void fprint_selected( char *options, char *fr, int verbose,
	int data_nr, int data_nc, char **genes, double **mat,
	double *obs_dist, int nb_rdn, double **rand_dist,
	int K, RANG_VAL **rv_triee,
	double cut_off, int kept_genes, double fdr,
	int core_nb, int *core_indices, int extra_nb, int *extra_indices
	)
/* Print out in a single file selected objects
   The options field defines the objects that will be printed out and the
   order. Options is splitted and the corresponding function is called.
   The output file is the assembly of DBF input file and an extra token.
   Parameters are the same as the parameters or objects in the DBF function.
*/
{
	char *delims = ",";
	char *option = NULL;
	char *outputFN = NULL;
	FILE *fd;
	char len = 0;

	// Output file naming and opening
	len = strlen( fr ) + 30;
	outputFN = (char *)malloc( (len+1) * sizeof(* outputFN));
	strncpy( outputFN, fr, len );
	fd = fopen( strncat( fr, "-dbfAll.txt", 30 ), "w" );
	if (fd == NULL) {
		fprintf(stderr, "File %s can\'t be written!!!\n", outputFN);
		return;
	}

	// Format parser and output
	option = strtok( options, delims );
	while ( option != NULL ) {
		if (!strcmp( option, "coreList" ))
			fprint_genes( fd, option, core_nb, core_indices, genes );
		else if (!strcmp( option, "coreRawVal" ))
			fprint_value( fd, option, core_nb, core_indices, genes, data_nc, mat );
		else if (!strcmp( option, "extraList" ))
			fprint_genes( fd, option, extra_nb, extra_indices, genes );
		else if (!strcmp( option, "extraRawVal" ))
			fprint_value( fd, option, extra_nb, extra_indices, genes, data_nc, mat );
		else if (!strcmp( option, "dists" ))
			fprint_dists( fd, option, data_nr, obs_dist, nb_rdn, rand_dist );
		else if (!strcmp( option, "knnDists" ))
			fprint_knn_dists( fd, option, data_nr, genes, K, rv_triee );
		else if (!strcmp( option, "thresholds" ))
			fprint_thresholds( fd, option, cut_off, kept_genes, fdr );
		else
			if (verbose>1)
				fprintf( stderr, "Unknown format option \'%s\'.\n", option );
		option = strtok( NULL, delims );
	}

	// Closing
	fclose( fd );
	free( outputFN );
}

/*================================= BEGIN DBF ==============================*/
/* DBF
mat         :   a matrix of double values in C format
nb_genes	:	le nombre de genes (= de lignes dans la matrice data)
nb_samples	:	le nombre d'echantillons (= de colonnes dans la matrice data)
list_genes	:	la liste des sondes (vecteur de chaines de caracteres en R)
list_samples	:	la liste des echantillons (vecteur de chaines de caracteres en R)
dm		:	le type de distance ("spearman", "euclidean" ou "pearson")
K		:	le nombre de voisins pour chaque gene
nbr		:	le nombre de randomisations sur les distances
verb		:	0 = silencieux ; 1 = affichage des progressions
memory	        :	la quantite de memoire, en Mo, allouee au sous-echantillon de distances
FDR		:	le % de faux positifs accepte
p_dist		:	1 = sauvegarde les distributions (observee et random) ; 0 = pas de sauvegarde
out		:	tableau de chaines de caracteres (selon R) pour stocker les genes conserves par DBF
gf    		:	nom du fichier de sortie pour le graphe (pour MCL)
set_seed        :       nombre entier pour positionner le generateur de nombres aléatoire */

void dbf_(double **mat, int *nb_genes, int *nb_samples,
		char **list_genes, char **list_samples,
		char **dm, int *knn, int *nbr, int *verb, int *memory, int *FDR,
		int *p_dist, char **out, char **gf, int *set_seed,
		int *lnn, char **out_format, char **file_root )
{

	// On dereference les pointeurs envoyes par R
	int data_nr = *nb_genes;
	int data_nc = *nb_samples;
	int K = *knn;
	int nb_rdn = *nbr;
	int verbose = *verb;
	int mem = *memory;
	int fdr = *FDR;
	char *dist_m = *dm;
	char *gr_file = *gf;
	int seed = *set_seed;

	// Extended options
	int L = *lnn;
	char *options = *out_format;
	char *fr = *file_root;

	// des variables bien utiles ...
	int i, j, k, q;

	// utilise pour compter les valeurs dans le sous-echantillon de distances random
	int N = 0;

	// distances min et max ; utilises pour le calcul de la correlation (0..1) dans le cas des distances euclidiennes
	int min, max;

	// le nombre de genes conserves pour le FDR demande
	int kept_genes = 0;

	// variables utilisees pour stocker le pourcentage d'avancement (verbose)
	double pc, pc0;

	// variables utilisees pour le calcul des distances
	double dist = 0, sum;

	// utilise pour le stockage temporaire des distances random
	double rvl;

	// tableau utilise pour le calcul du pearson
	double *carre = NULL;

	// tableau utilise pour trier les distances random de chaque gene
	double **tri;

	// clocking
	clock_t clock0 = -1;

	#if debugLevel > 2
	fprintf(stderr,"data_nr = %d\n",*nb_genes);
	fprintf(stderr,"data_nc = %d\n",*nb_samples);
	fprintf(stderr,"K = %d\n",*knn);
	fprintf(stderr,"nb_rdn = %d\n",*nbr);
	fprintf(stderr,"verbose = %d\n",*verb);
	fprintf(stderr,"mem = %d\n",*memory);
	fprintf(stderr,"fdr = %d\n",*FDR);
	fprintf(stderr,"dist_m = %s\n",*dm);
	fprintf(stderr,"gr_file = %s\n",*gf);
	fprintf(stderr,"seed = %d\n",*set_seed);
	fprintf(stderr,"gr_file = %s\n",*gf);
	fprintf(stderr,"out_format = %s\n",*out_format);
	fprintf(stderr,"file_root = %s\n",*file_root);
	#endif

	// matrice utilisee pour la transformation en rangs sur les lignes
	double **mat_rank = (double **)malloc(data_nr * sizeof(double *));

	// Calcul du nombre de distances dans le sous-echantillon
	int nb_randomval = mem * 1024 * (1024 / sizeof(dist));

	double max_randomval = ((data_nr - 2.0) * (data_nr - 1.0)) / 2;
	if (max_randomval < nb_randomval) nb_randomval = max_randomval;

	// Calcul du taux de sous-echantillonnage (x 1000)
	unsigned int NN = (double)(max_randomval) / nb_randomval * 1000.0;
	double rand_ratio = NN / (RAND_MAX + 1.0);

	// Reservation de la memoire pour le sous-echantillon
	double *randomval = (double *)malloc(nb_randomval * sizeof(double));

	if (verbose) {
		fprintf(stderr, "Randomization: %d ", nb_randomval);
		fprintf(stderr, "(1/%-8.3f ratio)\n", (double)(NN/1000.0));
	}

	// Demarrage du generateur de nombres aleatoires
	if (seed==0) {
	  time_t tt;
	  time(&tt);
	  srand(tt);
	  if (verbose) fprintf(stderr, "Seed = %ld\n", tt);
	} else {
	  srand(seed);
	  if (verbose) fprintf(stderr, "Seed = %d\n", seed);
	}

	// Calcul des distances

	// Affichage de la progression
	if (verbose) {
		fprintf(stderr, "Pre-computation for distances\n");
		clock0 = clock();
	}

	// Pearson pre-computing
	if (!strcmp("pearson", dist_m) || !strcmp("spm", dist_m) || !strcmp("spgm", dist_m)) {
		// Each gene is mean-centered and sum of squarre is computed
		double moy;
		carre = (double *)calloc(data_nr, sizeof(double));
		for (i = 0; i < data_nr; i++) {
			moy = 0.0;
			for(j = 0; j < data_nc; j++) moy += mat[i][j];
			moy /= data_nc;
			for (j = 0; j < data_nc; j++) {
				mat[i][j] -= moy;
				carre[i] += mat[i][j] * mat[i][j];
			}
		}
	}

	// Spearman pre-computing
	double spearman_dist_coeff = 6.0 / (data_nc * ((data_nc * data_nc) - 1.0));
	// Convert intensities to ranks
	if (!strcmp("spearman", dist_m) || !strcmp("spm", dist_m) || !strcmp("spgm", dist_m)) // transformation par les rangs sur chaque ligne ...
		for (i = 0; i < data_nr; i++) mat_rank[i] = rank_transform(mat[i], data_nc);

	// Euclidean pre-computing
	// Init min, max
	min = 5000000;
	max = 0;

	// Affichage de la progression
	if (verbose) {
		logClock( clock(), clock0, verbose );
		pc = pc0 = 0.0;
		fprintf(stderr, "Computing distances: %6.2f%%\r", pc);
		clock0 = clock();
	}

	// Reservation du tableau de RANG_VAL
	// Chaque ligne contiendra K valeurs et sera triee par ordre croissant de rang
	RANG_VAL ret, **rv_triee = (RANG_VAL **)malloc(data_nr * sizeof(RANG_VAL *));
	for (i = 0; i < data_nr; i++) rv_triee[i] = (RANG_VAL *)calloc(K, sizeof(RANG_VAL));

	// Compute distance between gene pairs
	double dist_diff = 0;
	for (i = 0; i < data_nr; i++) {
		for (j = i + 1; j < data_nr; j++) {

			// Compute distance
			if (!strcmp("spearman", dist_m)) {      // SPEARMAN  0 <= dist <= 2
				sum = 0;
				for (k = 0; k < data_nc; k++) { dist_diff = mat_rank[i][k] - mat_rank[j][k] ; sum += dist_diff * dist_diff; };
				dist = sum * spearman_dist_coeff;
			}
			else if (!strcmp("pearson", dist_m)) {  // PEARSON   0 <= dist <= 2
				sum = 0;
				for (k = 0; k < data_nc; k++) sum += mat[i][k] * mat[j][k];
				dist = 1 - sum / sqrt(carre[i] * carre[j]);
			}
			else if (!strcmp("euclidean", dist_m)) { // EUCLID
				sum = 0;
				for (k = 0; k < data_nc; k++) sum += (mat[i][k] - mat[j][k]) * (mat[i][k] - mat[j][k]);
				dist = sqrt(sum);
				// Update "min" et "max"
				if (dist > max) max = dist;
				if (dist < min) min = dist;
			}
			else if (!strcmp("spm", dist_m) || !strcmp("spgm", dist_m)){ // PEARMAN == mean (SPEARMAN, PEARSON)!!!!
			  //spearman
			  sum = 0;
			  for (k = 0; k < data_nc; k++) { dist_diff = mat_rank[i][k] - mat_rank[j][k] ; sum += dist_diff * dist_diff; };
			  dist = sum * spearman_dist_coeff;
			  //pearson
			  sum = 0;
			  for (k = 0; k < data_nc; k++) sum += mat[i][k] * mat[j][k];
			  if(!strcmp("spm", dist_m)){
			  //calcul de la moyenne arithmetique [ (s+p)/2 ]
			    dist += 1 - sum / sqrt(carre[i] * carre[j]);
			    dist /= 2.0;
			  }
			  else{
			  //calcul de la moyenne geometrique [ sqrt(s*p) ]
			    dist = sqrt(dist*(1 - sum / sqrt(carre[i] * carre[j])));
			  }
			}

			// Add the distance to the heap of gene i
			if ((j - 1) < K)
				// If heap is not full with K distances, then add
				add_value(rv_triee[i], j - 1, dist, j);
			else if (dist < rv_triee[i][0].val) {
				// If heap is already filled with K distances,
				// and the current distance is lower than the greatest value,
				// then remove root (greatest value) and add
				//remove_root(rv_triee[i], K, NULL);
				//add_value(rv_triee[i], K - 1, dist, j);
				replace_root(rv_triee[i], K, dist, j);
			}

			// Symmetrically, put the distance in the heap of gene j
			if (i < K)
				add_value(rv_triee[j], i, dist, i);
			else if (dist < rv_triee[j][0].val) {
				//remove_root(rv_triee[j], K, NULL);
				//add_value(rv_triee[j], K - 1, dist, i);
				replace_root(rv_triee[j], K, dist, i);
			}

			// Randomly put the distance into the sample
			if (N < nb_randomval) {
				k = (int)(rand() * rand_ratio);
				if (k < 1000) randomval[N++] = dist;
			}
		}

#define KEEP_ALL_K_SORTED_DISTANCES
		/* undefining this constant remove the sort of the array, but is
		incompatible with the lnn option */
#ifdef KEEP_ALL_K_SORTED_DISTANCES
		/* The heap is transformed in a sorted list; the root, is the greatest
		   value, is removed and the heap is rebuild */
		for (j = 0; j < K; j++) {
			remove_root(rv_triee[i], K - j, &ret);
			rv_triee[i][K - j - 1].val = ret.val;
			rv_triee[i][K - j - 1].i = ret.i;
		}
#else
		/* Focuss only on the greatest value, ie Kth value */
		remove_root(rv_triee[i], K, &ret);
		rv_triee[i][K - 1].val = ret.val;
#endif

		// Affichage de la progression
		if (verbose) {
			pc = (int)(i * 1000.0 / data_nr) / 10.0;
			if (pc != pc0) {
				fprintf(stderr, "Computing distances: %6.2f%%\r", pc);
				pc0 = pc;
			}
		}
	}

	// Fin du calcul des distances
	if (verbose) {
		fprintf(stderr, "Computing distances: 100.00%%\n");
		logClock( clock(), clock0, verbose );
		fprintf(stderr, "Randomization: %d obtained, %d asked\n", N, nb_randomval);
	}

	// Calcul du cut-off

	// Affichage de la progression
	if (verbose) {
		pc = pc0 = 0.0;
		fprintf(stderr, "Computing FDR: %6.2f%%\r", pc);
		clock0 = clock();
	}


	/* Extract the observed distance distribution
	   The observed distance vector consists in the Kth distance calculated for
	   each gene; this vector is sorted from lowest to greatest.
	*/
	double *obs_dist = (double *)malloc(data_nr * sizeof(double));
	for (i = 0; i < data_nr; i++) {   // Extract
		add_double(obs_dist, i, rv_triee[i][K - 1].val);
	}
	double tmp;
	for (j = 0; j < data_nr; j++) {   // Sort
		remove_root_double(obs_dist, data_nr - j, &tmp);
		obs_dist[data_nr - j - 1] = tmp;
	}

	/* Compute the simulated distribution from randomly selected distance
	   The simulated set consist in nb_rdn simulated sets. Each set consists in
	   data_nr distances, simulating one for each gene. Representative distance
	   for each gene is aimed to mimic the observed distance calculation.
	   The loops are similar to previous, except that only distance are needed,
	   not position.
	*/
	double **rand_dist = (double **)malloc(nb_rdn * sizeof(double *));
	for (i = 0; i < nb_rdn; i++)
		rand_dist[i] = (double *)malloc(data_nr * sizeof(double));

	// Setup an array of K distance for each gene (similar to rv_triee)
	tri = (double **)malloc(data_nr * sizeof(double *));
	for (i = 0; i < data_nr; i++)
		tri[i] = (double *)calloc(K, sizeof(double));

	rand_ratio = N / (RAND_MAX + 1.0);
	for (q = 0; q < nb_rdn; q++) {
		for (i = 0; i < data_nr; i++) {
			for (j = i + 1; j < data_nr; j++) {

				// The distance is randomly selected from the sample
				k = (int)(rand() * rand_ratio);
				rvl = randomval[k];

				// The distance is added to the heap of gene i
				if ((j - 1) < K)
					add_double(tri[i], j - 1, rvl);
				else if (rvl < tri[i][0]) {
					//remove_root_double(tri[i], K, NULL);
					//add_double(tri[i], K - 1, rvl);
					replace_double(tri[i], K, rvl);
				}

				// Symmetrically, put the distance in the heap of gene j
				if (i < K)
					add_double(tri[j], i, rvl);
				else if (rvl < tri[j][0]) {
					//remove_root_double(tri[j], K, NULL);
					//add_double(tri[j], K - 1, rvl);
					replace_double(tri[j], K, rvl);
				}
			}
#ifdef KEEP_ALL_K_SORTED_DISTANCES
			for (j = 0; j < K; j++) {
				remove_root_double(tri[i], K - j, &tmp);
				tri[i][K - j - 1] = tmp;
			}
#else
			remove_root_double(tri[i], K, &tmp);
			tri[i][K - 1] = tmp;
#endif

			// Affichage de la progression
			if (verbose) {
				pc = (int)((i + q * data_nr) * 1000.0 / (data_nr * nb_rdn)) / 10.0;
				if (pc != pc0) {
					fprintf(stderr, "Computing FDR: %6.2f%%\r", pc);
					pc0 = pc;
				}
			}
		}
		/* Extract the random distance distribution */
		for (i = 0; i < data_nr; i++) {   // Extract
			add_double(rand_dist[q], i, tri[i][K - 1]);
		}
		for (j = 0; j < data_nr; j++) {   // Sort
			remove_root_double(rand_dist[q], data_nr - j, &tmp);
			rand_dist[q][data_nr - j - 1] = tmp;
		}
	}
	if (verbose) {
		fprintf(stderr, "Computing FDR: 100.00%%\n");
		logClock( clock(), clock0, verbose );
		fprintf(stderr, "Computing cut-off\n");
		clock0 = clock();
	}

	/* Compute cut-off at selected FDR
	   Observed distribution is compared to distributions obtained by
	   randomization. All distributions are sorted by increasing distance.
	   For each observed distance, distances from random ditrib. up to the same
	   index that are lower than current observed distance are counted. Dividing
	   this count by the number of randomization leads a FDR value.
	   NB: performance could be improved by a double index loop, one index
	   running on obs_dit, the other one on rand_dist.
	*/
	double cnt_false;
	for (i = 0; i < data_nr; i++) {
		cnt_false = 0.0;
		for (j = 0; j < nb_rdn; j++) {
			q = 0;
			for (k = 0; k <= i; k++)
				if (rand_dist[j][k] <= obs_dist[i]) q++;
			cnt_false += q;
		}
		cnt_false /= (double)nb_rdn; // average to nb of random sets
		if (fdr < (cnt_false / (i + 1) * 100)) {
			kept_genes = i; // nb of genes at setup FDR
			break;
		}
	}
	double cut_off = obs_dist[kept_genes];
	if (verbose) {
	  //fprintf(stderr, "cut-off = %f genes = %d\n", cut_off, kept_genes);
		fprintf(stderr, "number of conserved genes = %d\n", kept_genes);
	}

	/* Compute the number of genes at selected cut-off
	   Take into account that observed distance could appear more than once,
	   so the number of kept genes must be exactly computed
	*/
	int k2 = 0;
	for (i = 0; i < data_nr; i++)
		//if (obs_dist[i] < cut_off) k2++ else break;
		if (rv_triee[i][K - 1].val < cut_off) k2++; // same test as the following loops
	k2 *= 1 + K;

	/* Build up a graph and a list of kept genes
	   NB: graph seems not to be oriented for MCL; auto-links do not change
	   MCL output, so could be safely removed. Duplicate links should not
	   influence MCL, but in fact it does; duplicate links could not be
	   replaced by a link whose distance is half the real distance.
	*/
	if (verbose) {
		fprintf(stderr, "Building graph\n");
	}

	char ***fwgtgraph = (char ***)malloc(k2 * sizeof(char **));
	for (i = 0; i < k2; i++) fwgtgraph[i] = (char **)malloc(2 * sizeof(char *));

	double *weight = (double *)malloc((K * data_nr) * sizeof(double));

	k = k2 = 0; // k total ap boucle = kept_genes
	for (i = 0; i < data_nr; i++)

		// If distance is lower than cut-off, gene is kept
		if (rv_triee[i][K - 1].val < cut_off) {
			out[k] = list_genes[i]; // gene name is kept
			k++;
			// link the gene to itself
			fwgtgraph[k2][0] = list_genes[i];
			fwgtgraph[k2][1] = list_genes[i];
			weight[k2] = 1;
			k2++;
			// link the gene to each its neighbours
			for (j = 0; j < K; j++) {
				fwgtgraph[k2][0] = list_genes[i];
				fwgtgraph[k2][1] = list_genes[rv_triee[i][j].i];
				if (strcmp("euclidean", dist_m))
					weight[k2] = 1 - rv_triee[i][j].val;
				else if(rv_triee[i][j].val != 0)
					weight[k2] = (rv_triee[i][j].val - min) / (max - min);
				else
					weight[k2] = 0;
				k2++;
			}
		}

	/* Extract hearts, ie filter the graph
	   Links are printed out if both nodes have passed the cut-off criteria
	*/
	/*	if (verbose) {
		fprintf(stderr, "Extracting hearts\n");
		}*/

	void *racine = NULL;
	void *val;
	FILE *sortie = fopen(gr_file, "w");

	// Build the binary tree with all kept genes
	// NB: tsearch adds a value if it can't find it
	for(i = 0; i < k; i++) {
		val = tsearch((void *)out[i], &racine, compare);
		if (val == NULL) exit(EXIT_FAILURE);
	}

	// Output a filtered graph
	for (i = 0; i < k2; i++)
		if (tfind((void *)fwgtgraph[i][0], &racine, compare) != NULL)
			if (tfind((void *)fwgtgraph[i][1], &racine, compare) != NULL)
				fprintf(sortie, "%s %s %g\n", fwgtgraph[i][0], fwgtgraph[i][1], weight[i]);
	fclose(sortie);


	/* */
	int core_nb = 0;
	int *core_indices = NULL;
	int extra_nb = 0;
	int *extra_indices = NULL;
	compute_core_genes( obs_dist, data_nr, cut_off, rv_triee, K,
		&core_nb, &core_indices );
	compute_extra_genes( obs_dist, data_nr, cut_off, rv_triee, K,
		L, &extra_nb, &extra_indices );
	if (verbose) {
		fprintf(stderr, "Genes   core = %d   extra = %d\n", core_nb, extra_nb);
	}


	/* Extra outputs
	*/
	if (verbose) {
		logClock( clock(), clock0, verbose );
		//fprintf(stderr, "Printing out\n");
		clock0 = clock();
	}
	if (*options)
	fprint_selected( options, fr, verbose,
		data_nr, data_nc, list_genes, mat,
		obs_dist, nb_rdn, rand_dist,
		K, rv_triee,
		cut_off, kept_genes, fdr,
		core_nb, core_indices, extra_nb, extra_indices
	);

	/* Free
	*/
	if (verbose) {
		logClock( clock(), clock0, verbose );
		//fprintf(stderr, "Freeing\n");
	}
	for (i = 0; i < data_nr; i++) {
	  free(mat[i]);
	  if (!strcmp("spearman", dist_m)||!strcmp("spm", dist_m)||!strcmp("spgm", dist_m)) {
	    free(mat_rank[i]);
	  }
	  free(rv_triee[i]);
	  free(tri[i]);
	}
	free(mat_rank);
	free(randomval);
	if (carre != NULL) free(carre);
	free(rv_triee);
	free(obs_dist);
	for (i = 0; i < nb_rdn; i++) free(rand_dist[i]);
	free(rand_dist);
	free(tri);
	for (i = 0; i < k2; i++) free(fwgtgraph[i]);
	free(fwgtgraph);
	free(weight);
	if (verbose) {
		fprintf(stderr, "DBF done\n");
       	}

}

/*
data		:	une matrice de double, selon R
*/
void DBF(double *data, int *nb_genes, int *nb_samples,
	char **list_genes, char **list_samples,
	char **dm, int *knn, int *nbr, int *verb, int *memory, int *FDR,
	int *p_dist, char **out, char **gf, int *set_seed)
{
	int lnn = 0;
	char *out_format = "";
	char *fr = "";

	// Convert R matrix in C format
	double **mat = read_double_matrix( *nb_genes, *nb_samples, data );
	// Real process
	dbf_( mat, nb_genes, nb_samples,
		list_genes, list_samples,
		dm, knn, nbr, verb, memory, FDR,
		p_dist, out, gf, set_seed,
		&lnn, &out_format, &fr );
	//
	free(mat);
}

/*================================ END DBF =================================*/

/*==========================================================================*/
/* Command-line helpers
*/

/* Fonction qui split une ligne en fonction d'un separateur.
   Renvoie le nombre de champs et alloue automatiquement la place
   necessaire pour les valeurs du tableau en memoire
*/
int split(char ***tab, char *s, char *delim) {
	char *tmp = strdup(s);
	int i, n, k, in_token;
	in_token = n = k = 0;

	for (i = 0; i < strlen(s); i++)
		if (strchr(delim, (int)(*(tmp + i))) != NULL) {
			*(tmp + i) = 0;
			in_token = 0;
		}
		else if (!in_token) {
			in_token = 1;
			n++;
		}
	*tab = (char **)calloc(n, sizeof(char *));
	for (i = 0; i < strlen(s); i++)
		if (*(tmp + i) != 0) {
			(*tab)[k++] = strdup(tmp + i);
			i += strlen(tmp + i);
		}
	free(tmp);
	return n;
}

#define BUFF_LEN 19999

/* Read data from file
*/
void read_file(const char *inputFN, int *nbGene, int *nbSample,
	char ***geneName, char ***sampleName, double ***dataVal)
{
	FILE *inputFH;
	char *buffer = (char *)malloc((BUFF_LEN + 1) * sizeof(char));
	int nb_sample, nb_gene;
	char **sample, **gene, **tab;
	double **data;
	int i, j;

	// Open data file
	inputFH = fopen(inputFN, "ro");
	if (inputFH == NULL) {
		fprintf(stderr, "File %s not found!!!\n", inputFN);
		exit(-1);
	}

	// Read first line and split it to get sample name
	fgets(buffer, BUFF_LEN, inputFH);
   	nb_sample = split(&sample, buffer, "\t") -1 ; // skip first column

	// Compute number of genes
	nb_gene = 0;
	while (fgets(buffer, BUFF_LEN, inputFH) != NULL)
		nb_gene++;

	// Allocate memory for gene names
	gene = (char **)malloc(nb_gene * sizeof(char *));

	// Allocate memory for data
	data = (double **)malloc(nb_gene * sizeof(double *));
	for (i = 0 ; i < nb_gene; i++)
		data[i] = (double *)malloc(nb_sample * sizeof(double));

	// Load data (rank or intensities)
	fseek(inputFH, 0, SEEK_SET); // reset file cursor
	fgets(buffer, BUFF_LEN, inputFH); // skip first line
	j = 0;
	while (fgets(buffer, BUFF_LEN, inputFH) != NULL) {
		split(&tab, buffer, "\t");
		gene[j] = strdup(tab[0]);
		for (i = 1; i <= nb_sample; i++)
			data[j][i - 1] = atof(tab[i]);
		for (i = 0; i <= nb_sample; i++)
			free(tab[i]);
		free(tab);
		j++;
	}

	// Done with loading
	fclose(inputFH);

	// Returned parameters
	*nbSample = nb_sample;
	*nbGene = nb_gene;
	*sampleName = sample;
	*geneName = gene;
	*dataVal = data;
}

// Convert C float read from data file to double matrix in R spirit
double *Cfloat2Rdouble(int nl, int nc, float **mat) {
	double *m = NULL;
	int l, c;

	// Memory allocation
	m = malloc( nl * nc * sizeof(* m) );

	// Convert
	for (l = 0; l < nl; l++)
		for (c = 0; c < nc; c++)
			*(m + c * nl + l) = mat[l][c];

	return m;
}

// Test functions
#define T_RAND_NB	20000
#define T_ARRAY_LENGTH 200
void test_replace_root() {
	// Define variables
	double distance[T_ARRAY_LENGTH];
	double rvl;
	int j;
	double rand_ratio = (RAND_MAX + 1.0);

	// Fill array with random numbers
	for (j=0 ; j<T_ARRAY_LENGTH ; j++) {
		rvl = rand() / rand_ratio;
		// Add random numbers
		if (j<T_ARRAY_LENGTH)
			add_double(distance, j, rvl);
		// Replace greatest value
		else if (rvl < distance[0]) {
			//remove_root_double(distance, T_ARRAY_LENGTH, NULL);
			//add_double(distance, T_ARRAY_LENGTH - 1, rvl);
			replace_double(distance, T_ARRAY_LENGTH, rvl);
		}
	}
	// Extract root recursively and check that the array is sorted
}

/*==========================================================================*/
