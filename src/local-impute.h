#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <getopt.h>
#include <string.h>
#include <math.h>

#define MAX_ID_LEN (1027)
#define MAX_PAIRS (2056)
#define DEF_P_ERR (1)
#define MAX_ALLELE_COUNTS (20)
#define ERROR_PER_THOUSAND (5)
#define XPAR1ST (60001) /* Start and stop of PAR on sec chrs */
#define XPAR1EN (2699520)
#define XPAR2ST (154931044)
#define XPAR2EN (155260560)
#define YPARST (59034050)
#define YPAREN (59363566)
#define MAX_BLACKLIST_LEN (10000)

/* (c) 2019 - Ed Green
   Astrea Forensics
*/

typedef struct impute_pair {
  char chr[MAX_ID_LEN+1];
  unsigned int L_pos;
  char L_rsID[MAX_ID_LEN+1];
  char la0; // linked site reference allele
  char la1; // linked site alternate allele
  short unsigned int T0L0; // observations of T=0; L=0 in panel
  short unsigned int T0L1;
  short unsigned int T1L0;
  short unsigned int T1L1;
  short unsigned int L0;   // counts of L=0 in data
  short unsigned int L1;   // counts of L=1 in data
} Impute_Pair;

typedef struct impute_set {
  char chr[MAX_ID_LEN+1];
  unsigned int T_pos;
  char T_rsID[MAX_ID_LEN+1];
  short unsigned int T0; // counts of T=0 in panel
  short unsigned int T1; // counts of T=1 in panel
  short unsigned int num_pairs;
  /* By convention impute_pairs[0] is the Target Site itself.
     It is in perfect LD with itself, of course */
  Impute_Pair impute_pairs[MAX_PAIRS];
  double epsilon;
} Impute_Set;

typedef struct blacklist {
  char** bl_ids;
  size_t num_bl_ids;
  char bl_ids_array[ MAX_BLACKLIST_LEN * MAX_ID_LEN ];
} Blacklist;

unsigned long long* fact; // array of factorial values

Blacklist* read_blacklist( const char* bl_fn );

Impute_Set* init_Impute_Set( void );

/* Add Target site data 
   Takes pointer to Impute_Set and data to populate
   the Target site
   Returns: 0 => everything copacetic
            non-zero => some problem
*/
int add_Target_to_Impute_Set( Impute_Set* isp,
			      char* chr,
			      unsigned int T_pos,
			      char* T_rsID,
			      short unsigned int t0,
			      short unsigned int t1 );

void clear_Impute_Set( Impute_Set* isp );

/* Takes pointer to an Impute_Set and information about
   a linked site.
   Adds this site to the impute_pairs[] iff
   there is room and iff
   this site has >= min_D linkage disequilibrium with
   the target site
   Returns:  0 => added it, everything copacetic
            -1 => no more room in impute_pairs
	    -2 => D < min_D
*/
int add_impute_pair( Impute_Set* isp,
		     char* chr,
		     unsigned int L_pos,
		     char* L_rsID,
		     char la0,
		     char la1,
		     short unsigned int t0l0,
		     short unsigned int t0l1,
		     short unsigned int t1l0,
		     short unsigned int t1l1,
		     double min_D,
		     int Target // boolean, true=> this is a target, add it regardless of D
		     );

/* Returns the LD value (D) given the
   input data in Impute_Pair* */
double calc_D( Impute_Pair* imp );

/* Add observed allele counts to the Linked Site at
   an Impute_Pair*
   Args: (1) Impute_Pair* imp => pointer to Impute_Pair 
         (2) short unsigned int l0 => counts of the 0 allele
                                      at the Linked site
   Args: (3) short unsigned int l1 => counts of the 1 allele
                                      at the Linked site
   Returns: nothing
*/
void add_linked_counts( Impute_Pair* imp,
			short unsigned int l0,
			short unsigned int l1 );

/* Returns the probability that the Target site genotype is
   0, 0 given the data at all Linked Sites in the input
   Impute_Set */
double p_allLgivenT00( Impute_Set* isp );
double p_allLgivenT01( Impute_Set* isp );
double p_allLgivenT11( Impute_Set* isp );


/* Find the probability of the observed alleles at site L given that
   the genotype at T is 0/0 and the haplotype frequencies in the panel 
   as described in the input Impute_Pair */
double p_LgivenT00( const Impute_Pair* imp, const double epsilon );
double p_LgivenT01( const Impute_Pair* imp, const double epsilon );
double p_LgivenT11( const Impute_Pair* imp, const double epsilon );

/* NchooseK 
   n choose k = n! / k!(n-k)!
*/
unsigned long long NchooseK( unsigned n, unsigned k,
			     unsigned long long* factorial );
unsigned long long* init_factorial( void );


/* These functions return the *PRIORS* for a particular genotype
   at the Target site, given the background populations frequencies
   of each allele. That is, these can be calculated without seeing
   any actual sequence data */

/* Returns the probability of the target genotype being 0,0 
   given the reference panel allele frequencies at T
 */
double p_T00( Impute_Set* isp );

/* Returns the probability of the target genotype being 0, 1
   i.e., heterozygous, given the reference panel allele
   frequencies at T => 2pq
   Hardy Freakin Weinberg
*/
double p_T01( Impute_Set* isp );

/* Returns the probability of the target genotype being 1,1
   given the reference panel allele frequencies at T
 */
double p_T11( Impute_Set* isp );

/* in_sex_haploid
   Determines if the Target site of an Impute_Set is in the
   special PAR (pseudo-autosomal region of a sex chromosome
   (X, Y, chrX, or chrY 
   Return 0 if not in PAR
   Return 1 if in PAR
*/
int in_sex_haploid( Impute_Set* isp );
  
/* isp2genotype 
   Determines the most likely genotype given data in the Impute Set
   Arguments: (1) Impute_Set* isp
              (2) double min_ratio the minimum odds ratio to call 
	          both alleles
	      (3) double min_ratio_zero_cov special ratio if no Target
                  site alleles are seen, usually more stringent
              (4) double min_ratio_target_only special ratio to call
                  genotype using only Target site data, usually more
                  stringent
              (5) double special_hom if non-zero, invoke special homozygous
                         call rule and use this ratio cutoff
              (6) char* a0 place to put the first allele call
              (7) char* a1 place to put the second allele call
	      (8) int use_priors
   Returns: nothing
   Populates the chars pointed to by a0 and a1 with the allele calls
*/
void isp2genotype( Impute_Set* isp, int Male, double min_ratio,
		   double min_ratio_zero_cov,
		   double min_ratio_target_only,
		   double special_hom, char* a0, char* a1,
		   int use_priors );

int unambiguous_genotype( const double pHomRef, const double pHet,
			  const double pHomAlt, const double min_ratio,
			  const char ref, const char alt,
			  char* a0, char* a1 );
int black_check( Blacklist* BL, const char* rsID );
static int cmpstringp(const void *p1, const void *p2);
FILE* fileOpen( const char* name, char access_mode[] );
int idCmp(const void* id1_, const void* id2_);
