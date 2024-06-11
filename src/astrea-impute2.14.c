#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include <pthread.h>
#include "hap.h"
#include "pileup.h"
#include "math.h"
#include "stats.h"
#include "time.h"

#define VERSION (14)
#define EPSILON (5) // error rate in 1000
#define EPSILON_D (6) // error rate in 1000 at C sites where alt is T
#define MAX_COV_DEF (24)
#define DOWN_COV_DEF (32) // if cov < MAX but higher than this, downsample
#define MQ_CUT (20)
#define BQ_CUT (20)
#define VERBOSE (0)
#define B_DEF (1)
#define BFC_DEFAULT (100) // Bayes Factor cutoff default
#define GLRC_DEFAULT (100) // Genotype likelihood ratio cutoff default
#define XPAR1ST (60001) /* Start and stop of PAR on sex chrs */
#define XPAR1EN (2699520)
#define XPAR2ST (154931044)
#define XPAR2EN (155260560)
#define YPAR1ST (10001)
#define YPAR1EN (2649520)
#define YPAR2ST (59034050)
#define YPAR2EN (59363566)
#define XPAR1ST_HG38 (10001) /* Start and stop of PAR on GRC38 sex chrs */
#define XPAR1EN_HG38 (2781479)
#define XPAR2ST_HG38 (155701383)
#define XPAR2EN_HG38 (156030895)
#define YPAR1ST_HG38 (10001)
#define YPAR1EN_HG38 (2781479)
#define YPAR2ST_HG38 (56887903)
#define YPAR2EN_HG38 (57217415)
#define MAX_HAP_BUFFER (512)

void help() {
  printf( "astrea-impute2 v%d\n", VERSION );
  printf( "(c) Astrea Forensics, Ed Green - 2019, 2020\n" );
  printf( "Generates genotype files from mpileup data and\n" );
  printf( "haplotype records.\n" );
  printf( "OUTPUT: SNP_ID\tchr\tpos\tAllele1\tAllele2\tcall_LR\tAlt_fr\tBayesFac\n" );
  printf( "OPTIONS:\n" );
  printf( "-H Input file of haplotypes\n" );
  printf( "-m Input file of mpileup data with base observations\n" );
  printf( "-M maximum coverage cutoff;            default = %d\n",
	  MAX_COV_DEF );
  printf( "   Ignores data at sites with observed coverage higher\n" );
  printf( "   than this value, if specified.\n" );
  printf( "-d downsample coverge to this;         default = %d\n",
	  DOWN_COV_DEF );
  printf( "   At sites whose observed coverage is under the cutoff\n" );
  printf( "   but higher than this number, down-sample the observed\n" );
  printf( "   data to this amount.\n" );
  printf( "-p percentile coverage cutoff;         default = none\n" );
  printf( "   If specified, calculate the empirical coverage of the\n" );
  printf( "   input mpileup data. Then, use all sites whose observed\n" );
  printf( "   coverage is less than this percentile. For example: \n" );
  printf( "   -p 95 would use all sites that are <= the observed\n" );
  printf( "   coverage at 95 percent of all sites in the mpileup data\n" );
  printf( "   Note that this option overrides the value of -M\n" );
  printf( "-e rate of sequencing error;           default = %.5f\n",
	  (double)EPSILON/1000 );
  printf( "-E rate of C->T sequencing error;      default = %.5f\n",
	  (double)EPSILON_D/1000 );
  printf( "Note that -E must be >= -e\n" );
  printf( "-b set to true if sample is from male; default = %d\n",
	  B_DEF );
  printf( "-g genome version (19 or 38); default=19\n" );
  printf( "-V make VCF output\n" );
  printf( "-I sampleID to assign for VCF file\n" );
  printf( "-B Bayes Factor ratio cutoff           default = %d\n",
	  BFC_DEFAULT );
  printf( "-G Genotype likelihood ratio cutoff    default = %d\n",
	  GLRC_DEFAULT );
  printf( "Sites with Bayes Factor and Genotype likelihood ratios\n" );
  printf( "better than the cutoff will be included in the output\n" );
  printf( "VCF file with PASS in the FILTER field.\n" );
  exit( 0 );
}

typedef struct gt_in_out {
  Hap* h;
  Pu_chr* puc;
  long unsigned int** nchoosek;
  char a0;
  char a1;
  double qual;
  double alt_freq;
  double bf;
} GIO;

void print_header( const double epsilon, const double epsilon_d );
void make_VCF_header(const char* invocation,
		     const char* sampleID );
int find_dynamic_max_cov( Pu_chr* puc, unsigned int perc_cov );
int autosomal_or_pseudo( const Hap* h );
void* thread_call_genotype( void* ptr );
int call_genotype( const Hap* h, const Pu_chr* puc,
		   long unsigned int** nchoosek,
		   char* a0p, char* a1p,
		   double* qualp, double* alt_freq,
		   double* bf );
void call_genotypes( GIO** gio_array, size_t n_gio );
void write_vcf_line( GIO* gi );
void make_vcf_gt_string( GIO* gi, char* gt_string );
/* Calculate the probability of the data observed in the input
   pileup line, with the input ref and alt alleles, given the
   input genotype, described as the number of alternative alleles
   0=homozygous ref, 1=het, 2=homozygous alternative alelle */
double pPul_givenG( const Pul* p,
		    // pointer to mpileup line with data
		    const char a0,
		    // reference allele at this position
		    const char a1,
		    // alternate allele at this position
		    const unsigned int n_ref,
		    // genotype (number of alt alleles:0,1,or 2
		    long unsigned int** nchoosek
		    // nchoosek array
		    );
int max_cov  = MAX_COV_DEF;
int down_cov = DOWN_COV_DEF;
int perc_cov = -1; // -1 = not set; anything else means user
                   // set a percential cut-off
double epsilon   = (double)EPSILON/(double)1000.0; // error rate
double epsilon_d = (double)EPSILON_D/(double)1000.0; // deam error rate
int verbose = VERBOSE;
int is_male = B_DEF;
int genome_version = 19;
int VCF = 0;
double BF_cutoff = (double)BFC_DEFAULT; // Bayes Factor cutoff for VCF filter PASS
double GLR_cutoff = (double)GLRC_DEFAULT; // Genotype likelihood ratio cutoff for VCF filter PASS

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  long unsigned int** nchoosek;

  File_Src* hs; // hap data - defined in file-io.h
  GIO** gio_array; // array of GIO input/output structs
  Hap* h;
  Pu_chr* puc; // mpileup data - defined in pileup.h
  char hap_fn[MAX_FN_LEN];
  char mp_fn[MAX_FN_LEN];
  char invocation[2048];
  char sampleID[2048] = "SampleID";
  size_t invoc_strlen;
  char a0, a1;
  int hap_set = 0;
  int mp_set  = 0;
  int make_header = 0;
  double qual, alt_freq, bf; // qual ratio, alternate allele freq, bayes factor
  size_t i = 0;
  size_t gio_inx;
  int hap_read_status = 1;
  int ich;

  while( (ich=getopt( argc, argv, "H:m:v:M:d:I:e:E:b:G:B:g:p:ahV" )) != -1 ) {
    switch(ich) {
    case 'a' :
      make_header = 1;
      break;
    case 'H' :
      strcpy( hap_fn, optarg );
      hap_set = 1;
      break;
    case 'm':
      strcpy( mp_fn, optarg );
      mp_set = 1;
      break;
    case 'M' :
      max_cov = atoi( optarg );
      break;
    case 'p' :
      perc_cov = atoi( optarg );
      break;
    case 'd' :
      down_cov = atoi( optarg );
      break;
    case 'g' :
      genome_version = atoi( optarg );
      break;
    case 'h':
      help();
      break;
    case 'e' :
      epsilon = strtod( optarg, NULL );
      break;
    case 'E' :
      epsilon_d = strtod( optarg, NULL );
      break;
    case 'b' :
      is_male = atoi( optarg );
      break;
    case 'v' :
      verbose = atoi( optarg );
      break;
    case 'V' :
      VCF = 1;
      break;
    case 'I' :
      strcpy( sampleID, optarg );
      break;
    case 'B' :
      BF_cutoff = strtod( optarg, NULL );
      break;
    case 'G' :
      GLR_cutoff = strtod( optarg, NULL );
      break;
    default :
      help();
    }
  }

  /* Reconstruct the invocation string */
  invocation[0] = '\0';
  for( i = 0; i < argc; i++ ) {
    strcat( invocation, argv[i] );
    invoc_strlen = strlen( invocation );
    invocation[invoc_strlen] = ' '; // add space after each element
    invocation[invoc_strlen+1] = '\0'; // null terminate
  }
  invoc_strlen = strlen( invocation ); // remove space after final
  invocation[invoc_strlen - 1] = '\0'; //        element
  
  /* Just make header? */
  if (make_header) {
    print_header( epsilon, epsilon_d );
  }

  /* Get input files? */
  if ( !(hap_set && mp_set) ) {
    help();
  }

  /* Read in the mpileup data, making a dictionary */
  fprintf( stderr, "Reading mpileup data in: %s\n", mp_fn );
  puc = init_Pu_chr( mp_fn );
  if ( puc == NULL ) {
    fprintf( stderr, "Problem reading mpileup data in: %s\n",
	     mp_fn );
    exit( 1 );
  }
  fprintf( stderr, "Read %lu mpileup lines\n", puc->n_puls );

  /* If user specific -p option, we need to compute the
     observed coverage distribution over the mpileup data
     in puc and set the max_cov from that */
  if ( perc_cov != -1 ) {
    max_cov = find_dynamic_max_cov( puc, perc_cov );
  }

  /* Set up nchoosek */
  nchoosek = init_nchoosek( max_cov );
  
  /* Set up the GIO array & malloc memories */
  gio_array = (GIO**)malloc(sizeof(GIO*) * MAX_HAP_BUFFER);
  for (gio_inx = 0; gio_inx < MAX_HAP_BUFFER; gio_inx++) {
    gio_array[gio_inx] = (GIO*)malloc(sizeof(GIO));
    gio_array[gio_inx]->h        = init_Hap();
    gio_array[gio_inx]->puc      = puc;
    gio_array[gio_inx]->nchoosek = nchoosek;
    gio_array[gio_inx]->qual     = 0.0;
    gio_array[gio_inx]->alt_freq = 0.0;
    gio_array[gio_inx]->bf       = 0.0;    
  }
  gio_inx = 0;
  fprintf( stderr, "Set up GIO array for buffer of %u records\n",
	   MAX_HAP_BUFFER );

  hs = init_FS( hap_fn );
  if ( hs == NULL ) {
    fprintf( stderr, "Problem reading haplotype records in: %s\n",
	     hap_fn );
  }

  if (VCF) {
    make_VCF_header(invocation, sampleID);
  }

  i = 0;
  while ( hap_read_status == 1 ) {
    if ( gio_inx == MAX_HAP_BUFFER ) {
      call_genotypes( gio_array, gio_inx );
      gio_inx = 0;
    }
    else {
      hap_read_status = read_next_Hap( gio_array[gio_inx]->h, hs );
      gio_inx++;
      i++;
    }
  }
  gio_inx--; // was incremented falsely on the last EOF read
  i--;
  call_genotypes( gio_array, gio_inx );

  fprintf( stderr, "Read %lu haplotype records\n", i );
  
  destroy_FS( hs );
  exit( 0 );
}

void make_VCF_header(const char* invocation, const char* sampleID ) {
  time_t current_time;
  current_time = time(NULL);
  printf("##fileformat=VCFv4.1\n");
  printf("##filedate=%s",ctime(&current_time));
  printf("##source=astrea-impute2 v%d\n", VERSION);
  printf("##astrea-impute2Invocation=%s\n", invocation);
  printf("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
  printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",
	 sampleID );
  return;
}

void print_header( const double epsilon, const double epsilon_d ) {
  printf( "# astrea-impute2 VERSION %d\n", VERSION );
  printf( "# (c) 2020 Astrea Forensics, Ed Green\n" );
  if ( genome_version == 38 ) {
    printf( "# Genome: hg38 / GRC38 / build 38.1\n" );
  }
  else {
    printf( "# Genome: hg19 / hs37d5 / build 37.1\n" );
  }
  printf( "# is_male: %d\n", is_male );
  printf( "# max_cov: %d\n", max_cov );
  printf( "# down_cov: %d\n", down_cov );
  printf( "# epsilon: %1.5f\n", epsilon );
  printf( "# epsilon_d: %1.5f\n", epsilon_d );
  printf( "rsID\tchromosome\tposition\tallele1\tallele2\todds-ratio\talt_freq\tbayes_factor\n" );
  exit( 0 );
}

int find_dynamic_max_cov( Pu_chr* puc, unsigned int perc_cov ) {
  size_t i;
  float cum_perc = 0.0;
  unsigned int* cov_dist;
  unsigned int cum_sites;
  
  if ( (perc_cov >=0) && (perc_cov <= 100) ) {
    // looks copacetic
    cov_dist =
      (unsigned int*)malloc(sizeof(unsigned int) * (MAX_COV+1));
			     
    /* initialize cov_dist to 0 */
      for( i = 0; i < MAX_COV; i++ ) {
	cov_dist[i] = 0;
      }
      i = 0; // reset this generic index variable
      /* Go through pileup data, incrementing counts
	 in cov_dist */
      while( i < puc->n_puls ) {
	cov_dist[ puc->puls[i]->cov ]++;
	i++;
      }
      /* Find maximum coverage level that is <= to
	 perc_cov over all sites */
      cum_sites = 0;
      for( i = 0; i < MAX_COV; i++ ) {
	cum_sites += cov_dist[i];
	cum_perc = 100.0 * (float)cum_sites/(float)puc->n_puls;
	if ( cum_perc >= (float)perc_cov ) {
	  fprintf( stderr,
		   "Setting maximum pileup coverage cutoff to %lu\n", i);
	  free( cov_dist );
	  return (int)i;
	}
      }
  }
  else {
    fprintf( stderr, "-p must be between 0 and 100\n" );
    exit( 2 );
  }
}

void* thread_call_genotype( void* ptr ) {
  GIO* gi;
  gi = (GIO*) ptr;
  call_genotype( gi->h, gi->puc, gi->nchoosek,
		 &gi->a0, &gi->a1, &gi->qual,
		 &gi->alt_freq, &gi->bf );
}

/* Args: Hap* definition of haplotype
         Pu_chr* with pileup data
         char* a0p, char* a1p, double* qualp, and double* alt_freq
	 are where we should put two allele calls, quality score,
	 and alternate allele haplotype frequencies
 */
int call_genotype( const Hap* h, const Pu_chr* puc,
		   long unsigned int** nchoosek,
		   char* a0p, char* a1p,
		   double* qualp, double* alt_freq,
		   double* bf ) {

  size_t pos_inx, target_inx, hap1_inx, hap2_inx;
  unsigned int num_a0;
  unsigned int num_ref_alleles;
  unsigned int target_num_ref_alleles;
  unsigned int pairwise_haps = 0;
  unsigned int total_pairwise_haps = 0;
  double pDgivenG[MAX_SITES][3]; // prob data at each site given
  //                                genotypes with 0,1,or2 alt alleles
  double best_p[3]; // best hap pair prob for genotypes
  //                   with 0, 1, or 2 ref alleles
  double tot_p[3]; // total prob of each target genotype, including prior
  double bayes_p[3]; // total prob of each target genotype, without prior
  double hap_pair_prior;
  best_p[0] = 0.0;
  best_p[1] = 0.0;
  best_p[2] = 0.0;
  tot_p[0] = 0.0;
  tot_p[1] = 0.0;
  tot_p[2] = 0.0;      
  bayes_p[0] = 0.0;
  bayes_p[1] = 0.0;
  bayes_p[2] = 0.0;
  
  char hap_string[MAX_SITES+1];
  size_t best_h1inx_g[3];
  size_t best_h2inx_g[3];
  double p;
  Pul* pul; // pileup line pointer

  total_pairwise_haps = find_total_pairwise_haps( h );
  
  /* For each position in the haplotype, calculate the likelihood
     of the data, given the 3 possible genotypes:
     homozygous 0, het, homozygous 1
     This saves time because the likelihood of each haplotype pair
     is simply a product of these. */
  for (pos_inx = 0; pos_inx < h->n_sites; pos_inx++) {
    if (h->pos[pos_inx] == h->target_pos) {
      target_inx = pos_inx;
    }
    /* For each possible number of ref alleles, num_a0 */
    pul = fetch_Pul( puc, h->pos[pos_inx] );
    for (num_a0 = 0; num_a0 <= 2; num_a0++ ) {
      pDgivenG[pos_inx][num_a0] =
	pPul_givenG( pul,
		     h->a0s[pos_inx],
		     h->a1s[pos_inx],
		     num_a0, nchoosek );
    }
    if (verbose) {
      printf("%d\t%c\t%c\t%s\t%1.2e\t%1.2e\t%1.2e\n",
	     h->pos[pos_inx],
	     h->a0s[pos_inx],
	     h->a1s[pos_inx],
	     &pul->bases[0],
	     pDgivenG[pos_inx][0],
	     pDgivenG[pos_inx][1],
	     pDgivenG[pos_inx][2] );
    }
  }

  /* For each pair of haplotypes, determine the implied genotype
     over each position. Use the pDgivenG over each position to
     find the aggregate probability of the data given this 
     haplotype pair. Keep track of the best (most probable) haplotype
     pair for each Target site genotype. */
  for( hap1_inx = 0; hap1_inx < h->n_haps; hap1_inx++ ) {
    for( hap2_inx = hap1_inx; hap2_inx < h->n_haps; hap2_inx++ ) {
      p = 1.0;
      for( pos_inx = 0; pos_inx < h->n_sites; pos_inx++ ) {
	num_ref_alleles = 2; // how many ref alleles does this hap pair
	                     // this haplotype pair have at this site?
	num_ref_alleles -= (h->hap_alleles[hap1_inx][pos_inx] +
			    h->hap_alleles[hap2_inx][pos_inx]);

	if ( pos_inx == target_inx ) {
	  target_num_ref_alleles = num_ref_alleles;
	}
	/* The prob of this hap pair *= prob of genotype at
	   this position given the mpilup data in pDgivenG */
	p *=  pDgivenG[pos_inx][num_ref_alleles];
      }
      /* We've summed p over all positions */
      /* Add prob of this hap pair to appropriate target genotype */
      bayes_p[target_num_ref_alleles] += p;
      
      /* Find prior for these two haplotypes */
      pairwise_haps = (h->hap_counts[hap1_inx] * h->hap_counts[hap2_inx]);
      if ( hap1_inx != hap2_inx ) { // x,y + y,x 
	pairwise_haps *= 2;
      }
      hap_pair_prior = (double)pairwise_haps / (double)total_pairwise_haps;
      /* prob of the target_num_ref_alleles genotype is 
	 prob D|hap pair x prob hap pair */
      p *= hap_pair_prior;

      tot_p[target_num_ref_alleles] += p;
      /* NEED TO SUM OVER ALL INSTEAD OF FINDING TWO BEST */
      if ( p > best_p[target_num_ref_alleles] ) {
	best_p[target_num_ref_alleles] = p;
	best_h1inx_g[target_num_ref_alleles] = hap1_inx;
	best_h2inx_g[target_num_ref_alleles] = hap2_inx;
      }
    }
  }

  if ( verbose >= 1 ) {
    for( target_num_ref_alleles = 0; target_num_ref_alleles < 3;
	 target_num_ref_alleles++ ) {
      printf( "Total prob for %u ref alleles: %4.3e\n",
	      target_num_ref_alleles, tot_p[target_num_ref_alleles] );
      printf( "Best haplotypes for %u ref alleles - probability = %4.3e\n",
	      target_num_ref_alleles, best_p[target_num_ref_alleles] );
      printf( "Hap1: %s %u\n",
	      stringify_hap( h, best_h1inx_g[target_num_ref_alleles],
			     hap_string ),
	      h->hap_counts[best_h1inx_g[target_num_ref_alleles]] );

      printf( "Hap2: %s %u\n",
	      stringify_hap( h, best_h2inx_g[target_num_ref_alleles],
			     hap_string ),
	      h->hap_counts[best_h2inx_g[target_num_ref_alleles]] );
    }
  }

  /* Set the alt_freq - this is from the haplotype table and doesn't
     depend on what our genotype call is */
  *alt_freq = target_alt_freq( h );
  
  /* Now, call genotype and its score from data in tot_p[] */
  if ( tot_p[0] >= tot_p[1] ) { // 0 > 1
    if ( tot_p[1] >= tot_p[2] ) { // 0 > 1 > 2
      /* homozygous alt */
      *a0p = h->Ta1;
      *a1p = h->Ta1;
      *qualp = tot_p[0]/tot_p[1];
      *bf    = bayes_p[0]/bayes_p[1];
      return 1;
    }
    else { // 0,2 > 1
      if ( tot_p[0] >= tot_p[2] ) { // 0 > 2 > 1
	*a0p = h->Ta1;
	*a1p = h->Ta1;
	*qualp = tot_p[0]/tot_p[2];
	*bf    = bayes_p[0]/bayes_p[2];
	return 1;
      }
      else { // 2 > 0 > 1
	*a0p = h->Ta0;
	*a1p = h->Ta0;
	*qualp = tot_p[2]/tot_p[0];
	*bf    = bayes_p[2]/bayes_p[0];	
	return 1;
      }
    }
  }
  else { // 1 > 0
    if ( tot_p[0] >= tot_p[2] ) { // 1 > 0 > 2
      /* Can't be het if this position is in the 
	 non-pseudoautsomal region of a sex chromosome
	 and the sample is a male */
      if ( is_male &&
	   !autosomal_or_pseudo(h) ) {
	*a0p = h->Ta0;
	*a1p = h->Ta0;
	*qualp = tot_p[0]/tot_p[2];
	*bf    = bayes_p[0]/bayes_p[2];
	return 1;
      }
      else {
	*a0p = h->Ta0;
	*a1p = h->Ta1;
	*qualp = tot_p[1]/tot_p[0];
	*bf    = bayes_p[1]/bayes_p[0];
	return 1;
      }
    }
    else { // 1,2 > 0
      if ( tot_p[1] >= tot_p[2] ) { // 1 > 2 > 0
      /* Check for sex chromosome coordinates */
	if ( is_male && !autosomal_or_pseudo(h) ) {
	  *a0p = h->Ta1;
	  *a1p = h->Ta1;
	  *qualp = tot_p[2]/tot_p[0];
	  *bf    = bayes_p[2]/bayes_p[0];
	  return 1;
	}
	else {
	  *a0p = h->Ta0;
	  *a1p = h->Ta1;
	  *qualp = tot_p[1]/tot_p[2];
	  *bf    = bayes_p[1]/bayes_p[2];
	  return 1;
	}
      }
      else { // 2 > 1 > 0
	*a0p = h->Ta0;
	*a1p = h->Ta0;
	*qualp = tot_p[2]/tot_p[1];
	*bf    = bayes_p[2]/bayes_p[1];
	return 1;
      }
    }
  }
  fprintf( stderr, "Cannot decide genotype %.5f %.5f %.5f",
	   tot_p[0], tot_p[1], tot_p[2] );
  return 0;
}

void call_genotypes( GIO** gio_array, size_t n_gio ) {
  GIO* gi;
  pthread_t thread_id[MAX_HAP_BUFFER];
  size_t i;
  for( i = 0; i < n_gio; i++ ) {
    pthread_create( &thread_id[i], NULL,
		    thread_call_genotype, gio_array[i] );
  }
  for( i = 0; i < n_gio; i++ ) {
    pthread_join( thread_id[i], NULL );
  }
		  
  /* Write output */
  for( i = 0; i < n_gio; i++ ) {
    gi = gio_array[i];
    if (VCF) {
      /* Implement */
      write_vcf_line( gi );
    }
    else {
      printf( "%s\t%s\t%u\t%c\t%c\t%4.3e\t%4.3e\t%4.3e\n",
	      gi->h->TrsID, gi->h->chr, gi->h->target_pos,
	      gi->a0, gi->a1, gi->qual, gi->alt_freq, gi->bf );
    }
  }
}

void write_vcf_line( GIO* gi ) {
  char filter_status[64];
  char gt_string[64];

  /* Apply VCF filters */
  if ( gi->bf < BF_cutoff ) {
    strcpy( filter_status, "BF" );
    return;
  }
  else {
    if ( gi->qual < GLR_cutoff ) {
      strcpy( filter_status, "GLR" );
      return;
    }
    else {
      strcpy( filter_status, "PASS" );
    }
  }

  make_vcf_gt_string( gi, gt_string );
  
  printf("%s\t", gi->h->chr);
  printf("%u\t", gi->h->target_pos);
  printf("%s\t", gi->h->TrsID);
  printf("%c\t", gi->h->Ta0);
  printf("%c\t", gi->h->Ta1);
  printf(".\t" );
  printf("%s\t", filter_status);
  printf("AF=%.3f\t", target_alt_freq(gi->h) );
  printf( "GT\t" );
  printf( "%s\n", gt_string );
}

void make_vcf_gt_string( GIO* gi, char* gt_string ) {
  gt_string[1] = '/';
  gt_string[3] = '\0';
  if (gi->a0 == gi->h->Ta0) {
    gt_string[0] = '0';
  }
  else {
    if (gi->a0 == gi->h->Ta1) {
      gt_string[0] = '1';
    }
    else {
      gt_string[0] = '.';
    }
  }
  if (gi->a1 == gi->h->Ta0) {
    gt_string[2] = '0';
  }
  else {
    if (gi->a1 == gi->h->Ta1) {
      gt_string[2] = '1';
    }
    else {
      gt_string[2] = '.';
    }
  }
  return;
}

/* Return TRUE if pos is in an autosome or the pseudoautosomal 
   region of a sex chromosome. We only call this if the sample is a
   male, so no need to check that */
int autosomal_or_pseudo( const Hap* h ) {
  if ((h->chr[0] == 'X') || (h->chr[4] == 'X')) {
    /* It's X. Is it in pseudoautosomal regions? */
    if ( genome_version == 19 ) {
      if ( ((h->target_pos >= XPAR1ST) && (h->target_pos <= XPAR1EN))
	   ||
	   ((h->target_pos >= XPAR2ST) && (h->target_pos <= XPAR2EN)) ) {
	return 1;
      }
      else {
	return 0;
      }
    }
    if ( genome_version == 38 ) {
      if ( ((h->target_pos >= XPAR1ST_HG38) && (h->target_pos <= XPAR1EN_HG38))
	   ||
	   ((h->target_pos >= XPAR2ST_HG38) && (h->target_pos <= XPAR2EN_HG38)) ) {
	return 1;
      }
      else {
	return 0;
      }
    }
  }
  if  ((h->chr[0] == 'Y') || (h->chr[4] == 'Y')) {
    /* It's Y. Is it in pseudoautosomal regions? */
    if ( genome_version == 19 ) {
      if ( ((h->target_pos >= YPAR1ST) && (h->target_pos <= YPAR1EN))
	   ||
	   ((h->target_pos >= YPAR2ST) && (h->target_pos <= YPAR2EN)) ) {
	return 1;
      }
      else {
	return 0;
      }
    }
    if ( genome_version == 38 ) {
      if ( ((h->target_pos >= YPAR1ST_HG38) && (h->target_pos <= YPAR1EN_HG38))
	   ||
	   ((h->target_pos >= YPAR2ST_HG38) && (h->target_pos <= YPAR2EN_HG38)) ) {
	return 1;
      }
      else {
	return 0;
      }
    }
  }
  return 1;   /** Neither X nor Y **/
}


/* Calculate the probability of the data observed in the input
   pileup line, with the input ref and alt alleles, given the
   input genotype, described as the number of alternative alleles
   0=homozygous ref, 1=het, 2=homozygous alternative alelle */
double pPul_givenG( const Pul* p,
		    // pointer to mpileup line with data
		    const char a0,
		    // reference allele at this position
		    const char a1,
		    // alternate allele at this position
		    const unsigned int n_ref,
		    // genotype (number of alt alleles:0,1,or 2
		    long unsigned int** nchoosek
		    // nchoosek array
		    ) {
  unsigned int obs_ref = 0;
  unsigned int obs_alt = 0;
  double p_ref         = 1.0; // prob of a ref allele, given genotype
  // epsilon & epsilon_d are now globals
  //double epsilon       = (double)EPSILON/(double)1000.0; // error rate
  //double epsilon_d     = (double)EPSILON_D/(double)1000.0;
                         // deamination error rate
  unsigned int mq_cut  = MQ_CUT;
  unsigned int bq_cut  = BQ_CUT;

  double deam_risk_obs    = 0;
  double no_deam_risk_obs = 0;

  size_t i;

  int CT_site = 0;
  int GA_site = 0;
  int valid_allele = 0;

  if ( p == NULL ) {
    return 1.0; // no sequence data, this is fully probable
    // regardless of model
  }
  if ( p->cov > max_cov ) {
    return 1.0; // too much sequence data to classify, these
    // data are fully probable under any model
  }
  
  if ( ((a0 == 'C') && (a1 == 'T')) ||
       ((a0 == 'T') && (a1 == 'C')) ) {
    CT_site = 1;
  }
  if ( ((a0 == 'G') && (a1 == 'A')) ||
       ((a0 == 'A') && (a1 == 'G')) ) {
    GA_site = 1;
  }

  /* Go through each base in the Pul and classify it */
  for( i = 0; (i < p->cov) && (i < down_cov); i++ ) {
    valid_allele = 0;
    if ( (p->base_quals[i] > bq_cut) &&  /* Check filters */
	 (p->map_quals[i]  > mq_cut) &&
	 (p->cov <= max_cov) ) {
      /* Observation that could be affected by deamination?
	 These are positive strand alignments at C/T sites
	 and negative strand alignments at G/A sites */
      if ( p->bases[i] == a0 ) {
	valid_allele = 1;
	obs_ref++;
      }
      else {
	if ( p->bases[i] == a1 ) {
	  valid_allele = 1;
	  obs_alt++;
	}
      }
      if ( valid_allele ) {
	if ( (CT_site && (p->strands[i] == 1)) ||
	     (GA_site && (p->strands[i] == -1)) ) {
	  deam_risk_obs++;
	}
	else {
	  no_deam_risk_obs++;
	}
      }
    }
  }
    
  /* Calculate the overall probability of a reference base, given
     the strand observations we've seen and the genotype */
  if ( n_ref == 2 ) { // We're calculating probability of homozygous reference
    if ( (CT_site && (a0 == 'C')) ||
	 (GA_site && (a0 == 'G')) ) {
      /* The deamination risk is making fewer reference observations */
      p_ref =
	( ((1.0 - epsilon) * no_deam_risk_obs) +
	  ((1.0 - (epsilon_d - epsilon)) * deam_risk_obs) )
	/
	( no_deam_risk_obs + deam_risk_obs );
    }
    else { /* Deamination is making more reference observations.
	      But we're calculating what we should see under the model
	      of homozygous reference, so ignore deamination */
      p_ref =
	(1.0 - epsilon);
    }
  }

  if (n_ref == 1 ) { // We're calculating probability of het
    if ( CT_site || GA_site ) {
      if ( (CT_site && (a0 == 'C')) ||
	   (GA_site && (a0 == 'G')) ) {
	/* The deamination risk is making fewer reference observations */
	p_ref = 
	  ( (0.5 * no_deam_risk_obs) +
	    ((0.5 - (epsilon_d - epsilon)) * deam_risk_obs) )
	  /
	  ( no_deam_risk_obs + deam_risk_obs );
      }
      else { 
	/* The deamination risk is making more reference observations */
	p_ref =
	  ( (0.5 * no_deam_risk_obs) +
	    ((0.5 + (epsilon_d - epsilon)) * deam_risk_obs) )
	  /
	  ( no_deam_risk_obs + deam_risk_obs );
      }
    }
    else { /* Transversion. Allele balance is only affected by loss of
	      C on plus or G on minus alignments */
      p_ref = 0.5; // THIS NEEDS FURTHER LOGIC
    }
  }
  if (n_ref == 0 ) { // We're calculating probability of homozygous alternate
    if ( (CT_site && (a1 == 'C')) ||
	 (GA_site && (a1 == 'G')) ) {
      /* The deamination risk is making more reference observations */
      p_ref =
	( ((epsilon) * no_deam_risk_obs) +
	  (((epsilon_d - epsilon)) * deam_risk_obs) )
	/
	(no_deam_risk_obs + deam_risk_obs);
    }
    else {
      /* The deamination risk is making fewer reference observations */
      p_ref = epsilon;
    }

  }

  if ( verbose >= 2 ) {
    fprintf( stderr, "P(num ref alleles = %u) obs_ref %u obs_alt %u %0.10f\n",
	     n_ref, obs_ref, obs_alt,
	     binomial_p( obs_ref + obs_alt, obs_ref, p_ref, nchoosek ) );
  }
  if ( (obs_ref + obs_alt) == 0 ) {
    return 1.0;
  }
  if ( isfinite( p_ref ) &&
       isnormal( p_ref ) ) {
    return binomial_p( obs_ref + obs_alt, obs_ref, p_ref, nchoosek );
  }
  else {
    fprintf( stderr, "Won't call binomial_p with p = %4.4e\n", p_ref );
    return 1.0;
  }
}


