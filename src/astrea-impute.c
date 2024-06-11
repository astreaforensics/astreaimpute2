#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "local-impute.h"
#include "file-io.h"

#define VERSION (18)
#define TABLE_TYPE (4)
#define LOW_DEF (0)
#define HIGH_DEF (256)

void help( double min_D, double min_ratio ) {
  printf( "astrea-impute v%d\n", VERSION );
  printf( "(c) Astrea Forensics 2019, 2020\n" );
  printf( "OPTIONS:\n" );
  printf( "-f Input file of linked sites and allelic observations; see below\n" );
  printf( "-t table type; 1(default)=just genotype likelihods\n" );
  printf( "               2         =full likelihoods from each site>\n" );
  printf( "               3         =Genotype calls\n" );
  printf( "-m minimum D for linked site; default = %.3f\n", min_D );
  printf( "   Any loci whose D (linkage disequilibrium) compared to the\n" );
  printf( "   target site is less than this information will not be used\n" );
  printf( "   when determining the Target site genotype.\n" );
  printf( "-r minimum ratio of best call to second best call to\n" );
  printf( "   make a full genotype call; default = %.5f\n", min_ratio );
  printf( "-R minimum ratio, as above, at zero coverage sites\n" );
  printf( "   to make a full genotype call; default = %.5f\n", min_ratio );
  printf( "   This cutoff is used for sites where the Target site has\n" );
  printf( "   no sequence data.\n" );
  printf( "-u Minimum genotype likelihood ratio cutoff for only using\n" );
  printf( "   Target Site data. If the genotype is clear already, using\n" );
  printf( "   only data from the Target Site, do not consider data from\n" );
  printf( "   Linked Sites. This is useful for high coverage data.\n" );
  printf( "   default = %.5f\n", min_ratio );
  printf( "   This ratio is used to determine whether a genotype can be\n" );
  printf( "   unambiguously called from the Target site data alone without\n" );
  printf( "   considering information from linked sites.\n" );
  printf( "-v ratio of pHet vs pHom for unobserved allele to call as\n" );
  printf( "   homozygous for the observed allele. This can be used to\n" );
  printf( "   increase the called homozygosity in a principled way.\n" );
  printf( "   If not set, this rule does is not used.\n" );
  printf( "-L low coverage cutoff to include data at a linked site. If the\n" );
  printf( "   total observed coverage at a linked site is less than this\n" );
  printf( "   amount, then this linked site data will be skipped.\n" );
  printf( "   default = %d\n", LOW_DEF );
  printf( "-H high coverage cutoff to include data at a linked site. If the\n" );
  printf( "   total observed coverage at a linked site is higher than this\n" );
  printf( "   amount, then this linked site data will be skipped.\n" );
  printf( "   default = %d\n", HIGH_DEF );
  printf( "-T minimum Target Site coverage to make output, regardless\n" );
  printf( "   of the ratio values. If the Target Site coverage is less\n" );
  printf( "   than this value NO OUTPUT will be made for this site.\n" );
  printf( "   default = 0\n" );
  printf( "-a maximum Target Site coverage to make output, like -T above\n" );
  printf( "   default = HIGH_DEF\n" ); 
  printf( "-Y input is male data. Heterozygosity is impossible on sex\n" );
  printf( "   chromosomes except in psuedo-autosomal regions.\n" );
  printf( "-B filename of IDs of blacklisted Targets to exclude from output\n" );
  printf( "-p if set, use the reference panel allele frequencies as priors\n" );
  printf( "   when calculating genotype probabilities\n" );
  printf( "Makes a genotype file of the Target Markers listed in the input\n" );
  printf( "file using the patterns of linkage and observed alleles at target\n" );
  printf( "and linked sites in the input file.\n" );
  exit( 0 );
}

void output_impute_set( Impute_Set* isp, const int Male,
			const int table_type,
			const double min_ratio_zero_cov,
			const double min_ratio_target_only,
			const double special_hom,
			const int use_priors );
double min_ratio = 10.0;

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  extern unsigned long long* fact;

  double min_D = 0.0001;
  double min_ratio_zero_cov    = min_ratio;
  double min_ratio_target_only = min_ratio;
  double special_hom = 0.0;
  int table_type = TABLE_TYPE;
  int min_Targ_cov = 0;
  int max_Targ_cov = HIGH_DEF;
  int Male = 0;
  int blacklist_set = 0;
  int use_priors;
  File_Src* fs; // for input haplotype table
  char input_table_fn[MAX_FN_LEN+1];
  char table_line[MAX_LINE_LEN+1];
  char bl_fn[MAX_FN_LEN+1];
  Blacklist* blacklist;
  Impute_Set* impute_set;
  int ich;
  char* tchr;
  char* lchr;
  unsigned int tpos, lpos;
  char* tid;
  char* lid;
  char la0; // linked site reference allele
  char la1; // linked site alternate allele
  short unsigned int t0l0, t0l1, t1l0, t1l1, L0, L1;
  size_t ipn;
  int input_set = 0;
  int low_cov = LOW_DEF;
  int high_cov = HIGH_DEF;
  while( (ich=getopt( argc, argv, "f:m:t:r:R:u:v:T:L:H:B:a:Yp" )) != -1 ) {
    switch(ich) {
    case 'f' :
      strcpy( input_table_fn, optarg );
      input_set = 1;
      break;
    case 'm' :
      min_D = atof( optarg );
      break;
    case 't' :
      table_type = atoi( optarg );
      break;
    case 'r' :
      min_ratio = atof( optarg );
      break;
    case 'R' :
      min_ratio_zero_cov = atof( optarg );
      break;
    case 'u' :
      min_ratio_target_only = atof( optarg );
      break;
    case 'v' :
      special_hom = atof( optarg );
      break;
    case 'T' :
      min_Targ_cov = atoi( optarg );
      break;
    case 'a' :
      max_Targ_cov = atoi( optarg );
      break;
    case 'h' :
      help(min_D, min_ratio);
    case 'L' :
      low_cov = atoi( optarg );
      break;
    case 'H' :
      high_cov = atoi( optarg );
      break;
    case 'Y' :
      Male = 1;
      break;
    case 'B' :
      strcpy( bl_fn, optarg );
      blacklist_set = 1;
      blacklist = read_blacklist( bl_fn );
      break;
    case 'p' :
      use_priors = 1;
      break;
    default :
      help(min_D, min_ratio);
    }
  }
  if ( !input_set ) {
    help(min_D, min_ratio);
  }

  /* Initialization */
  tchr =(char*)malloc(sizeof(char)*(MAX_ID_LEN+1));
  lchr =(char*)malloc(sizeof(char)*(MAX_ID_LEN+1));
  tid  =(char*)malloc(sizeof(char)*(MAX_ID_LEN+1));
  lid  =(char*)malloc(sizeof(char)*(MAX_ID_LEN+1));
  
  impute_set = init_Impute_Set();
  fact = init_factorial();

  /* Set up the input file */
  fs = init_FS( input_table_fn );
  if ( fs == NULL ) {
    fprintf( stderr, "EXIT: Problem reading %s\n", input_table_fn );
    exit( 35 );
  }
  
  /* Parse input file */
  while( get_line_FS( fs, table_line ) != NULL ) {
    if ( table_line[0] != '#' ) { // skip comment lines that begin with #
      sscanf( table_line,
	      "%s %u %s %s %u %s %s %s %hu %hu %hu %hu %hu %hu\n",
	      tchr, &tpos, tid, lchr, &lpos, lid, &la0, &la1,
	      &t0l0, &t0l1, &t1l0, &t1l1, &L0, &L1 );
      /* Special line indicating a Target site? */
      if ( strcmp( tid, lid ) == 0 ) {
	/* This is a new site */
	/* First, process the old site */
	if ( ((impute_set->impute_pairs[0].L0 + impute_set->impute_pairs[0].L1)
	      >= min_Targ_cov)
	     &&
	     ((impute_set->impute_pairs[0].L0 + impute_set->impute_pairs[0].L1)
	      <= max_Targ_cov) ) {
	  if ( !blacklist_set ||
	       black_check( blacklist, impute_set->T_rsID ) ) {
	    output_impute_set( impute_set, Male, table_type,
			       min_ratio_zero_cov,
			       min_ratio_target_only,
			       special_hom,
			       use_priors );
	  }
	}
	  /* Now, clear the impute_set and add data for this new guy */
	clear_Impute_Set( impute_set );
	if ( add_Target_to_Impute_Set( impute_set, tchr, tpos, tid,
				       (t0l0 + t0l1),
				       (t1l0 + t1l1) ) ) {
	  fprintf( stderr, "Problem parsing:\n%s %s\n", tid, lid );
	}
	if ( add_impute_pair( impute_set, lchr, lpos, lid, la0, la1,
			      t0l0, t0l1, t1l0, t1l1, min_D, 1 ) ) {
	  fprintf( stderr, "Problem adding Target: %s\n", lid );
	}
	else {
	  add_linked_counts( &impute_set->impute_pairs[0], L0, L1 );
	}
      }
      /* Linked site, try to add it as impute_pair, if there is any data */
      else {
	if ( ((L0 + L1) >= low_cov) && ((L0 + L1) <= high_cov) ) {
	  if ( !add_impute_pair( impute_set, lchr, lpos, lid, la0, la1,
				 t0l0, t0l1, t1l0, t1l1, min_D, 0 ) ) {
	    ipn = impute_set->num_pairs - 1;
	    add_linked_counts( &impute_set->impute_pairs[ipn], L0, L1 );
	  }
	}
      }
    }
  }

  /* Process the last one. It wasn't triggered by the "next" one since
     this is EOF */
  if ( ((impute_set->impute_pairs[0].L0 + impute_set->impute_pairs[0].L1)
	>= min_Targ_cov)
       &&
       ((impute_set->impute_pairs[0].L0 + impute_set->impute_pairs[0].L1)
	<= max_Targ_cov) ) {
    if ( !blacklist_set ||
	 black_check( blacklist, impute_set->T_rsID ) ) {
      output_impute_set( impute_set, Male, table_type,
			 min_ratio_zero_cov,
			 min_ratio_target_only,
			 special_hom,
			 use_priors );
    }
  }
  
  destroy_FS( fs );
  exit( 0 );
}

/* Takes a pointer to Impute_Set
   Writes output to STDOUT about this site.
   Output format depends on global table_type */
void output_impute_set( Impute_Set* isp, const int Male,
			const int table_type,
			const double min_ratio_zero_cov,
			const double min_ratio_target_only,
			const double special_hom,
			const int use_priors ) {
  Impute_Pair* imp;
  char a0 = '0';
  char a1 = '0';
  size_t i;
  if ( table_type == 1 ) {
    if ( isp->num_pairs > 0 ) {
      printf( "%s %u %s %e %e %e %.5f %.5f %.5f\n",
	      isp->chr, isp->T_pos, isp->T_rsID,
	      p_allLgivenT00( isp ),
	      p_allLgivenT01( isp ),
	      p_allLgivenT11( isp ),
	      p_T00( isp ),
	      p_T01( isp ),
	      p_T11( isp ) );
    }
  }
  if ( table_type == 2 ) {
    for( i = 0; i < isp->num_pairs; i++ ) {
      imp = &isp->impute_pairs[i];
      printf( "%s %u %s %s %u %s %c %c %u %u %u %u %u %u %.12f %.12f %.12f\n",
	      isp->chr, isp->T_pos, isp->T_rsID,
	      imp->chr, imp->L_pos, imp->L_rsID, imp->la0, imp->la1,
	      imp->T0L0, imp->T0L1, imp->T1L0, imp->T1L1,
	      imp->L0, imp->L1,
	      p_LgivenT00( imp, isp->epsilon ),
	      p_LgivenT01( imp, isp->epsilon ),
	      p_LgivenT11( imp, isp->epsilon ) );
    }
  }

  if ( table_type == 3 ) {
    if ( isp->num_pairs > 0 ) {
      isp2genotype( isp, Male, min_ratio, min_ratio_zero_cov,
		    min_ratio_target_only, special_hom,
		    &a0, &a1, use_priors );
      if ( a0 != '0' ) {
	printf( "%s\t%s\t%u\t%c\t%c\n",
		isp->T_rsID,
		isp->chr,
		isp->T_pos,
		a0, a1 );
      }
    }
  }
  return;
}
