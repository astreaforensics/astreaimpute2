#include <stdlib.h>                                                          
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "local-impute.h"

void help( void ) {
  printf( "test-local-impute1\n" );
  printf( "Runs some tests of the local-impute code\n" );
  exit( 0 );
}

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  extern unsigned long long* fact;
  Impute_Set* impute_set;
  int ich;
  size_t i;
  
  while( (ich=getopt( argc, argv, "k:f:h" )) != -1 ) {
    switch(ich) {
    case 'h' :
      help();
    default :
      help();
    }
  }

  /* Initialization */
  impute_set = init_Impute_Set();
  fact = init_factorial();
  
  
  /* Add some Target Site data */
  if( add_Target_to_Impute_Set( impute_set,
				"21",
				100,
				"rs123",
				90,
				110 ) ) {
    printf( "Problem adding Target site data\n" );
  }

  /* Analyze Target Site */
  printf( "Probability of genotypes at Target Site" );
  printf( "Given count(0) = %d and count(1) = %d\n",
	  impute_set->T0, impute_set->T1 );
  printf( "p(0,0) = %.4f\n", p_T00(impute_set) );
  printf( "p(0,1) = %.4f\n", p_T01(impute_set) );
  printf( "p(1,1) = %.4f\n", p_T11(impute_set) );


  /* Add some data from Linked Sites sites */
  if ( add_impute_pair( impute_set,
			"21",
			200,
			"rs124",
			99, // T0L0
			10,  // T0L1
			9, // T1L0
			99, // T1L1
			0 ) ) {
    printf( "Problem adding Linked site data\n" );
  }
  printf( "D of linked site1: %.4f\n",
	  calc_D(&impute_set->impute_pairs[0]) );
  
  if ( add_impute_pair( impute_set,
			"21",
			300,
			"rs125",
			98, // T0L0
			21,  // T0L1
			23, // T1L0
			97, // T1L1
			0 ) ) {
    printf( "Problem adding Linked site data\n" );
  }
  printf( "D of linked site2: %.4f\n",
	  calc_D(&impute_set->impute_pairs[1]) );
  
  if ( add_impute_pair( impute_set,
			"21",
			400,
			"rs125",
			99, // T0L0
			10,  // T0L1
			19, // T1L0
			99, // T1L1
			0 ) ) {
    printf( "Problem adding Linked site data\n" );
  }
  printf( "D of linked site1: %.4f\n",
	  calc_D(&impute_set->impute_pairs[0]) );

  /* Add some allele observations at linked sites */
  add_linked_counts( &impute_set->impute_pairs[0], 4, 5 );
  add_linked_counts( &impute_set->impute_pairs[1], 3, 2 );
  add_linked_counts( &impute_set->impute_pairs[2], 1, 1 );

  /* Calculate the p_L given each Target site genotype */
  for(i = 0; i < impute_set->num_pairs; i++ ) {
    printf( "For Linked site %lu\n", i );
    printf( "  %.8f = p(L) | T=0,0\n",
	    p_LgivenT00( &impute_set->impute_pairs[i] ) );
    printf( "  %.8f = p(L) | T=0,1\n",
	    p_LgivenT01( &impute_set->impute_pairs[i] ) );
    printf( "  %.8f = p(L) | T=1,1\n",
	    p_LgivenT11(&impute_set->impute_pairs[i] ) );
  }

  printf( "P(T=0,0), P(T=0,1), P(T=1,1) | all linked sites data:\n" );
  printf( "%.15f, %.15f, %.15f\n",
	  p_allLgivenT00( impute_set ),
	  p_allLgivenT01( impute_set ),
	  p_allLgivenT11( impute_set ) );
	  
  exit( 0 );

}
