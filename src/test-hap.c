#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "hap.h"

void help( void ) {
  printf( "test-hap -h <astrea-impute2 haplotype file -S -H>\n" );
  printf( "Tests parsability of an ai2 haplotype file.\n" );
  printf( "-S <report number of sites for each entry>\n" );
  printf( "-H <report number of haplotypes for each entry>\n" );
  exit( 0 );
}

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  int ich;
  int input_set = 0;
  int report_num_sites = 0;
  int report_num_haps  = 0;
  char hap_fn[MAX_FN_LEN];
  char max_sites_rsID[MAX_ID_LEN];
  char max_haps_rsID[MAX_ID_LEN];
  Hap* h;
  File_Src* hs;
  size_t i         = 0;
  size_t max_sites = 0;
  size_t max_haps  = 0;
  int read_status  = 0;
  unsigned int last_pos = 0;
  
  while( (ich=getopt( argc, argv, "h:SH" )) != -1 ) {
    switch(ich) {
    case 'h' :
      strcpy( hap_fn, optarg );
      input_set = 1;
      break;
    case 'S' :
      report_num_sites = 1;
      break;
    case 'H' :
      report_num_haps = 1;
      break;
    default :
      help();
    }
  }
  if ( !input_set ) {
    help();
  }

  h  = init_Hap();
  hs = init_FS( hap_fn );

  read_status = read_next_Hap( h, hs );
  while( read_status == 1) {
    i++;
    if ( h->target_pos <= last_pos ) {
      fprintf( stderr, "Target sites not sorted: %u %u\n",
	       h->target_pos, last_pos );
      exit(1);
    }
    else {
      last_pos = h->target_pos;
    }

    if ( report_num_sites ) {
      printf( "%lu\n", h->n_sites );
    }
    if ( report_num_haps ) {
      printf( "%lu\n", h->n_haps );
    }

    if ( h->n_sites > max_sites ) {
      max_sites = h->n_sites;
      strcpy( max_sites_rsID, h->TrsID );
    }
    if ( h->n_haps > max_haps ) {
      max_haps = h->n_haps;
      strcpy( max_haps_rsID, h->TrsID );
    }
    read_status = read_next_Hap( h, hs );
  }

  if ( read_status != 0 ) {
    fprintf( stderr, "%s encountered haplotype parse error\n", hap_fn );
    exit( 1 );
  }

  if ( !report_num_sites && !report_num_haps ) {
    printf( "Parsed %lu haplotype records\n", i );
    if ( h->n_haps == 0 ) {
      exit( 0 );
    }
    
    printf( "Target %s had the most sites in haplotype: %lu\n",
	    max_sites_rsID, max_sites );
    printf( "Target %s had the most haps in haplotype: %lu\n",
	    max_haps_rsID, max_haps );  
  }
    exit( 0 );
}
