#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "hap.h"
void help( void ) {
  printf( "hap2af-table -h <astrea-impute2 haplotype file>\n" );
  printf( "Parses a hap table and generates output for each\n" );
  printf( "Target sites in the following tab-delimited format:\n" );
  printf( "CHR POS ID REF ALT AF\n" );
  printf( "where AF is the alternate allele frequency.\n" );
  exit( 0 );
}
int main ( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  int ich;
  int input_set = 0;
  char hap_fn[MAX_FN_LEN];
  Hap* h;
  File_Src* hs;
  int read_status  = 0;
  while( (ich=getopt( argc, argv, "h:" )) != -1 ) {
    switch(ich) {
    case 'h' :
      strcpy( hap_fn, optarg );
      input_set = 1;
      break;
    default :
      help();
    }
  }
  if ( !input_set ) {
    help();
  }
  h = init_Hap();
  hs = init_FS( hap_fn );
  read_status = read_next_Hap( h, hs );
  while( read_status == 1) {
    printf( "%s\t%u\t%s\t%c\t%c\t%.3f\n",
            h->chr,
            h->target_pos,
            h->TrsID,
            h->Ta0,
            h->Ta1,
            target_alt_freq( h ) );
    read_status = read_next_Hap( h, hs );
  }
  destroy_FS( hs );
  exit( 0 );
}
