#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "pileup.h"

void help( void ) {
  printf( "test-pileup -m <mpileup file>\n" );
  exit( 0 );
}

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  int ich;
  int input_set = 0;
  char pu_fn[MAX_FN_LEN];
  char pu_str[MAX_FIELD_WIDTH*2];
  Pu_chr* puc;
  size_t i = 0;
  Pul* pl;
  
  while( (ich=getopt( argc, argv, "m:" )) != -1 ) {
    switch(ich) {
    case 'm' :
      strcpy( pu_fn, optarg );
      input_set = 1;
      break;
    default :
      help();
    }
  }
  if ( !input_set ) {
    help();
  }
  
  puc = init_Pu_chr( pu_fn );
  
  /* Say something about puc */
  printf( "Parsed %lu mpileup lines\n", puc->n_puls );

  /* */
  pl = fetch_Pul( puc, 9411409 );

  exit( 0 );
}
