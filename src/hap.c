#include "hap.h"

Hap* init_Hap( void ) {
  Hap* h;
  h = (Hap*)malloc(sizeof( Hap ));
  h->n_sites = 0;
  h->n_haps  = 0;
  return h;
}

/* Reads next hap entry from File_Src* f, putting
   result into Hap* h
   Returns 1 => copacetic, read a record
           0 => EOF
          -1 => problem
 */
int read_next_Hap( Hap* h, File_Src* f) {
  char l[MAX_LINE_LEN];
  char hap_string[MAX_SITES];
  int hap_section = 0; // boolean, set true when in hap section
                       // of record
  int end_of_record = 0; // boolean, set true when end-of-record
                         // token is observed (## E_R)
  size_t site_number = 0;
  size_t hap_number  = 0;
  /* Check first line to make sure it looks like a hap record */
  if (get_line_FS( f, l ) == NULL) {
    return 0;
  }
  if ( strncmp( l, "## TARGET", 9 ) != 0 ) {
    fprintf( stderr,
	     "%s: does not look like beginning of hap record\n",
	     l );
    return -1;
  }
  
  /* Get Target Site info */
  if ( sscanf( l, "## TARGET %d %c %c %s",
	       &h->target_pos,
	       &h->Ta0,
	       &h->Ta1,
	       h->TrsID ) == 4 ) {
    while (!hap_section) {
      get_line_FS( f, l );
      if ( strncmp( l, "## Haplotypes", 13 ) == 0 ) {
	hap_section = 1;
	h->n_sites = site_number;
      }
      else {
	sscanf( l, "%s\t%d\t%c\t%c\t%s",
		h->chr,
		&h->pos[site_number],
		&h->a0s[site_number],
		&h->a1s[site_number],
		h->rsIDs[site_number] );
	site_number++;
      }
    }
    get_line_FS( f, l ); // skip ##_____*_____ line
    while (!end_of_record) {
      get_line_FS( f, l );
      if ( strncmp( l, "## E_R", 6 ) == 0 ) {
	end_of_record = 1;
	h->n_haps     = hap_number;
	return 1;
      }
      else {
	if ( hap_number >= MAX_HAPS ) {
	  fprintf( stderr,
		   "%s has more than %d haplotypes\n",
		   h->TrsID, MAX_HAPS );
	}
	sscanf( l, "%u\t%s",
		&h->hap_counts[hap_number],
		hap_string );
	if ( string2ui( hap_string, 
			h->hap_alleles[hap_number],
			h->n_sites ) ) {
	  fprintf( stderr,
		   "Problem parsing %lu haplotype string on Target: %s\n",
		   hap_number,
		   h->TrsID );
	}
	hap_number++;
      }
    }
  }
  else {
    fprintf( stderr,
	     "%s: could not parse all parts\n", l );
    return -1;
  }
}

int string2ui( char* string, unsigned int* ha, size_t n ) {
  size_t i;
  size_t promblem = 0;
  for( i = 0; i < n; i++ ) {
    switch( string[i] ) {
    case '0' :
      ha[i] = 0;
      break;
    case '1' :
      ha[i] = 1;
      break;
    default :
      ha[i] = 2;
      fprintf( stderr, "Problem reading a haplotype string:\n\"%s\"\n",
	       string );
      promblem = 1;
    }
  }
  return promblem;
}

char* stringify_hap( const Hap* h, const size_t hap_inx, char* s ) {
  size_t i;
  for( i = 0; i < h->n_sites; i++ ) {
    if (h->hap_alleles[hap_inx][i] == 0) {
      s[i] = '0';
    }
    if (h->hap_alleles[hap_inx][i] == 1) {
      s[i] = '1';
    }
  }
  s[i] = '\0';
  return s;
}
    
unsigned int find_total_pairwise_haps( const Hap* h ) {
  size_t h1;
  unsigned int t = 0;
  for( h1 = 0; h1 < h->n_haps; h1++ ) {
    t += h->hap_counts[h1];
  }
  return t*t;
}

size_t get_target_inx( const Hap* h ) {
  size_t i;
  for( i = 0; i < h->n_sites; i++ ) {
    if ( strcmp( h->rsIDs[i], h->TrsID ) == 0 ) {
      return i;
    }
  }
  fprintf( stderr, "Cannot find target %s\n", h->TrsID );
  return 0;
}

double target_alt_freq( const Hap* h ) {
  size_t target_inx;
  size_t h_inx;
  unsigned int n_ref = 0;
  unsigned int n_alt = 0;
  target_inx = get_target_inx( h );
  for( h_inx = 0; h_inx < h->n_haps; h_inx++ ) {
    if ( h->hap_alleles[h_inx][target_inx] == 0 ) {
      n_ref += h->hap_counts[h_inx];
    }
    else {
      n_alt += h->hap_counts[h_inx];
    }
  }
  if ( (n_ref + n_alt) == 0 ) {
    return 0;
  }
  return (double)n_alt / (n_ref + n_alt);
}

void destroy_Hap( Hap* h ) {
  free( h );
}
