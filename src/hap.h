#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "file-io.h"

#define MAX_SITES (255)
#define MAX_HAPS (6400)
#define MAX_RSID_LEN (24)
#define MAX_ID_LEN (256)

typedef struct hap {
  char chr[MAX_ID_LEN];
  unsigned int target_pos;
  char Ta0;
  char Ta1;
  char TrsID[MAX_RSID_LEN];
  unsigned int pos[MAX_SITES];
  char a0s[MAX_SITES];
  char a1s[MAX_SITES];
  char rsIDs[MAX_SITES][MAX_RSID_LEN];
  unsigned int hap_counts[MAX_HAPS];
  unsigned int hap_alleles[MAX_HAPS][MAX_SITES]; // 0001011000100, e.g
  size_t n_sites;
  size_t n_haps;
} Hap;

Hap* init_Hap( void );
int read_next_Hap( Hap* h, File_Src* f );
int string2ui( char* string, unsigned int* ha, size_t n );
char* stringify_hap( const Hap* h, const size_t hap_inx, char* s );
unsigned int find_total_pairwise_haps( const Hap* h );
size_t get_target_inx( const Hap* h );
double target_alt_freq( const Hap* h );
void destroy_Hap( Hap* h );
