#ifndef PILEUP_H
#define PILEUP_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "file-io.h"

#define MAX_FIELD_WIDTH (10240)
#define MAX_COV (128)
#define MAX_ID_LEN (256)
// Can be reset once we find the max from actual
// tables

typedef enum { false, true } bool;

typedef struct pul {
  char chr[ MAX_ID_LEN ];
  unsigned int pos;
  char ref;
  unsigned int cov;
  char bases[MAX_COV];
  size_t best_alt_inx;
  unsigned int base_quals[MAX_COV];
  unsigned int map_quals[MAX_COV];
  int strands[MAX_COV];
} Pul;
typedef struct pul* PulP;

typedef struct pu_chr {
  Pul** puls; /* array of Pul* */
  size_t n_puls;        /* number of pileup lines */
  Pul* pul_arr; /* memories where we'll put actually puls */
} Pu_chr;

char revcom_base( const char base );
size_t base_inx( char b );
int best_base_from_pul( PulP pp,
                        unsigned int mqc,
			unsigned int covc, bool weighted );
int rand_good_base_from_pul( PulP pp,
                             unsigned int mqc,
			     unsigned int covc );

int line2pul( char* line, PulP pp );
int base_inx_from_pul( PulP pp,
                       unsigned int mqc,
		       unsigned int covc  );
Pul* fetch_Pul( const Pu_chr* puc, const size_t pos );
static int cmp_Pul( const void *pul1, const void *pul2 );
Pu_chr* init_Pu_chr( const char* fn );

#endif /* PILEUP_H */
