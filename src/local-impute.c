#include "local-impute.h"

/* (c) 2019 - Ed Green
   Astrea Forensics
*/

Impute_Set* init_Impute_Set( void ) {
  Impute_Set* isp;
  isp = (Impute_Set*)malloc(sizeof(Impute_Set));
  isp->epsilon = (double)ERROR_PER_THOUSAND / 1000;
  isp->chr[0] = '\0';
  isp->chr[1] = '\0';
  isp->chr[2] = '\0';
  isp->chr[3] = '\0';
  isp->chr[4] = '\0';
  
  isp->T_rsID[0] = '\0';
  isp->T0 = 0;
  isp->T1 = 0;
  isp->num_pairs = 0;
  isp->impute_pairs[0].L0 = 0;
  isp->impute_pairs[0].L1 = 0;
  return isp;
}

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
			      short unsigned int t1 ) {
  if (isp == NULL) {
    return -1;
  }
  if ( strlen( T_rsID ) <= MAX_ID_LEN ) {
    strcpy( isp->T_rsID, T_rsID );
  }
  else {
    strncpy( isp->T_rsID, T_rsID, MAX_ID_LEN );
    isp->T_rsID[MAX_ID_LEN] = '\0';
  }
  strcpy( isp->chr, chr );
  isp->T_pos = T_pos;
  isp->T0 = t0;
  isp->T1 = t1;
  return 0;
}

void clear_Impute_Set( Impute_Set* isp ) {
  isp->num_pairs = 0;
  return;
}

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
		     char la0, // reference allele
		     char la1, // alternate allele
		     short unsigned int t0l0,
		     short unsigned int t0l1,
		     short unsigned int t1l0,
		     short unsigned int t1l1,
		     double min_D,
		     int Target
		     ) {
  Impute_Pair* imp;

  /* If this Impute_Set is full, don't add more.
     Just return -1 */
  if( isp->num_pairs >= MAX_PAIRS ) {
    return -1;
  }

  /* There's room. Get a pointer to the next one to fill it up */
  imp = &(isp->impute_pairs[isp->num_pairs]);

  /* Add the data for this pair */
  strcpy( imp->chr, chr );
  strcpy( imp->L_rsID, L_rsID );
  imp->L_pos = L_pos;
  imp->la0   = la0;
  imp->la1   = la1;
  imp->T0L0  = t0l0;
  imp->T0L1  = t0l1;
  imp->T1L0  = t1l0;
  imp->T1L1  = t1l1;

  /* Add this if it's a Target site or a Linked site with 
     >= min_D LD to Target */
  if ( Target ||
       (calc_D( imp ) >= min_D) ) {
    isp->num_pairs++;
    return 0;
  }
  else {
    return -2;
  }
}

/* Returns the LD value (D) given the
   input data in Impute_Pair* 
   D = pAB - (pA * pB)
*/
double calc_D( Impute_Pair* imp ) {
  double D;
  /* D = pAB - (pA * pB) */
  D = ((double)imp->T0L0 / (imp->T0L0 + imp->T0L1 +
			    imp->T1L0 + imp->T1L1))
    -
    ( (double)(imp->T0L0 + imp->T0L1) / (imp->T0L0 +
					 imp->T0L1 +
					 imp->T1L0 +
					 imp->T1L1)
      *
      (double)(imp->T0L0 + imp->T1L0) / (imp->T0L0 +
					 imp->T0L1 +
					 imp->T1L0 +
					 imp->T1L1) );

  return fabs( D );
}

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
			short unsigned int l1 ) {
  if ( imp != NULL ) {
    imp->L0 = l0;
    imp->L1 = l1;
    return;
  }
  return;
}

/* Returns probability of all Linked sites data given
   that the Target genotype is 0, 0 (and the haplotype
   frequencies in this Impute_Set */
double p_allLgivenT00( Impute_Set* isp ) {
  size_t i;
  double p = 1.0;
  for( i = 0; i < isp->num_pairs; i++ ) {
    p *= p_LgivenT00( &isp->impute_pairs[i], isp->epsilon );
  }
  return p;
}

double p_allLgivenT01( Impute_Set* isp ) {
  size_t i;
  double p = 1.0;
  for( i = 0; i < isp->num_pairs; i++ ) {
    p *= p_LgivenT01( &isp->impute_pairs[i], isp->epsilon );
  }
  return p;
}

double p_allLgivenT11( Impute_Set* isp ) {
  size_t i;
  double p = 1.0;
  for( i = 0; i < isp->num_pairs; i++ ) {
    p *= p_LgivenT11( &isp->impute_pairs[i], isp->epsilon );
  }
  return p;
}

/* Find the probability of the observed alleles at this Linked site (L),
   given that the genotype at T is 0/0 and the haplotype frequencies in 
   the panel as described in the input Impute_Pair and the error rate
   epsilon
*/
double p_LgivenT00( const Impute_Pair* imp, const double epsilon ) {
  double pL0;
  /* It's weird, but sometimes the reference allele is *never* seen
     it the 1KG (or other panel) data */
  if ( (imp->T0L0 + imp->T0L1) == 0 ) {
    pL0 = 0.0;
  }
  else {
    /* Find Probability L=0 | T=0,0 from the haplotype frequencies */
    pL0 =
      (double)imp->T0L0
      /
      (double)(imp->T0L0 + imp->T0L1);
  }
  
  /* Now, adjust the pL0 by the error rate, epsilon */
  pL0 =
    (pL0 * (1 - epsilon))  // probability of it being 0 & no error
    +                      // plus
    ((1 - pL0) * epsilon); // probability of it not being 0 but error
  
  /* Now, the probability of the given number of L=0 (successes)
     given the pL0 (which comes from the haplotype frequencies)
     can be calculated from the binomial distribution:
     p(k, n, p) = n choose k * p**k * (1-p)**(n-k)
  */
  return ((double)NchooseK((imp->L0 + imp->L1), imp->L0, fact) *
	  pow(pL0, imp->L0) *
	  pow(1-pL0, imp->L1));
}

/* */
double p_LgivenT01( const Impute_Pair* imp, const double epsilon ) {
  /* It's weird, but sometimes the reference allele is *never* seen
     it the 1KG (or other panel) data */
  double pL0;
  if ( ((imp->T0L0 + imp->T0L1) > 0) &&
       ((imp->T1L0 + imp->T1L1) > 0) ) {
    /* Saw both alleles at least once */
    /* Probability of L=0 | T=0,1 
       is half the chance of seeing L0 if T=0 &
       half the chance of seeing L0 if T=1 */
    pL0 =
      (0.5 * (double)imp->T0L0/(imp->T0L0 + imp->T0L1)) +
      (0.5 * (double)imp->T1L0/(imp->T1L0 + imp->T1L1));    
  }
  else { /* either ref or alt is all the alleles! WTF? Target
	    site is not polymorphic in panel.
	    rs104893968, e.g., on chr6 has 0 ALT allele in 1KG
	    In this case, pL0 | T01 must be 0.5 */
    pL0 = 0.5;
  }
  
  /* Adjust pL0 by the error rate, epsilon */
  pL0 = (pL0 * (1-epsilon)) + ((1-pL0) * epsilon);
  
  return ((double)NchooseK((imp->L0 + imp->L1), imp->L0, fact) *
	  pow(pL0, imp->L0) *
	  pow(1-pL0, imp->L1));
}

/* Like p_LgivenT00, but for the T=1,1 genotype instead
 of the T=0,0 genotype
*/
double p_LgivenT11( const Impute_Pair* imp, const double epsilon ) {
  double pL0;
  if ( (imp->T1L0 + imp->T1L1) == 0 ) {
    pL0 = 0.0;
  }
  else {
    /* Find Probability L=0 | T=1,1 from the haplotype frequencies:
       T=0, L=0 / (T=0, L=0 + T=0, L=1) */
    pL0 = (double)imp->T1L0/(double)(imp->T1L0 + imp->T1L1);
  }
  
  /* Adjust pL0 by the error rate, epsilon */
  pL0 = (pL0 * (1-epsilon)) + ((1-pL0) * epsilon);
  
  /* Now, the probability of the given number of L=0 (successes)
     given the pL0 (which comes from the haplotype frequencies)
     can be calculated from the binomial distribution:
     p(k, n, p) = n choose k * p**k * (1-p)**(n-k)
  */
  return ((double)NchooseK((imp->L0 + imp->L1), imp->L0, fact) *
	  pow(pL0, imp->L0) *
	  pow(1-pL0, imp->L1));
}

/* NchooseK 
   n choose k = n! / k!(n-k)!

   If input n or k are > MAX_ALLELE_COUNTS, crudely down-sample
   to avoid catching on fire
*/
unsigned long long NchooseK( unsigned n, unsigned k,
			     unsigned long long* factorial ) {

  double d_n, d_k;
  if ( (n > MAX_ALLELE_COUNTS) ||
       (k > MAX_ALLELE_COUNTS) ) {
    d_n = (double)n;
    d_k = (double)k;
    d_k = (d_k / d_n ) * MAX_ALLELE_COUNTS;
    d_n = (double)MAX_ALLELE_COUNTS;

    n = (unsigned)round(d_n);
    k = (unsigned)round(d_k);
  }

  return ( factorial[n]
	   /
	   ( factorial[k] * factorial[n-k] ) ); 
}

unsigned long long* init_factorial( void ) {
  unsigned long long* f;
  size_t i;
  f = (unsigned long long*)malloc(sizeof(unsigned long long) *
				  (MAX_ALLELE_COUNTS +1));
  f[0] = 1;
  for( i = 1; i <= MAX_ALLELE_COUNTS; i++ ) {
    f[i] = f[i-1] * i;
  }
  return f;
}

/* Returns the probability of the target genotype being 0,0 
   given the reference panel allele frequencies at T
 */
double p_T00( Impute_Set* isp ) {
  return (pow( (double)isp->T0 / (isp->T0 + isp->T1),
	       2 ));
}

/* Returns the probability of the target genotype being 0, 1
   i.e., heterozygous, given the reference panel allele
   frequencies at T => 2pq
   Hardy Freakin Weinberg
*/
double p_T01( Impute_Set* isp ) {
  return ( 2 *
	   (double)isp->T0 / (isp->T0 + isp->T1) *
	   (double)isp->T1 / (isp->T0 + isp->T1) );
}

/* Returns the probability of the target genotype being 1,1
   given the reference panel allele frequencies at T
 */
double p_T11( Impute_Set* isp ) {
  return (pow( (double)isp->T1 / (isp->T0 + isp->T1),
	       2 ));
}

int in_sex_haploid( Impute_Set* isp ) {
  if ( (isp->chr[0] == 'X') ||
       (isp->chr[4] == 'X') ) { // on X
    if ( (isp->T_pos >= XPAR1ST) &&
	 (isp->T_pos <= XPAR1EN) ) {
      return 0;
    }
    if ( (isp->T_pos >= XPAR2ST) &&
	 (isp->T_pos <= XPAR2EN) ) {
      return 0;
    }
    return 1; // on X but not in either PAR
  }
  if ( (isp->chr[0] == 'Y') ||
       (isp->chr[4] == 'Y') ) {
    if ( (isp->T_pos >= YPARST) &&
	 (isp->T_pos <= YPAREN) ) {
      return 0;
    }
    return 1;
  }
  return 0;
}
  
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
              (5) char* a0 place to put the first allele call
              (6) char* a1 place to put the second allele call
	      (7) int use_priors - if true, use the genotype priors for
                  calling likelihoods
   Returns: nothing
   Populates the chars pointed to by a0 and a1 with the allele calls
*/
void isp2genotype( Impute_Set* isp, int Male, double min_ratio,
		   double min_ratio_zero_cov,
		   double min_ratio_target_only,
		   double special_hom, char* a0, char* a1,
		   int use_priors ) {
  double pHomRef, pHet, pHomAlt;
  double pBest;
  double pSecond;
  double ratio_to_use;
  char ref;
  char alt;
  ref = isp->impute_pairs[0].la0; // first "pair" is always itself
  alt = isp->impute_pairs[0].la1; // first "pair" is always itself

  /* First check if the answer is clear already from the Target
     site data alone. If so, no need to look at linked sites
     For this, use the min_ratio_target_only
 */
  pHomRef = p_LgivenT00( &isp->impute_pairs[0], isp->epsilon );
  if ( Male && in_sex_haploid( isp ) ) {
    pHet = 0.0;
  }
  else {
    pHet    = p_LgivenT01( &isp->impute_pairs[0], isp->epsilon );
  }
  pHomAlt = p_LgivenT11( &isp->impute_pairs[0], isp->epsilon );

  if ( use_priors ) {
    pHomRef *= p_T00(isp);
    pHet    *= p_T01(isp);
    pHomAlt *= p_T11(isp);
  }
  
  /* First check if the Target site data alone make a clear answer
     If so, just set the alleles and be done */
  if ( unambiguous_genotype( pHomRef, pHet, pHomAlt,
			     min_ratio_target_only, ref, alt,
			     a0, a1 ) ) {
    return;
  }
  /* No clear answer from Target Site only data.
     Use all the linked sites data, too */
  else {
    pHomRef = p_allLgivenT00(isp);
    if ( Male && in_sex_haploid( isp ) ) {
      pHet = 0.0;
    }
    else {
      pHet  = p_allLgivenT01(isp);
    }
    pHomAlt = p_allLgivenT11(isp);
  }

  if ( use_priors ) {
    pHomRef *= p_T00(isp);
    pHet    *= p_T01(isp);
    pHomAlt *= p_T11(isp);
  }
  
  /* pHomRef >= pHet >= pHomAlt */
  if ( (pHomRef >= pHet) && (pHet >= pHomAlt) ) {
    pBest   = pHomRef;
    pSecond = pHet;
    *a0 = ref;
    *a1 = ref;
  }
  
  /* pHomRef >= pHomAlt >= pHet */
  if ( (pHomRef >= pHomAlt) && (pHomAlt >= pHet) ) {
    pBest   = pHomRef;
    pSecond = pHomAlt;
    *a0 = ref;
    *a1 = ref;
  }

  /* pHet >= pHomRef >= pHomAlt */
  if ( (pHet >= pHomRef) && (pHomRef >= pHomAlt) ) {
    pBest   = pHet;
    pSecond = pHomRef;
    *a0 = ref;
    *a1 = alt;
  }
  
  /* pHet >= pHomAlt >= pHomRef */
  if ( (pHet >= pHomAlt) && (pHomAlt >= pHomRef) ) {
    pBest   = pHet;
    pSecond = pHomAlt;
    *a0 = alt;
    *a1 = ref;
  }
  
  /* pHomAlt >= pHet >= pHomRef */
  if ( (pHomAlt >= pHet) && (pHet >= pHomRef) ) {
    pBest   = pHomAlt;
    pSecond = pHet;
    *a0 = alt;
    *a1 = alt;
  }
  
  /* pHomAlt >= pHomRef >= pHet*/
  /* Not sure the can actually happen, but check anyway */
  if ( (pHomAlt >= pHomRef) && (pHomRef >= pHet) ) {
    pBest   = pHomAlt;
    pSecond = pHomRef;
    *a0 = alt;
    *a1 = alt;
  }

  /* Now, *a0 and *a1 are set with most probably alleles, in order */
  
  /* Is this a site with no coverage at the Target site?
     If so, use the min_ratio_zero_cov */
  if ( (isp->impute_pairs[0].L0 == 0) &&
       (isp->impute_pairs[0].L1 == 0) ) {
    ratio_to_use = min_ratio_zero_cov;
  }
  else {
    ratio_to_use = min_ratio;
  }
  
  /* Check the ratio of the best genotype versus the second best.
     If the ratio is >= min_ratio, then we call it. Don't compare
     if the pSecond == 0! */
  if ( pSecond > 0 ) {
    if ( (pBest / pSecond) >= ratio_to_use ) {
      ; // we're done, fall through
    }
    else { /* best isn't distinguishable from second best.
      /* special_hom is set? */
      if ( special_hom > 0 ) {
	if ( (pHomRef >= pHet) && (pHet >= pHomAlt) &&
	     (isp->impute_pairs[0].L0 > 0) &&
	     (isp->impute_pairs[0].L1 == 0) &&
	     ((pHet / pHomAlt) > special_hom) &&
	     (pBest / pSecond) > 10 ) {
	  /* homozygous ref by -v option rule, alleles are set already */
	  return ;
	}
	if ( (pHomAlt >= pHet) && (pHet >= pHomRef) &&
	     (isp->impute_pairs[0].L0 == 0) &&
	     (isp->impute_pairs[0].L1 > 0) &&
	     ((pHet/pHomRef) > special_hom) &&
	     (pBest / pSecond) > 10 ) {
	  /* homozygous alt by -v option rule, alleles are set already */
	  return;
	}
      }
      /* either didn't specify special_hom or didn't get used
	 , maybe a one-call */
      *a1 = '0';
      /* First check Target site. If there is data for one allele
	 and not the other, then use that */
      if ( (isp->impute_pairs[0].L0 > 0) &&
	   (isp->impute_pairs[0].L1 == 0) ) {
	*a0 = ref;
      }
      if ( (isp->impute_pairs[0].L1 > 0) &&
	   (isp->impute_pairs[0].L0 == 0) ) {
	*a0 = alt;
      }
      if ( (isp->impute_pairs[0].L0 > 0) &&
	   (isp->impute_pairs[0].L1 > 0) &&
	   (pHet > pHomRef) &&
	   (pHet > pHomAlt) ) {
	/* Both Target site alleles observed and pHet is most likely,
	   call it a het regardless */
	*a0 = ref;
	*a1 = alt;
      }
      /* Must have been a oo coverage site? 
	 We'll just live with the half-call */
    }
  }
  else { // pSecond == 0
    ; // we're done, fall through
  }

  /* Check to make sure we're not calling a homozygote of one allele
     but SAW the other allele! This can happen in really weird cases.
     Change homozygous to heterozygous
  */
  if ( (*a0 == alt) &&
       (*a1 == alt) ) { // If this is gonna be called homozygous alt
    if ( isp->impute_pairs[0].L0 > 0 ) { // and we SAW ref alleles!
      *a0 = ref; // set one allele to ref
    }
  }
  if ( (*a0 == ref) &&
       (*a1 == ref) ) { // If this is gonna be called homozygous ref
    if ( isp->impute_pairs[0].L1 > 0 ) { // and we SAW alt alleles!
      *a1 = alt; // set one allele to alt
    }
  }

  /*
    For non-autosomal regions, low confidence or het read data could
    result in heterozygous calls despite pHet == 0.0. So we will do a
    no-call if this happens
   */
  if (pHet == 0.0 && (a1 != a0)) {
    *a0 = '0';
    *a1 = '0';
  }

  return;
}

/* unambiguous_genotype
   Args: probabilities of Homozygous Ref, 
                          Heterozygous,
                          Homozygous Alternate,
                          minimum ratio to be unambiguous
   Returns: True IFF the ratio of the probability of the best
            genotype >= min_ratio the probability of the second
            best genotype
*/
int unambiguous_genotype( const double pHomRef, const double pHet,
			  const double pHomAlt, const double min_ratio,
			  const char ref, const char alt,
			  char* a0, char* a1 ) {

  double pBest, pSecond;
  /* pHomRef >= pHet >= pHomAlt */
  if ( (pHomRef >= pHet) && (pHet >= pHomAlt) ) {
    pBest   = pHomRef;
    pSecond = pHet;
    *a0 = ref;
    *a1 = ref;
  }
  
  /* pHomRef >= pHomAlt >= pHet */
  if ( (pHomRef >= pHomAlt) && (pHomAlt >= pHet) ) {
    pBest   = pHomRef;
    pSecond = pHomAlt;
    *a0 = ref;
    *a1 = ref;
  }

  /* pHet >= pHomRef >= pHomAlt */
  if ( (pHet >= pHomRef) && (pHomRef >= pHomAlt) ) {
    pBest   = pHet;
    pSecond = pHomRef;
    *a0 = ref;
    *a1 = alt;
  }
  
  /* pHet >= pHomAlt >= pHomRef */
  if ( (pHet >= pHomAlt) && (pHomAlt >= pHomRef) ) {
    pBest   = pHet;
    pSecond = pHomAlt;
    *a0 = ref;
    *a1 = alt;
  }
  
  /* pHomAlt >= pHet >= pHomRef */
  if ( (pHomAlt >= pHet) && (pHet >= pHomRef) ) {
    pBest   = pHomAlt;
    pSecond = pHet;
    *a0 = alt;
    *a1 = alt;
  }
  
  /* pHomAlt >= pHomRef >= pHet*/
  /* Not sure the can actually happen, but check anyway */
  if ( (pHomAlt >= pHomRef) && (pHomRef >= pHet) ) {
    pBest   = pHomAlt;
    pSecond = pHomRef;
    *a0 = alt;
    *a1 = alt;
  }

  if ( pSecond > 0 ) {
    if ( (pBest / pSecond) >= min_ratio ) {
      return 1;
    }
    else {
      return 0;
    }
  }
  else { // pSecond is 0
    if (pBest > 0) {
      return 1;
    }
  }
  return 0;
}

Blacklist* read_blacklist( const char* bl_fn ) {
  Blacklist* BL;
  size_t i;
  FILE* bl_f;
  char blacklist_line[ MAX_ID_LEN + 1 ];
  size_t string_len;
  /* Initialize the blacklist */
  BL = (Blacklist*)malloc(sizeof(Blacklist));
  BL->bl_ids = (char**)malloc(sizeof(char*) * MAX_BLACKLIST_LEN);
  for( i = 0; i < MAX_BLACKLIST_LEN; i++ ) {
    BL->bl_ids[i] = &BL->bl_ids_array[i * MAX_ID_LEN ];
  }
  BL->num_bl_ids = 0;
  
  /* Read in the blacklist IDs */
  bl_f = fileOpen( bl_fn, "r" );
  if ( bl_f == NULL ) {
    fprintf( stderr, "Problem reading blacklist %s\n", bl_fn );
    return NULL;
  }

  while( fgets( blacklist_line, MAX_ID_LEN, bl_f ) != NULL ) {
    // need to strip newline
    string_len = strlen( blacklist_line );
    if ( blacklist_line[ string_len - 1 ] == '\n' ) {
      blacklist_line[ string_len - 1 ] = '\0';
    }
    strcpy( BL->bl_ids[ BL->num_bl_ids ], blacklist_line );
    BL->num_bl_ids++;
  }
  fclose( bl_f );

  /* Sort the blacklist */
  qsort( &BL->bl_ids[0], BL->num_bl_ids, sizeof(char*), cmpstringp );
  
  return BL;
}

/* Return true if not in blacklist
   return false (0) if in blacklist */
int black_check( Blacklist* BL, const char* rsID ) {
  if (bsearch( &rsID, BL->bl_ids,
	       BL->num_bl_ids, sizeof(char*), cmpstringp )
      == NULL ) {
    return 1;
  }
  return 0;
}

static int cmpstringp(const void *p1, const void *p2) {
  /* The actual arguments to this function are "pointers to
     pointers to char", but strcmp(3) arguments are "pointers
     to char", hence the following cast plus dereference */
  
  return strcmp(* (char * const *) p1, * (char * const *) p2);
}
int idCmp(const void* id1_, const void* id2_) {
        char** id1p = (char**) id1_;
        char** id2p = (char**) id2_;
        char* id1 = *id1p;
        char* id2 = *id2p;
        return strcmp(id1, id2);
}
