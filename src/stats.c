#include "stats.h"

long unsigned int** init_nchoosek( size_t N ) {
  /* n+1 choose k = (n choose k) x (1/k+1) x (n-k) */
  long unsigned int** nchoosek;
  long unsigned int* choosek;
  size_t k, n;

  nchoosek = (long unsigned int**)malloc((N+1)*sizeof(long unsigned int*));  
  choosek = (long unsigned int*)malloc((N+1)*(N+1) * sizeof(long unsigned int));

  for( n = 0; n <= N; n++ ) {
    nchoosek[n] = &choosek[n*(N+1)];
  }

  /* 0 choose anything is 1, for some reason */
  for( k=0; k <= N; k++ ) {
    nchoosek[0][k] = 1;
  }

  for( n = 1; n <= N; n++ ) {
    nchoosek[n][0] = 1;
    nchoosek[n][1] = n;
    for( k = 2; k <= N; k++ ) {
      nchoosek[n][k] =
	nchoosek[n][k-1] * (n-k+1)/k;
    }
  }
  return nchoosek;
}

/* Danger: no range checking for nchoosek!
   (n choose k) x p**k x (1-p)**(n-k)

 */
double binomial_p( unsigned int n, unsigned int k, double p,
		   long unsigned int** nchoosek ) {
  if ( k > n ) {
    fprintf( stderr, "stats.c:binomial_p problem: k > n\n" );
    return 1;
  }
  if ( (p > 1.0) || (p < 0.0) ) {
    fprintf( stderr,
	     "stats.c: binomial_p problem: p not in [0,1]; p = %4.4e\n",
	     p );
    return 1;
  }
  return pow(p, (double)k) * pow(1.0-p, (double)n-k) * (double)nchoosek[n][k];
}
