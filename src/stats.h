#ifndef STATS_H
#define STATS_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "file-io.h"

long unsigned int** init_nchoosek( size_t n );

double binomial_p( unsigned int n, unsigned int k, double p,
		   long unsigned int** nchoosek );
#endif
