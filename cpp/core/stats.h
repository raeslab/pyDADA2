// Statistical helpers replacing Rcpp/R nmath functions in the ported DADA2 core (spec §5.4).

#ifndef DADA2_STATS_H
#define DADA2_STATS_H

// Upper tail of the Poisson distribution, P(X > n) for X ~ Poisson(lambda).
// Equivalent to R's ppois(n, lambda, lower.tail = FALSE) / Rcpp::ppois(n, lambda, false).
// Returns 1 for n < 0 and 0 for lambda == 0 (n >= 0).
double ppois_upper(int n, double lambda);

#endif  // DADA2_STATS_H
