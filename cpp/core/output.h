// Output builders for dada_uniques (was the Rcpp DataFrame/Matrix builders in DADA2 src/error.cpp).

#ifndef DADA2_OUTPUT_H
#define DADA2_OUTPUT_H

#include "api.h"
#include "dada.h"

void b_make_clustering_df(B *b, Sub **subs, Sub **birth_subs, bool has_quals, dada2::DadaOutput &out);
void b_make_transition_by_quality_matrix(B *b, Sub **subs, bool has_quals, int ncol, dada2::DadaOutput &out);
void b_make_cluster_quality_matrix(B *b, Sub **subs, bool has_quals, unsigned int maxlen, dada2::DadaOutput &out);
void b_make_birth_subs_df(B *b, Sub **birth_subs, bool has_quals, dada2::DadaOutput &out);

#endif  // DADA2_OUTPUT_H
