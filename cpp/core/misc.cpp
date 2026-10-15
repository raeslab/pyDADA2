// Adapted from DADA2 src/misc.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * Rprintf/Rcpp::stop -> dada2_log/dada2_stop.

#include "dada.h"

void err_print(double err[4][4]) {
  for(int i=0;i<4;i++) {
    if(i==0)  { dada2_log("{"); }
    else      { dada2_log(" "); }
    dada2_log("{");
    for(int j=0;j<4;j++) {
      dada2_log("%.6f", err[i][j]);
      if(j<3) { dada2_log(", "); }
    }
    if(i<3) { dada2_log("},\n"); }
    else    { dada2_log("}}\n"); }
  }
}

void align_print(char **al) {
  char *al0, *al1;
  if(!al) { return; }
  al0 = al[0]; al1 = al[1];
  
  dada2_log("%s\n", ntstr(al0));
  for(int i=0;i<strlen(al0);i++) {
    if(al0[i] == al1[i]) { dada2_log("|"); }
    else { dada2_log(" "); }
  }
  dada2_log("\n%s\n", ntstr(al1));
  
}

/* nt2int takes an input string (iseq) and converts it to numeric indices,
  which are stored in the char output array (oseq).
 Need memory at oseq equal to strlen of iseq, which is not checked currently.
  nt2int(seq,seq) gives previous functionality.
*/

void nt2int(char *oseq, const char *iseq) {
  int i, len = strlen(iseq);

  for (i = 0; i < len; i++, iseq++, oseq++) {
    switch (*iseq) {
    case 'A':
      *oseq = 1;
      break;
    case 'C':
      *oseq = 2;
      break;
    case 'G':
      *oseq = 3;
      break;
    case 'T':
      *oseq = 4;
      break;
    case 'N':
      *oseq = 5;
      break;
    case '-':
      *oseq = '-';
      break;
    default:
      dada2_log("invalid character in input:%c.\n",*iseq);
      *oseq = '\0';
    }
  }
  *oseq = '\0';
  return;
}

/* Converts seq in index form back to nts. */
void int2nt(char *oseq, const char *iseq) {
  int i, len = strlen(iseq);
  for (i = 0; i < len; i++, iseq++, oseq++) {
    switch (*iseq) {
    case 1:
      *oseq = 'A';
      break;
    case 2:
      *oseq = 'C';
      break;
    case 3:
      *oseq = 'G';
      break;
    case 4:
      *oseq = 'T';
      break;
    case 5:
      *oseq = 'N';
      break;
    case '-':
      *oseq = '-';
      break;
    default:
      break;
    }
  }
  *oseq = '\0';
  return;
}

/* Convenience function for diagnostic output. */
void ntcpy(char *oseq, const char *iseq) {
  strcpy(oseq, iseq);
  int2nt(oseq, oseq);
}

/* Convenience function for diagnostic output. */
char *ntstr(const char *iseq) {
  char *oseq = (char *) malloc(strlen(iseq)+1); //E
  if (oseq == NULL)  dada2_stop("Memory allocation failed!\n");
  
  ntcpy(oseq, iseq);
  return oseq;
}

/* Convenience function for diagnostic input. */
char *intstr(const char *iseq) {
  char *oseq = (char *) malloc(strlen(iseq)+1); //E
  if (oseq == NULL)  dada2_stop("Memory allocation failed!\n");
  
  strcpy(oseq, iseq);
  nt2int(oseq, oseq);
  return oseq;
}
