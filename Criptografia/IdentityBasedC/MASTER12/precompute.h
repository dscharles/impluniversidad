#include <pbc.h>
#include <pbc_poly.h>
#include "lpoly.h"

#ifndef __PBC_PRECOMP_H__
#define __PBC_PRECOMP_H__

typedef struct {
  int inf_flag;    // inf_flag == 1 means O, the point at infinity.
  element_t x, y;  // Otherwise we have the finite point (x, y).
} *point_ptr;

void precompute(lpoly *list, mpz_t q, element_t P, element_t Q);

void change_NMAX(int n);
void print_contador();
void compute_miller(element_ptr out, lpoly *list, element_ptr in2, pairing_t pairing);
void compute_millers(element_ptr out, lpoly *list, element_ptr in2[], int n_prod, pairing_t pairing);

void eval_miller(element_ptr out, element_ptr in1, element_ptr in2,
    pairing_t pairing);

void count_elem(lpoly *list);

#endif
