#include <pbc.h>
#include "darray.h"
#include <stdbool.h>

#ifndef __PBC_LPOLY_H__
#define __PBC_LPOLY_H__

typedef struct type_lnodepoly {
    struct type_lnodepoly *next;
    int n;
    int stps;
    element_ptr polynomial;
    element_ptr polynomialY;
} lnodepoly;

typedef struct type_lpoly {
    struct type_lnodepoly *first;
    struct type_lnodepoly *last;
    int MAXD;
} lpoly;

lpoly* lpoly_init();

void lpoly_free(lpoly *list);

void lpoly_add(lpoly *list, element_ptr polynomial, element_ptr polynomialY, int n);

element_ptr get_polynomial(lnodepoly *element);
element_ptr get_polynomialY(lnodepoly *element);
int get_steps(lnodepoly *list);

lnodepoly* first(lpoly *list);
lnodepoly* next(lnodepoly *list);
lnodepoly* last(lpoly *list);
bool islast(lnodepoly *list);

// DE GPT

typedef struct {
  // The coefficients are held in a darray which is resized as needed.
  // The last array entry represents the leading coefficient and should be
  // nonzero. An empty darray represents 0.
  darray_t coeff;
} *peptr;

typedef struct {
  field_ptr field;  // Ring where coefficients live.
  fieldmap mapbase; // Map element from underlying field to constant term.
} *pfptr;

typedef struct {
  field_ptr field;   // Base field.
  fieldmap mapbase;  // Similar to mapbase above.
  int n;             // Degree of extension.
  element_t poly;    // Polynomial of degree n.
  element_t *xpwr;   // x^n,...,x^{2n-2} mod poly
} *mfptr;

void poly_alloc(element_ptr e, int n);
void poly_remove_leading_zeroes(element_ptr e);
void poly_set_coeff(element_ptr e, element_ptr a, int n);
element_ptr poly_get_coeff(element_ptr p, int n);
int polymod_field_degree(field_t f);
#endif
