#include <pbc.h>
#include "lpoly.h"
#include <stdbool.h>

lpoly* lpoly_init(){
    lpoly *list = malloc( sizeof(lpoly) );
    list->MAXD = 0;

    list->first = malloc( sizeof(lnodepoly) );
    list->first->n = -1;
    list->first->next = NULL;

    list->last = list->first;

    return(list);
}

void lpoly_free(lpoly *list){
    lnodepoly *tmp;
    tmp = list->first;
    while(tmp != NULL){
        list->first = list->first->next;
        free(tmp);
        tmp = list->first;
    }
    free(list);
}

void lpoly_free2(lpoly list){
    lnodepoly *tmp;
    tmp = list.first;
    while(tmp != NULL){
        list.first = list.first->next;
        free(tmp);
        tmp = list.first;
    }
}

void lpoly_add(lpoly *list, element_ptr polynomial, element_ptr polynomialY, int n){
    list->last->next = malloc( sizeof(lnodepoly) );
    list->last = list->last->next;

    list->last->next = NULL;
    list->last->polynomial = polynomial;
    list->last->polynomialY = polynomialY;
    list->last->n = n;
    list->last->stps = 1;
}

element_ptr get_polynomial(lnodepoly *element){
    return element->next->polynomial;
}

element_ptr get_polynomialY(lnodepoly *element){
    return element->next->polynomialY;
}

int get_steps(lnodepoly *element){
    return element->next->n;
}

lnodepoly* first(lpoly *list){
    return list->first;
}

lnodepoly* next(lnodepoly *element){
    return element->next;
}

lnodepoly* last(lpoly *list){
    return list->last;
}

bool islast(lnodepoly *list){
    if(list->next == NULL)
    {
        return true;
    }
    return false;
}

// DE GPT

void poly_alloc(element_ptr e, int n) {
  pfptr pdp = e->field->data;
  peptr p = e->data;
  element_ptr e0;
  int k = p->coeff->count;
  while (k < n) {
    e0 = pbc_malloc(sizeof(element_t));
    element_init(e0, pdp->field);
    darray_append(p->coeff, e0);
    k++;
  }
  while (k > n) {
    k--;
    e0 = darray_at(p->coeff, k);
    element_clear(e0);
    pbc_free(e0);
    darray_remove_last(p->coeff);
  }
}

void poly_remove_leading_zeroes(element_ptr e) {
  peptr p = e->data;
  int n = p->coeff->count - 1;
  while (n >= 0) {
    element_ptr e0 = p->coeff->item[n];
    if (!element_is0(e0)) return;
    element_clear(e0);
    pbc_free(e0);
    darray_remove_last(p->coeff);
    n--;
  }
}

void poly_set_coeff(element_ptr e, element_ptr a, int n) {
  peptr p = e->data;
  if (p->coeff->count < n + 1) {
    poly_alloc(e, n + 1);
  }
  element_ptr e0 = p->coeff->item[n];
  element_set(e0, a);
  if (p->coeff->count == n + 1 && element_is0(a)) poly_remove_leading_zeroes(e);
}

element_ptr poly_get_coeff(element_ptr e, int n){
    peptr p = e->data;
    return darray_at(p->coeff, n);
}

int polymod_field_degree(field_t f) {
  mfptr p = f->data;
  return p->n;
}
