#include <pbc.h>
#include <pbc_poly.h>
#include "darray.h"
#include "precompute.h"
#include "lpoly.h"
#include <stdbool.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*************************************************************************
/ Fixed argument optimization
/ The first argument of e(P,Q) is fixed, so we can precompute
/ a lot of info.
/
/ Implementation for d-type curves
/ P in E(F_q)
/ Q in E(F_qd) --> E(F_q2d)
/ TWIST: (x,y) --> ( xv^{-1}, xv^{-2} sqrt{v} )
/
/ INFO: When computing g_n(X, Y) all X^i lie in F_qk (with Im = 0)
/ and all X^iY lie in F_qk (with Re = 0)
**************************************************************************/

element_ptr curve_equation;
int STEPMAX = 4;

void change_NMAX(int n){
    STEPMAX = n;
}

int mult1, add1, doublepoint, addpoint, multd, multk, addk, squarek;

void init_contador(){
    mult1 = 0;
    add1 = 0;
    doublepoint = 0;
    addpoint = 0;
    multd = 0;
    multk = 0;
    addk = 0;
    squarek = 0;
}

void print_contador(){
    printf("Mult Fp: %i \n", mult1);
    printf("Add Fp: %i \n", add1);
    printf("Double Point: %i \n", doublepoint);
    printf("Add Point: %i \n", addpoint);
    printf("Mult Fpd*Fpk: %i \n", multd);
    printf("Mult Fpk: %i \n", multk);
    printf("Add Fpk: %i \n", addk);
    printf("Square Fpk: %i \n", squarek);
}

void KSET0(element_t out){
    element_set0(out);
    element_ptr re_out = element_x(out);
    element_set0(element_item(re_out,0));
}

struct timeval tvvBegin, tvvEnd;
void timeval_sub(struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);

    printf("%lu\n", diff);
}



/*****************************
/  DATA STRUCTURE OF D-curves
*****************************/

typedef struct {
  field_ptr field; // The field where the curve is defined.
  element_t a, b;  // The curve is E: Y^2 = X^3 + a X + b.
  // cofac == NULL means we're using the whole group of points.
  // otherwise we're working in the subgroup of order #E / cofac,
  // where #E is the number of points in E.
  mpz_ptr cofac;
  // A generator of E.
  element_t gen_no_cofac;
  // A generator of the subgroup.
  element_t gen;
  // A non-NULL quotient_cmp means we are working with the quotient group of
  // order #E / quotient_cmp, and the points are actually coset
  // representatives. Thus for a comparison, we must multiply by quotient_cmp
  // before comparing.
  mpz_ptr quotient_cmp;
} *curve_data_ptr;


typedef struct {
  field_t Fq, Fqx, Fqd, Fqk;  // The fields F_q, F_q[x], F_q^d, F_q^k.
  field_t Eq, Etwist;         // The curves E(F_q) and E'(F_q^d).
  // Let v be the quadratic nonresidue used to construct F_q^k from F_q^d,
  // namely Fqk = Fqd[sqrt(v)].
  element_t nqrinv, nqrinv2;  // The constants v^-1 and v^-2.
  mpz_t tateexp;              // The Tate exponent,
                              // to standardize coset representatives.
  int k;                      // The embedding degree, usually 6.
  // Let x be the element used to build Fqd from Fq, i.e. Fqd = Fq[x].
  element_t xpowq, xpowq2;    // x^q and x^{2q} in F_q^d.
} *pptr;


inline element_ptr curve_a_coeff(element_t e) {
  return ((curve_data_ptr) e->field->data)->a;
}

inline element_ptr curve_b_coeff(element_t e) {
  return ((curve_data_ptr) e->field->data)->b;
}

element_ptr curve_x_coord(element_t P){
    point_ptr tmp = P->data;
    return tmp->x;
}

element_ptr curve_y_coord(element_t P){
    point_ptr tmp = P->data;
    return tmp->y;
}









/**************************
/     ORIGINAL MILLER
**************************/

void d_miller_evalfn(element_t e0,
    element_t a, element_t b, element_t c, element_t Qx, element_t Qy) {
  element_ptr re_out = element_x(e0);
  element_ptr im_out = element_y(e0);

  int i;
  int d = polymod_field_degree(re_out->field);
  for (i = 0; i < d; i++) {
    element_mul(element_item(re_out, i), element_item(Qx, i), a); mult1++;
    element_mul(element_item(im_out, i), element_item(Qy, i), b); mult1++;
  }
  element_add(element_item(re_out, 0), element_item(re_out, 0), c); add1++;
}

void miller(element_t res, mpz_t q, element_t P, element_ptr Qx, element_ptr Qy) {
  int m;
  element_t v;
  element_t Z;
  element_t a, b, c;
  element_t t0;
  element_t e0;
  const element_ptr cca = curve_a_coeff(P);
  const element_ptr Px = curve_x_coord(P);
  const element_ptr Py = curve_y_coord(P);
  element_ptr Zx, Zy;

  void do_tangent(void) {
    // a = -(3 Zx^2 + cc->a)
    // b = 2 * Zy
    // c = -(2 Zy^2 + a Zx);

    element_square(a, Zx); mult1++;
    element_mul_si(a, a, 3); add1++; add1++; add1++;
    element_add(a, a, cca); add1++;
    element_neg(a, a);

    element_add(b, Zy, Zy); add1++;

    element_mul(t0, b, Zy); mult1++;
    element_mul(c, a, Zx); mult1++;
    element_add(c, c, t0); add1++;
    element_neg(c, c);

    d_miller_evalfn(e0, a, b, c, Qx, Qy);
    element_mul(v, v, e0); multk++;
  }

  void do_line(void) {
    // a = -(B.y - A.y) / (B.x - A.x);
    // b = 1;
    // c = -(A.y + a * A.x);
    // but we multiply by B.x - A.x to avoid division.

    element_sub(b, Px, Zx); add1++;
    element_sub(a, Zy, Py); add1++;
    element_mul(t0, b, Zy); mult1++;
    element_mul(c, a, Zx); mult1++;
    element_add(c, c, t0); add1++;
    element_neg(c, c);

    d_miller_evalfn(e0, a, b, c, Qx, Qy);
    element_mul(v, v, e0); multk++;
  }

  element_init(a, Px->field);
  element_init(b, a->field);
  element_init(c, a->field);
  element_init(t0, a->field);
  element_init(e0, res->field);

  element_init(v, res->field);
  element_init(Z, P->field);

  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_y_coord(Z);

  element_set1(v);
  m = mpz_sizeinbase(q, 2) - 2;

  for(;;) {
    do_tangent();

    if (!m) break;

    element_double(Z, Z); doublepoint++;

    if (mpz_tstbit(q, m)) {
      do_line();

      element_add(Z, Z, P); addpoint++;
    }
    m--;

    element_square(v, v); squarek++;
  }

  element_set(res, v);

  element_clear(v);
  element_clear(Z);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(t0);
  element_clear(e0);

}











/******************************
/       PRECOMP ALGORITHM
******************************/

// Add a new step in the data structure
//
// list -> list of polynomials
// a, b, c -> coeff of the new step: aX + bY + C
// aa -> 1 if the new step comes from Doubling; 0 from the Add part.
//

void add_poly(lpoly *list, element_t a, element_t b, element_t c, int aa){
    // Init the new step
    field_ptr tbase = pbc_malloc(sizeof(*tbase));
    field_init_poly(tbase, a->field);

    element_ptr ppx, ppy;
    ppx = pbc_malloc( sizeof(*ppx));
    ppy = pbc_malloc( sizeof(*ppy));
    element_init(ppx, tbase);
    element_init(ppy, tbase);

    // p = ax + c + y(b)
    poly_set_coeff(ppy, b, 0);
    poly_set_coeff(ppx, c, 0);
    poly_set_coeff(ppx, a, 1);

    if(last(list)->n != -1 && (aa == 0 || STEPMAX > last(list)->stps)){
        // If there is some step in the list adn we can make it bigger

        // INIT

        element_ptr polyx = last(list)->polynomial;
        element_ptr polyy = last(list)->polynomialY;
        /*printf("\n MODIFICANDO \n");
        printf("POLY X: "); element_out_str(stdout, 10, polyx);
        printf("\nPOLY Y: "); element_out_str(stdout, 10, polyy);
        printf("\ncon cg: "); element_out_str(stdout, 10, curve_equation);
        printf("\npx : "); element_out_str(stdout, 10, ppx);
        printf("\npy : "); element_out_str(stdout, 10, ppy);*/

        element_ptr newx, newy, tmp1, tmp2;
        newx = pbc_malloc( sizeof(*newx));
        newy = pbc_malloc( sizeof(*newy));
        tmp1 = pbc_malloc( sizeof(*tmp1));
        tmp2 = pbc_malloc( sizeof(*tmp2));
        element_init(newx, tbase);
        element_init(newy, tbase);
        element_init(tmp1, tbase);
        element_init(tmp2, tbase);
        element_set0(newx); element_set0(newy); element_set0(tmp1); element_set0(tmp2);

        if(aa == 1){
            // If Doubling, we should do the square

            element_mul(tmp1, polyx, polyy);
            element_double(tmp1, tmp1);
            element_mul(tmp2, tmp1, ppy);
            element_mul(tmp2, tmp2, curve_equation);
            element_add(newx, newx, tmp2);
            element_mul(tmp2, tmp1, ppx);
            element_add(newy, newy, tmp2);

            element_square(tmp1, polyy);
            element_mul(tmp2, tmp1, ppy);
            element_mul(tmp2, tmp2, curve_equation);
            element_add(newy, newy, tmp2);
            element_mul(tmp2, tmp1, ppx);
            element_mul(tmp2, tmp2, curve_equation);
            element_add(newx, newx, tmp2);

            element_square(tmp1, polyx);
            element_mul(tmp2, tmp1, ppy);
            element_add(newy, newy, tmp2);
            element_mul(tmp2, tmp1, ppx);
            element_add(newx, newx, tmp2);
        }
        else{
            // If Add, we only need to multiply the new step

            element_mul(tmp1, polyx, ppx);
            element_add(newx, newx, tmp1);
            element_mul(tmp1, polyy, ppy);
            element_mul(tmp1, tmp1, curve_equation);
            element_add(newx, newx, tmp1);

            element_mul(tmp1, polyx, ppy);
            element_add(newy, newy, tmp1);
            element_mul(tmp1, polyy, ppx);
            element_add(newy, newy, tmp1);
        }

        //printf("\n RESULTADO X : "); element_out_str(stdout, 10, newx);
        //printf("\n RESULTADO Y : "); element_out_str(stdout, 10, newy); printf("\n");


        last(list)->polynomial = newx;
        last(list)->polynomialY = newy;

        last(list)->n = last(list)->n + aa;
        last(list)->stps = last(list)->stps + 1;

        // Update the maximum degree, later on we need it to precompute some elements.
        int d = poly_degree(last(list)->polynomial);
        if(d > list->MAXD){
            list->MAXD = d;
        }
        d = poly_degree(last(list)->polynomialY);
        if(d > list->MAXD){
            list->MAXD = d;
        }
    }else{
        // First step or previous step is full.
        lpoly_add(list, ppx, ppy, aa);
        last(list)->stps = 1;

        if(list->MAXD < 1){
            list->MAXD = 1;
        }

        //printf("ADD \n");
        //printf("POLY X: "); element_out_str(stdout, 10, ppx);
        //printf("\nPOLY Y: "); element_out_str(stdout, 10, ppy); printf("\n");
    }
}

// Modified Miller algorithm:
// It doesn't compute the pairing, but the lines that are needed.

void precompute(lpoly *list, mpz_t q, element_t P, element_t Q) {
  int m;
  element_t Z;
  element_t a, b, c;
  element_t t0;
  const element_ptr cca = curve_a_coeff(P);
  const element_ptr Px = curve_x_coord(P);
  const element_ptr Py = curve_y_coord(P);
  element_ptr Zx, Zy;

  curve_equation = pbc_malloc(sizeof(*curve_equation));
  field_ptr tbase = pbc_malloc(sizeof(*tbase));
  field_init_poly(tbase, curve_a_coeff(P)->field);
  element_t one;
  element_init(one, curve_a_coeff(P)->field);
  element_set1(one);
  element_init(curve_equation, tbase);
  poly_set_coeff(curve_equation, one, 3);
  poly_set_coeff(curve_equation, curve_a_coeff(P), 1);
  poly_set_coeff(curve_equation, curve_b_coeff(P), 0);

  void do_tangent(void) {
    // a = -(3 Zx^2 + cc->a)
    // b = 2 * Zy
    // c = -(2 Zy^2 + a Zx);

    element_square(a, Zx);
    element_mul_si(a, a, 3);
    element_add(a, a, cca);
    element_neg(a, a);

    element_add(b, Zy, Zy);

    element_mul(t0, b, Zy);
    element_mul(c, a, Zx);
    element_add(c, c, t0);
    element_neg(c, c);

    add_poly(list, a, b, c, 1);
  }

  void do_line(void) {
    // a = -(B.y - A.y) / (B.x - A.x);
    // b = 1;
    // c = -(A.y + a * A.x);
    // but we multiply by B.x - A.x to avoid division.

    element_sub(b, Px, Zx);
    element_sub(a, Zy, Py);
    element_mul(t0, b, Zy);
    element_mul(c, a, Zx);
    element_add(c, c, t0);
    element_neg(c, c);

    add_poly(list, a, b, c, 0);
  }

  element_init(a, Px->field);
  element_init(b, a->field);
  element_init(c, a->field);
  element_init(t0, a->field);

  element_init(Z, P->field);

  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_y_coord(Z);

  m = mpz_sizeinbase(q, 2) - 2;

  for(;;) {
    do_tangent();

    if (!m) break;

    element_double(Z, Z);
    if (mpz_tstbit(q, m)) {
      do_line();
      element_add(Z, Z, P);
    }
    m--;
  }

  element_clear(Z);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(t0);
}

// Computes y*X when y in F_qd and X in F_qk
// The output is in F_qk

void point_mult(element_ptr e0, element_ptr x, element_ptr y) {
  element_ptr re_out = element_x(x);
  element_ptr im_out = element_y(x);

  element_mul(re_out, re_out, y); multd++;
  element_mul(im_out, im_out, y); multd++;
}

// Add an element x in F_q to a big element in F_qk
// Used when a polynomial is evaluated

void point_add(element_ptr e0, element_ptr y, element_t x){
   element_ptr re_out = element_x(e0); // e0 = x + iy
   element_ptr y_out = element_x(y); // e0 = x + iy
   //printf("SUMA ");
   // element_out_str(stdout, 10, element_item(y_out, 0));
   // printf(" + ");
   // element_out_str(stdout, 10, x);
   element_add(element_item(re_out, 0), element_item(y_out, 0), x);
   add1++;
   //printf(" = ");
   //element_out_str(stdout, 10, element_item(re_out, 0));
}

// Computes x*Y when x in F_q and y in F_qk (but Im = 0).
// Used when a polynomial in X is evaluated

void basic_mult(element_ptr out, element_t x, element_ptr y){
    element_ptr re_out = element_x(out);

    element_ptr re_y = element_x(y);

    int i;
    int d = polymod_field_degree(re_out->field);

    for (i = 0; i < d; i++) {
      element_mul(element_item(re_out, i), element_item(re_y, i), x); mult1++;
    }
}

// Computes x*Y when x in F_q and y in F_qk (but Re = 0).
// Used when a polynomial in XY is evaluated

void basic_mult2(element_ptr out, element_ptr x, element_ptr y){
    element_ptr im_out = element_y(out);

    element_ptr im_y = element_y(y);

    int i;
    int d = polymod_field_degree(im_out->field);

    for (i = 0; i < d; i++) {
      element_mul(element_item(im_out, i), element_item(im_y, i), x); mult1++;
    }
}

// Computes x + y when x, y in F_qk (but Im = 0)
// Use when a polynomial in X is evaluated

void basic_add(element_ptr out, element_ptr x, element_ptr y){
    element_ptr re_out = element_x(out);
    element_ptr re_x = element_x(x);
    element_ptr re_y = element_x(y);

    element_add(re_out, re_x, re_y); add1 += polymod_field_degree(re_x->field) + 1;
}

// Computes x + y when x, y in F_qk (but Re = 0)
// Use when a polynomial in XY is evaluated

void basic_add2(element_ptr out, element_ptr x, element_ptr y){
    element_ptr im_out = element_y(out);
    element_ptr im_x = element_y(x);
    element_ptr im_y = element_y(y);

    element_add(im_out, im_x, im_y); add1 += polymod_field_degree(im_x->field) + 1;
}

// Computes two polynomials: one in X and the oter in XY.
void compute_polynomial(element_t out, element_t out2, lnodepoly *poly, element_ptr *point_precomp, element_ptr *point_precomp2)
{
    element_ptr p = get_polynomial(poly);
    int i;
    element_t tmp;
    element_init(tmp, out->field);

    element_ptr re_out = element_x(out);
    int d = polymod_field_degree(re_out->field);
    for (i = 0; i < d; i++) {
        element_mul(element_item(re_out, i), element_item(element_x(point_precomp[0]), i), poly_get_coeff(p, 1)); mult1++;
    }

    d = poly_degree(p)+1;

    for(i=2; i < d; i++){
        KSET0(tmp);
        basic_mult(tmp, poly_get_coeff(p, i), point_precomp[i-1]);

        basic_add(out, out, tmp);
    }
    point_add(out, out, poly_get_coeff(p, 0));

    p = get_polynomialY(poly);
    d = poly_degree(p)+1;

    for(i=0; i < d; i++){
        element_set0(tmp);
        basic_mult2(tmp, poly_get_coeff(p, i), point_precomp2[i]);
        basic_add2(out2, out2, tmp);
    }
}

void compute_polynomialN(element_t out, element_t out2, lnodepoly *poly, int row, int MAXD, element_ptr point_precomp[][MAXD], element_ptr point_precomp2[][MAXD])
{
    element_ptr p = get_polynomial(poly);
    int d = poly_degree(p)+1;

    int i;
    element_t tmp;
    element_init(tmp, out->field);

    for(i=1; i < d; i++){
        KSET0(tmp);
        basic_mult(tmp, poly_get_coeff(p, i), point_precomp[row][i-1]);

        basic_add(out, out, tmp);

    }
    point_add(out, out, poly_get_coeff(p, 0));

    p = get_polynomialY(poly);
    d = poly_degree(p)+1;

    element_init(tmp, out->field);
    for(i=0; i < d; i++){
        element_set0(tmp);
        basic_mult2(tmp, poly_get_coeff(p, i), point_precomp2[row][i]);
        basic_add2(out2, out2, tmp);
    }
}

void precomp_miller(element_t res, lpoly *list, element_ptr Qx, element_ptr Qy){
    //struct timeval tvBegin, tvEnd;
    //gettimeofday(&tvBegin, NULL);

    lnodepoly *actual = first(list);

    element_ptr *point_precomp, *point_precomp2;

    point_precomp = (element_ptr *)malloc((list->MAXD)*sizeof(element_ptr));
    point_precomp2 = (element_ptr *)malloc((list->MAXD)*sizeof(element_ptr));

    point_precomp[0] = Qx;
    point_precomp2[0] = Qy;

    int i = 0;
    for(i = 1; i < list->MAXD; i++){
        point_precomp[i] = pbc_malloc(sizeof(element_ptr));
        element_init(point_precomp[i], point_precomp[0]->field);
        element_mul(point_precomp[i], point_precomp[i-1], Qx); multk++;
        //element_mul(tmp, tmp, Qx);

        point_precomp2[i] = pbc_malloc(sizeof(element_ptr));
        element_init(point_precomp2[i], point_precomp[0]->field);
        element_mul(point_precomp2[i], point_precomp2[i-1], Qx); multk++;
        //element_mul(tmp2, tmp2, Qx);
    }

    //gettimeofday(&tvEnd, NULL);
    //printf("Precalc: "); timeval_sub(&tvEnd, &tvBegin);

    element_set1(res);
    element_t out, out2;
    element_init(out, res->field);
    element_init(out2, res->field);

    while(islast(actual) != true){

        int steps = get_steps(actual);
        while(steps > 0 && !element_is1(res)){
            element_square(res,res); squarek++;
            steps--;
        }
        //gettimeofday(&tvEnd, NULL);
        //printf("Step: "); timeval_sub(&tvEnd, &tvBegin);

        //gettimeofday(&tvBegin, NULL);
        KSET0(out);
        KSET0(out2);

        compute_polynomial(out, out2, actual, point_precomp, point_precomp2);

        //gettimeofday(&tvEnd, NULL);
        //printf("Poly: "); timeval_sub(&tvEnd, &tvBegin);

        //gettimeofday(&tvBegin, NULL);

        element_ptr im_out = element_y(out);
        element_set(im_out, element_y(out2));

        element_mul(res, res, out); multk++;

        actual = next(actual);
    }

    for(i=1; i < list->MAXD; i++){
        pbc_free(point_precomp[i]);
        pbc_free(point_precomp2[i]);
    }
    free(point_precomp);
    free(point_precomp2);
    element_clear(out); element_clear(out2);
}

void precomp_millers(element_t res, lpoly list[], element_t Qx[], element_t Qy[], int n_prod) {

    int MAXD = 0;
    int i = 0;
    for(i=0; i < n_prod; i++){
        if(list[i].MAXD > MAXD){
            MAXD = list[i].MAXD;
        }
    }

    element_ptr point_precomp[n_prod][MAXD];
    element_ptr point_precomp2[n_prod][MAXD];

    for(i = 0; i < n_prod; i++){
        point_precomp[i][0] = Qx[i];
        point_precomp2[i][0] = Qy[i];
    }

    int j = 0;
    for(j = 0; j < n_prod; j++){
        for(i = 1; i < MAXD; i++){
            point_precomp[j][i] = pbc_malloc(sizeof(element_ptr));
            element_init(point_precomp[j][i], Qx[0]->field);
            element_mul(point_precomp[j][i], point_precomp[j][i-1], Qx[j]);

            point_precomp2[j][i] = pbc_malloc(sizeof(element_ptr));
            element_init(point_precomp2[j][i], Qx[0]->field);
            element_mul(point_precomp2[j][i], point_precomp2[j][i-1], Qx[j]);
        }
    }

    element_set1(res);
    element_t out, out2;
    element_init(out, res->field);
    element_init(out2, res->field);

    lnodepoly *actuales[n_prod];
    for(i=0; i<n_prod; i++){
        actuales[i] = first(&list[i]);
    }

    int something_more(){
        int j = 0;
        for(j = 0; j < n_prod; j++){
            if(islast(actuales[j]) != true) return 1;
        }
        return 0;
    }

    while(something_more() == 1){
        int steps = get_steps(actuales[0]);
        while(steps > 0){
            element_square(res,res);
            steps--;
        }

        for(i = 0; i < n_prod; i++){
            if(islast(actuales[i]) != true){
                KSET0(out);
                KSET0(out2);

                compute_polynomialN(out, out2, actuales[i], i, MAXD, point_precomp, point_precomp2);

                element_ptr im_out = element_y(out);
                element_set(im_out, element_y(out2));
                element_mul(res, res, out);
            }
        }

        for(i=0; i<n_prod; i++){
            if(islast(actuales[i]) != true) actuales[i] = next(actuales[i]);
        }
    }
}



/******************************
/       PUBLIC FUNCTIONS
*******************************/

void compute_miller(element_ptr out, lpoly *list, element_ptr in2, pairing_t pairing) {
  element_ptr Qbase = in2;
  element_t Qx, Qy;
  pptr p = pairing->data;

  init_contador();

  element_init(Qx, p->Fqd);
  element_init(Qy, p->Fqd);
  // Twist: (x, y) --> (v^-1 x, v^-(3/2) y)
  // where v is the quadratic nonresidue used to construct the twist.
  element_mul(Qx, curve_x_coord(Qbase), p->nqrinv);
  // v^-3/2 = v^-2 * v^1/2
  element_mul(Qy, curve_y_coord(Qbase), p->nqrinv2);

  element_t tmp;
  element_init(tmp, out->field);
  element_set0(tmp);
  KSET0(tmp);
  element_ptr im_out = element_y(tmp); // e0 = x + iy
  element_set(im_out, Qy);


  element_t tmp2;
  element_init(tmp2, out->field);
  //KSET0(tmp2);
  element_set0(tmp2);
  im_out = element_x(tmp2); // e0 = x + iy
  element_set(im_out, Qx);

  precomp_miller(out, list, tmp2, tmp);

  pairing->finalpow(out);

  element_clear(Qx);
  element_clear(Qy);
  element_clear(tmp);
  element_clear(tmp2);
}

void compute_millers(element_ptr out, lpoly *list, element_ptr in2[], int n_prod, pairing_t pairing) {
  element_t Qx[n_prod], Qy[n_prod];
  element_t tmp[n_prod], tmp2[n_prod];
  pptr p = pairing->data;

  init_contador();

  int i = 0;
  for(i = 0; i < n_prod; i++){
    element_init(Qx[i], p->Fqd);
    element_init(Qy[i], p->Fqd);

    element_mul(Qx[i], curve_x_coord(in2[i]), p->nqrinv);
    element_mul(Qy[i], curve_y_coord(in2[i]), p->nqrinv2);

    element_init(tmp[i], out->field);
    element_init(tmp2[i], out->field);
    element_set0(tmp[i]); KSET0(tmp2[i]);

    element_ptr im_out = element_y(tmp[i]); element_set0(im_out);
    element_ptr re_out = element_x(tmp[i]); element_set0(re_out);
    element_add(im_out, im_out, Qy[i]);
    im_out = element_x(tmp2[i]);
    element_add(im_out, im_out, Qx[i]);
  }

  precomp_millers(out, list, tmp2, tmp, n_prod);

  pairing->finalpow(out);
}

void eval_miller(element_ptr out, element_ptr in1, element_ptr in2,
    pairing_t pairing) {
  init_contador();
  element_ptr Qbase = in2;
  element_t Qx, Qy;
  pptr p = pairing->data;

  element_init(Qx, p->Fqd);
  element_init(Qy, p->Fqd);
  // Twist: (x, y) --> (v^-1 x, v^-(3/2) y)
  // where v is the quadratic nonresidue used to construct the twist.
  element_mul(Qx, curve_x_coord(Qbase), p->nqrinv);
  // v^-3/2 = v^-2 * v^1/2
  element_mul(Qy, curve_y_coord(Qbase), p->nqrinv2);

  miller(out, pairing->r, in1, Qx, Qy);

  pairing->finalpow(out);

  element_clear(Qx);
  element_clear(Qy);
}

void count_elem(lpoly *list){
    lnodepoly *actual = first(list);

    int total;
    total = 0;

    while(islast(actual) != true){
        total = total + poly_degree(get_polynomial(actual)) + poly_degree(get_polynomialY(actual)) + 2;

        actual = next(actual);
    }

    printf("%i\n", total);
}
