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
int STEPMAX = 1;

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

/*int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}*/

void KSET0(element_t out){
    element_set0(out);
    element_ptr re_out = element_x(out);
    element_set0(element_item(re_out,0));
}





/*****************************
/  DATA STRUCTURE OF F-curves
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


struct f_pairing_data_s {
  field_t Fq, Fq2, Fq2x, Fq12;
  field_t Eq, Etwist;
  element_t negalpha;
  element_t negalphainv;
  mpz_t tateexp;

  //for tate exponentiation speedup:
  //x^{q^k} for various k
  element_t xpowq2, xpowq6, xpowq8;
};
typedef struct f_pairing_data_s f_pairing_data_t[1];
typedef struct f_pairing_data_s *f_pairing_data_ptr;


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

void miller(element_t res, mpz_t q, element_t P, element_ptr Qx, element_ptr Qy, element_t negalpha) {
  int m;
  element_t v;
  element_t Z;
  element_t a, b, c;
  element_t t0;
  element_t e0, e1;
  element_ptr Zx, Zy;
  const element_ptr Px = curve_x_coord(P);
  const element_ptr Py = curve_y_coord(P);

  void f_miller_evalfn(void)
  {
    // a, b, c lie in Fq
    // Qx, Qy lie in Fq^2
    // Qx is coefficient of x^4
    // Qy is coefficient of x^3
    //
    // computes v *= (a Qx x^4 + b Qy x^3 + c)
    //
    // recall x^6 = -alpha thus
    // x^4 (u0 + u1 x^1 + ... + u5 x^5) =
    // u0 x^4 + u1 x^5
    // - alpha u2 - alpha u3 x - alpha u4 x^2 - alpha u5 x^3
    // and
    // x^4 (u0 + u1 x^1 + ... + u5 x^5) =
    // u0 x^3 + u1 x^4 + u2 x^5
    // - alpha u3 - alpha u4 x - alpha u5 x^2
    element_ptr e2;
    inline void do_term(int i, int j, int k, int flag)
    {
      e2 = element_item(e0, i);
      element_mul(e1, element_item(v, j), Qx); mult1++;
      if (flag == 1) {element_mul(e1, e1, negalpha); mult1++; }
      element_mul(element_x(e1), element_x(e1), a); mult1++;
      element_mul(element_y(e1), element_y(e1), a); mult1++;
      element_mul(e2, element_item(v, k), Qy); mult1++;
      element_mul(element_x(e2), element_x(e2), b); mult1++;
      element_mul(element_y(e2), element_y(e2), b); mult1++;
      element_add(e2, e2, e1); add1++;
      if (flag == 2) {element_mul(e2, e2, negalpha); mult1++; }
      element_mul(element_x(e1), element_x(element_item(v, i)), c); mult1++;
      element_mul(element_y(e1), element_y(element_item(v, i)), c); mult1++;
      element_add(e2, e2, e1); add1++;
    }

    do_term(0, 2, 3, 2);
    do_term(1, 3, 4, 2);
    do_term(2, 4, 5, 2);
    do_term(3, 5, 0, 1);
    do_term(4, 0, 1, 0);
    do_term(5, 1, 2, 0);

    element_set(v, e0);

    /*
    element_ptr e1;

    e1 = element_item(e0, 4);

    element_mul(element_x(e1), element_x(Qx), a);
    element_mul(element_y(e1), element_y(Qx), a);

    e1 = element_item(e0, 3);

    element_mul(element_x(e1), element_x(Qy), b);
    element_mul(element_y(e1), element_y(Qy), b);

    element_set(element_x(element_item(e0, 0)), c);

    element_mul(v, v, e0);
    */
  }

  void do_tangent(void)
  {
    //a = -3 Zx^2 since cc->a is 0 for D = 3
    //b = 2 * Zy
    //c = -(2 Zy^2 + a Zx);
    element_square(a, Zx); mult1++;
    element_mul_si(a, a, 3); add1+=3;
    element_neg(a, a);

    element_add(b, Zy, Zy); add1++;

    element_mul(t0, b, Zy); mult1++;
    element_mul(c, a, Zx); mult1++;
    element_add(c, c, t0); add1++;
    element_neg(c, c);

    f_miller_evalfn();
  }

  void do_line(void)
  {
    //a = -(B.y - A.y) / (B.x - A.x);
    //b = 1;
    //c = -(A.y + a * A.x);
    //but we'll multiply by B.x - A.x to avoid division

    element_sub(b, Px, Zx); add1++;
    element_sub(a, Zy, Py); add1++;
    element_mul(t0, b, Zy); mult1++;
    element_mul(c, a, Zx); mult1++;
    element_add(c, c, t0); add1++;
    element_neg(c, c);

    f_miller_evalfn();
  }

  element_init(a, Px->field);
  element_init(b, a->field);
  element_init(c, a->field);
  element_init(t0, a->field);
  element_init(e0, res->field);
  element_init(e1, Qx->field);

  element_init(v, res->field);
  element_init(Z, P->field);

  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_y_coord(Z);

  element_set1(v);
  m = mpz_sizeinbase(q, 2) - 2;

  //TODO: sliding NAF
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
  element_clear(e1);
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

    //if(last(list)->n != -1 && aa == 1 && (2*STEPMAX - last(list)->stps -1) >= aa){
    if(last(list)->n != -1 && aa == 1 && (STEPMAX > last(list)->stps)){
        // If there is some step in the list
        // and there is room for a new polynomial
        //    ( we count a Doubling step as 2, and an Add step as 1 )
        //    ( because Doubling has an square, and Add not )

        // INIT

        element_ptr polyx = last(list)->polynomial;
        element_ptr polyy = last(list)->polynomialY;
        /*printf("\n MODIFICANDO \n");
        printf("aa: %i, n: %i \n", aa, last(list)->n);
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

        //printf("\n ADD \n");
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

void point_mult(element_t e0, element_t x, element_t y) {
  element_ptr re_x = element_x(x);
  element_ptr im_x = element_y(x);
  element_ptr re_out = element_x(e0);
  element_ptr im_out = element_y(e0);

  element_mul(re_out, re_x, y); mult1++;
  element_mul(im_out, im_x, y); mult1++;
}

// Add an element x in F_q to a big element in F_qk
// Used when a polynomial is evaluated

void point_add(element_t e0, element_t y, element_t x){
   element_ptr re_out = element_x(e0); // e0 = x + iy
   element_ptr y_out = element_x(y); // e0 = x + iy
   //printf("SUMA ");
   // element_out_str(stdout, 10, element_item(y_out, 0));
   // printf(" + ");
   // element_out_str(stdout, 10, x);
   element_add(element_item(re_out, 0), element_item(y_out, 0), x); add1+=12;
   //printf(" = ");
   //element_out_str(stdout, 10, element_item(re_out, 0));
}

// Computes x*Y when x in F_q and y in F_qk (but Im = 0).
// Used when a polynomial in X is evaluated

void basic_mult(element_t out, element_t x, element_ptr y){
    if(!element_is0(element_item(y,0))) point_mult(element_item(out, 0), element_item(y, 0), x);
    if(!element_is0(element_item(y,2))) point_mult(element_item(out, 2), element_item(y, 2), x);
    if(!element_is0(element_item(y,4))) point_mult(element_item(out, 4), element_item(y, 4), x);
}

// Computes x*Y when x in F_q and y in F_qk (but Re = 0).
// Used when a polynomial in XY is evaluated

void basic_mult2(element_t out, element_t x, element_ptr y){
    if(!element_is0(element_item(y,1))) point_mult(element_item(out, 1), element_item(y, 1), x);
    if(!element_is0(element_item(y,3))) point_mult(element_item(out, 3), element_item(y, 3), x);
    if(!element_is0(element_item(y,5))) point_mult(element_item(out, 5), element_item(y, 5), x);
}

// Computes x + y when x, y in F_qk (but Im = 0)
// Use when a polynomial in X is evaluated

void basic_add(element_t out, element_t x, element_t y){
    element_add(element_item(out,0), element_item(x,0), element_item(y,0));
    element_add(element_item(out,1), element_item(x,1), element_item(y,1));
    element_add(element_item(out,2), element_item(x,2), element_item(y,2));
    element_add(element_item(out,3), element_item(x,3), element_item(y,3));
    element_add(element_item(out,4), element_item(x,4), element_item(y,4));
    element_add(element_item(out,5), element_item(x,5), element_item(y,5));
    add1+=6;
}

// Computes x + y when x, y in F_qk (but Re = 0)
// Use when a polynomial in XY is evaluated

void basic_add2(element_t out, element_t x, element_t y){
    basic_add(out, x, y);
}

// Computes two polynomials: one in X and the oter in XY.
void compute_polynomial(element_t out, element_t out2, lnodepoly *poly, element_ptr *point_precomp, element_ptr *point_precomp2)
{
    element_ptr p = get_polynomial(poly);
    int d = poly_degree(p)+1;

    int i;
    element_t tmp;
    element_init(tmp, out->field);

    basic_mult(out, poly_get_coeff(p,1), point_precomp[0]);

    for(i=2; i < d; i++){
        KSET0(tmp);
        //struct timeval tvBegin, tvEnd, tvDiff;
        //gettimeofday(&tvBegin, NULL);

        basic_mult(tmp, poly_get_coeff(p, i), point_precomp[i-1]);

        //gettimeofday(&tvEnd, NULL);
        //timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
        //printf("MULT %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);

        //gettimeofday(&tvBegin, NULL);

        basic_add(out, out, tmp);

        //gettimeofday(&tvEnd, NULL);
        //timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
        //printf("ADD %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);
    }
    point_add(out, out, poly_get_coeff(p, 0));

    p = get_polynomialY(poly);
    d = poly_degree(p)+1;

    element_init(tmp, out->field);

    basic_mult2(out2, poly_get_coeff(p,0), point_precomp2[0]);

    for(i=1; i < d; i++){
        KSET0(tmp);
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

    basic_mult(out, poly_get_coeff(p,1), point_precomp[row][0]);

    for(i=2; i < d; i++){
        KSET0(tmp);
        basic_mult(tmp, poly_get_coeff(p, i), point_precomp[row][i-1]);
        //struct timeval tvBegin, tvEnd, tvDiff;
        //gettimeofday(&tvBegin, NULL);

        basic_add(out, out, tmp);

        //gettimeofday(&tvEnd, NULL);
        //timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
        //printf("ADD %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);
    }
    point_add(out, out, poly_get_coeff(p, 0));

    p = get_polynomialY(poly);
    d = poly_degree(p)+1;

    basic_mult(out2, poly_get_coeff(p,0), point_precomp2[row][0]);

    element_init(tmp, out->field);
    for(i=1; i < d; i++){
        KSET0(tmp);
        basic_mult2(tmp, poly_get_coeff(p, i), point_precomp2[row][i]);
        basic_add2(out2, out2, tmp);
    }
}

void advc_mult(element_ptr out, element_t x, element_t y, element_t negalpha){
    void mult3(element_ptr tmp, element_t n0, element_t n1, element_t n2, element_t a0, element_t a1, element_t a2){
        element_mul(element_item(tmp,4), n2, a2);

        element_t d0, d1, tmp2;
        element_init(d0, n0->field); element_mul(d0, n0, a0);
        element_init(d1, n0->field); element_mul(d1, n1, a1);
        element_init(tmp2, n0->field);

        element_add(tmp2, n1, n2);
        element_add(element_item(tmp, 3), a1, a2);
        element_mul(element_item(tmp,3), element_item(tmp,3), tmp2);
        element_sub(element_item(tmp,3), element_item(tmp,3), d0);
        element_sub(element_item(tmp,3), element_item(tmp,3), d1);

        element_add(tmp2, n0, n2);
        element_add(element_item(tmp,2), a0, a2);
        element_mul(element_item(tmp,2), element_item(tmp,2), tmp2);
        element_sub(element_item(tmp,2), element_item(tmp,2), element_item(tmp,4));
        element_sub(element_item(tmp,2), element_item(tmp,2), d0);
        element_add(element_item(tmp,2), element_item(tmp,2), d1);

        element_add(tmp2, n0, n1);
        element_add(element_item(tmp, 1), a0, a1);
        element_mul(element_item(tmp,1), element_item(tmp,1), tmp2);
        element_sub(element_item(tmp,1), element_item(tmp,1), d1);
        element_sub(element_item(tmp,1), element_item(tmp,1), d0);

        element_add(element_item(tmp,0), element_item(tmp,0), d0);
    }

    element_ptr D0, D1, D2;
    D0 = malloc(sizeof(element_ptr)); D1 = malloc(sizeof(element_ptr)); D2 = malloc(sizeof(element_ptr));
    element_init(D0, out->field); element_init(D1, out->field); element_init(D2, out->field);
    element_ptr a[6];
    element_ptr n[6];
    int i = 0;
    for(i = 0; i < 6; i++){
        a[i] = element_item(x, i);
        n[i] = element_item(y, i);
    }
    element_t sa[3];
    element_t sn[3];
    for(i = 0; i<3; i++){
        element_init(sa[i], a[0]->field);
        element_add(sa[i], a[i], a[i+3]);

        element_init(sn[i], n[0]->field);
        element_add(sn[i], n[i], n[i+3]);
    }
    mult3(D0, a[0], a[1], a[2], n[0], n[1], n[2]);
    mult3(D1, a[3], a[4], a[5], n[3], n[4], n[5]);
    mult3(D2, sa[0], sa[1], sa[2], sn[0], sn[1], sn[2]);

    // D1x^6 + (D2-D0-D1)X^3 + D0
    for(i=0; i < 6; i++){
        element_set(element_item(out, i), element_item(D0, i));
    }

    // Coef 0
    element_set(a[0], element_item(D2, 3));
    element_sub(a[0], a[0], element_item(D0, 3));
    element_sub(a[0], a[0], element_item(D1, 3));
    element_sub(a[0], a[0], element_item(D1, 0));
    element_mul(a[0], a[0], negalpha);
    element_add(element_item(out,0), element_item(out,0), a[0]);

    // Coef 1
    element_set(a[0], element_item(D2, 4));
    element_sub(a[0], a[0], element_item(D0, 4));
    element_sub(a[0], a[0], element_item(D1, 4));
    element_sub(a[0], a[0], element_item(D1, 1));
    element_mul(a[0], a[0], negalpha);
    element_add(element_item(out,1), element_item(out,1), a[0]);

    // Coef 2
    element_set(a[0], element_item(D1, 2));
    element_mul(a[0], a[0], negalpha);
    element_add(element_item(out,2), element_item(out,2), a[0]);

    // Coef 3
    element_set(a[0], element_item(D1, 3));
    element_mul(a[0], a[0], negalpha);
    element_add(a[0], a[0], element_item(D2, 0));
    element_sub(a[0], a[0], element_item(D1, 0));
    element_sub(a[0], a[0], element_item(D0, 0));
    element_add(element_item(out,3), element_item(out,3), a[0]);

    // Coef 4
    element_set(a[0], element_item(D1, 4));
    element_mul(a[0], a[0], negalpha);
    element_add(a[0], a[0], element_item(D2, 1));
    element_sub(a[0], a[0], element_item(D1, 1));
    element_sub(a[0], a[0], element_item(D0, 1));
    element_add(element_item(out,4), element_item(out,4), a[0]);

    // Coef 5
    element_set(a[0], element_item(D2, 2));
    element_sub(a[0], a[0], element_item(D1, 2));
    element_sub(a[0], a[0], element_item(D0, 2));
    element_add(element_item(out,3), element_item(out,3), a[0]);
}

void good_mult(element_ptr out, element_t x, element_t y, element_t negalpha){
    //advc_mult(out, x, y, negalpha);
    //return;
    element_t tmp, tmp2, tmp3;
    element_init(tmp, element_item(x,0)->field);
    element_init(tmp2, element_item(x,0)->field);
    element_init(tmp3, x->field); KSET0(tmp3);

    int i = 0; int k;
    for(i = 0; i < 6; i++){
        element_set0(tmp); element_set0(element_item(tmp, 0));
        int j;
        k = (i+1)%6;
        for(j = 5; j >= 0; j--){
            if(!element_is0(element_item(x,j)) && !element_is0(element_item(y,k))){
                element_set0(tmp2); element_set0(element_item(tmp2, 0));
                element_mul(tmp2, element_item(x, j), element_item(y, k)); mult1++;

                element_add(tmp, tmp, tmp2); add1++;
            }

            if(j == i+1 && !element_is0(tmp)){
                element_mul(tmp, tmp, negalpha); mult1++;
            }

            k++;
            if(k>5) k = k - 6;
        }
        element_add(element_item(tmp3, i), element_item(tmp3,i), tmp); add1++;
    }

    element_set(out, tmp3);
}

void precomp_miller(element_t res, lpoly *list, element_ptr Qx, element_ptr Qy, element_t negalpha){
    lnodepoly *actual = first(list);

    //struct timeval tvBegin, tvEnd, tvDiff;

    element_ptr *point_precomp, *point_precomp2;

    point_precomp = (element_ptr *)malloc((list->MAXD)*sizeof(element_ptr));
    point_precomp2 = (element_ptr *)malloc((list->MAXD)*sizeof(element_ptr));

    point_precomp[0] = Qx;
    point_precomp2[0] = Qy;

    int i = 0;
    for(i = 1; i < list->MAXD; i++){
        point_precomp[i] = pbc_malloc(sizeof(element_ptr));
        element_init(point_precomp[i], Qx->field);
        element_mul(point_precomp[i], point_precomp[i-1], Qx); multk++;
        //element_mul(tmp, tmp, Qx);

        point_precomp2[i] = pbc_malloc(sizeof(element_ptr));
        element_init(point_precomp2[i], Qx->field);
        element_mul(point_precomp2[i], point_precomp2[i-1], Qx); multk++;
        //element_mul(tmp2, tmp2, Qx);
    }

    element_set1(res);
    element_t out, out2;
    element_init(out, res->field);
    element_init(out2, res->field);

    while(islast(actual) != true){
        //gettimeofday(&tvBegin, NULL);

        int steps = get_steps(actual);
        while(steps > 0){
            element_square(res,res); squarek++;
            steps--;
        }
        //gettimeofday(&tvEnd, NULL);
        //timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
        //printf("STEPS %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);

        KSET0(out);
        KSET0(out2);

        //gettimeofday(&tvBegin, NULL);

        compute_polynomial(out, out2, actual, point_precomp, point_precomp2);

        //gettimeofday(&tvEnd, NULL);
        //timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
        //printf("EVAL %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);

        basic_add(out, out, out2);

        //gettimeofday(&tvBegin, NULL);

        // "Good" multiplication
        /*element_t tmp, tmp2, tmp3;
        element_init(tmp, element_item(Qx,0)->field);
        element_init(tmp2, element_item(Qx,0)->field);
        element_init(tmp3, Qx->field); KSET0(tmp3);
        for(i = 0; i < 6; i++){
           element_set0(tmp); element_set0(element_item(tmp, 0));
           int j;
           for(j = 5; j >= 0; j--){
               element_set0(tmp2); element_set0(element_item(tmp2, 0));

               if(i < j) { element_mul(tmp2, element_item(res, j), element_item(out, i-j+6)); }
               if(i >= j) { element_mul(tmp2, element_item(res, j), element_item(out, i-j)); }

               element_add(tmp, tmp, tmp2);

               if(j == i+1){
                   element_mul(tmp, tmp, negalpha);
               }
           }
           element_add(element_item(tmp3, i), element_item(tmp3,i), tmp);
        }

        element_set(res, tmp3);*/
        good_mult(res, res, out, negalpha);

        //gettimeofday(&tvEnd, NULL);
        //timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
        //printf("MULT %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);

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

void precomp_millers(element_t res, lpoly list[], element_t Qx[], element_t Qy[], int n_prod, element_t negalpha) {

    int MAXD = 0;
    int i = 0;
    for(i=0; i < n_prod; i++){
        if(list[i].MAXD > MAXD){
            MAXD = list[i].MAXD;
        }
    }

    element_ptr point_precomp[n_prod][MAXD];
    element_ptr point_precomp2[n_prod][MAXD];
    // point_precomp[ i*list.MAXD + j ] <- Qx[j]**i
    //point_precomp = (element_ptr *)malloc((n_prod*list[0].MAXD)*sizeof(element_ptr));
    //point_precomp2 = (element_ptr *)malloc((n_prod*list[0].MAXD)*sizeof(element_ptr));

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

    while(islast(actuales[0]) != true){
        int steps = get_steps(actuales[0]);
        while(steps > 0){
            element_square(res,res);
            steps--;
        }

        for(i = 0; i < n_prod; i++){
            KSET0(out);
            KSET0(out2);

            compute_polynomialN(out, out2, actuales[i], i, MAXD, point_precomp, point_precomp2);

            basic_add(out, out, out2);

            //element_mul(res, res, out);
            good_mult(res, res, out, negalpha);
        }

        for(i=0; i<n_prod; i++){
            actuales[i] = next(actuales[i]);
        }
    }
}



/******************************
/       PUBLIC FUNCTIONS
*******************************/

void compute_miller(element_ptr out, lpoly *list, element_ptr in2, pairing_t pairing) {
  init_contador();

  element_ptr Qbase = in2;
  element_t Qx, Qy;
  f_pairing_data_ptr p = pairing->data;

  element_init(Qx, p->Fq2);
  element_init(Qy, p->Fq2);
  // Twist: (x, y) --> (v^-1 x, v^-(3/2) y)
  // where v is the quadratic nonresidue used to construct the twist.
  element_mul(Qx, curve_x_coord(Qbase), p->negalphainv);
  // v^-3/2 = v^-2 * v^1/2
  element_mul(Qy, curve_y_coord(Qbase), p->negalphainv);

  element_t tmp;
  element_init(tmp, out->field);
  KSET0(tmp);
  element_ptr im_out = element_item(tmp, 3);
  element_set(im_out, Qy);

  element_t tmp2;
  element_init(tmp2, out->field);
  KSET0(tmp2);
  im_out = element_item(tmp2, 4); // e0 = x + iy
  element_set(im_out, Qx);

  precomp_miller(out, list, tmp2, tmp, p->negalpha);

  pairing->finalpow(out);

  element_clear(Qx);
  element_clear(Qy);
  element_clear(tmp);
  element_clear(tmp2);
}

void compute_millers(element_ptr out, lpoly *list, element_ptr in2[], int n_prod, pairing_t pairing) {
  init_contador();

  element_t Qx[n_prod], Qy[n_prod];
  element_t tmp[n_prod], tmp2[n_prod];
  f_pairing_data_ptr p = pairing->data;

  int i = 0;
  for(i = 0; i < n_prod; i++){
    element_init(Qx[i], p->Fq2);
    element_init(Qy[i], p->Fq2);


    element_mul(Qx[i], curve_x_coord(in2[i]), p->negalphainv);
    element_mul(Qy[i], curve_y_coord(in2[i]), p->negalphainv);

    element_init(tmp[i], out->field);
    element_init(tmp2[i], out->field);
    KSET0(tmp[i]); KSET0(tmp2[i]);

    element_ptr im_out = element_item(tmp[i], 3); element_set0(im_out);
    element_add(im_out, im_out, Qy[i]);
    im_out = element_item(tmp2[i], 4);
    element_add(im_out, im_out, Qx[i]);
  }

  precomp_millers(out, list, tmp2, tmp, n_prod, p->negalpha);

  pairing->finalpow(out);
}

void eval_miller(element_ptr out, element_ptr in1, element_ptr in2,
    pairing_t pairing) {
  init_contador();

  element_ptr Qbase = in2;
  element_t Qx, Qy;
  f_pairing_data_ptr p = pairing->data;

  element_init(Qx, p->Fq2);
  element_init(Qy, p->Fq2);
  // Twist: (x, y) --> (v^-1 x, v^-(3/2) y)
  // where v is the quadratic nonresidue used to construct the twist.
  element_mul(Qx, curve_x_coord(Qbase), p->negalphainv);
  // v^-3/2 = v^-2 * v^1/2
  element_mul(Qy, curve_y_coord(Qbase), p->negalphainv);

  miller(out, pairing->r, in1, Qx, Qy, p->negalpha);

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
