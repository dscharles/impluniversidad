#include <pbc.h>
#include <pbc_poly.h>
#include <pbc_curve.h>
#include <pbc_f_param.h>
#include <stdio.h>
#include "lpoly.h"
#include "precompute.h"
#include <time.h>
#include <sys/time.h>

struct f_param_s {
    mpz_t q; // Curve defined over F_q.
    mpz_t r; // The order of the curve.
    mpz_t b; // E: y^2 = x^3 + b
    mpz_t beta; //beta is a quadratic nonresidue in Fq
        //we use F_q^2 = F_q[sqrt(beta)]
    mpz_t alpha0, alpha1;
        //the polynomial x^6 + alpha0 + alpha1 sqrt(beta)
        //is irreducible over F_q^2[x], so
        //we can extend F_q^2 to F_q^12 using the
        //sixth root of -(alpha0 + alpha1 sqrt(beta))
};
typedef struct f_param_s f_param_t[1];
typedef struct f_param_s *f_param_ptr;

void timeval_subtract(struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);

    printf("%lu\n", diff);
}

int main(void)
{

	pairing_t pairing;
    char param[5000];
    size_t count = fread(param, 1, 5000, stdin);
    if (!count) pbc_die("input error");
    pairing_init_set_buf(pairing, param, count);

    struct timeval tvBegin, tvEnd;

    element_t g, h;
    element_t public_key, secret_key;
    element_t sig;
    element_t temp1, temp2;

    element_init_G2(g, pairing);
    element_init_G2(public_key, pairing);
    element_init_G1(h, pairing);
    element_init_G1(sig, pairing);
    element_init_GT(temp1, pairing);
    element_init_GT(temp2, pairing);
    element_init_Zr(secret_key, pairing);

    // Generating key
    element_random(g);
    element_random(secret_key);
    element_pow_zn(public_key, g, secret_key);

    // Generating message
    element_from_hash(h, "ABCDEF", 6);

    element_pow_zn(sig, h, secret_key);


    int contador = 0;
    // RANDOM TESTS

    // Fp

    element_t p1, p2;
    element_init(p1, element_x(h)->field);
    element_init(p2, p1->field);
    element_random(p1);
    element_random(p2);

    // multiplication

    gettimeofday(&tvBegin, NULL);

    for(contador = 0; contador < 1000; contador++) element_mul(p1, p1, p2);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);

    //square

    gettimeofday(&tvBegin, NULL);

    for(contador = 0; contador < 1000; contador++) element_square(p1, p2);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);

    // add

    gettimeofday(&tvBegin, NULL);

    for(contador = 0; contador < 1000; contador++) element_add(p1, p1, p2);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);

    // invers

    gettimeofday(&tvBegin, NULL);

    for(contador = 0; contador < 1000; contador++) element_invert(p1, p1);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);






    // Fpk

    element_t q1, q2;
    element_init_GT(q1, pairing);
    element_init_GT(q2, pairing);
    element_random(q1);
    element_random(q2);

    // multiplication

    gettimeofday(&tvBegin, NULL);

    for(contador = 0; contador < 1000; contador++) element_mul(q1, q1, q2);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);

    //square

    gettimeofday(&tvBegin, NULL);

    for(contador = 0; contador < 1000; contador++) element_square(q1, q2);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);

    // add

    gettimeofday(&tvBegin, NULL);

    for(contador = 0; contador < 1000; contador++) element_add(q1, q1, q2);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);

    // invers

    gettimeofday(&tvBegin, NULL);

    for(contador = 0; contador < 1000; contador++) element_invert(q1, q1);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);






    // CURVE OPERATIONS

    element_t punto, punto2;
    element_init(punto, h->field); element_random(punto);
    element_init(punto2, h->field); element_random(punto2);

    // add

    gettimeofday(&tvBegin, NULL);

    element_mul(punto, punto, punto2);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);

    // double

    gettimeofday(&tvBegin, NULL);

    element_double(punto, punto2);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);





   // SIZE GROUP
    int m = mpz_sizeinbase(pairing->r, 2) - 2;
    printf("%i\n",  m);
    contador = 0;
    for(;;){
        if(!m) break;
        if(mpz_tstbit(pairing->r,m)) contador++;
        m--;
    }
    printf("%i\n", contador);






    // One pairing
    gettimeofday(&tvBegin, NULL);

    eval_miller(temp1, sig, g, pairing);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);







    // One pairing (with precomputed values)

    // Original method


    pairing_pp_t pp;
    // Precomp
    gettimeofday(&tvBegin, NULL);

    pairing_pp_init(pp, sig, pairing);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);

    // Eval
    gettimeofday(&tvBegin, NULL);

    pairing_pp_apply(temp1, g, pp);

    gettimeofday(&tvEnd, NULL);
    timeval_subtract(&tvEnd, &tvBegin);


    pairing_pp_clear(pp);

    void do_precomp(){
        lpoly *list;

        // precomputation
        gettimeofday(&tvBegin, NULL);

        list = lpoly_init();
        precompute(list, pairing->r, sig, g);

        gettimeofday(&tvEnd, NULL);
        timeval_subtract(&tvEnd, &tvBegin);

        // DMAX
        printf("%i\n", list->MAXD);

        // eval
        gettimeofday(&tvBegin, NULL);

        compute_miller(temp2, list, g, pairing);

        gettimeofday(&tvEnd, NULL);
        timeval_subtract(&tvEnd, &tvBegin);

        lpoly_free(list);
    }

    // n = 1
    change_NMAX(1);
    do_precomp();

    // n = 2
    change_NMAX(2);
    do_precomp();

    // n = 3
    change_NMAX(3);
    do_precomp();

    // n = 4
    change_NMAX(4);
    do_precomp();

    // Multipairing

    void do_multi(int m){
        int i = 0;
        lpoly list[m];
        lpoly *tmp_list;

        // prevalues

        element_t gg[m];
        element_ptr ggg[m];
        element_t hh[m];
        for(i = 0; i < m; i++){
            element_init_G2(gg[i], pairing);
            element_random(gg[i]);
            ggg[i] = malloc(sizeof(element_ptr));
            element_init(ggg[i], gg[i]->field);
            element_set(ggg[i], gg[i]);

            element_init_G1(hh[i], pairing);
            element_random(hh[i]);
        }


        // precomputation
        gettimeofday(&tvBegin, NULL);

        for(i = 0; i < m; i++){
            tmp_list = lpoly_init();
            precompute(tmp_list, pairing->r, hh[i], gg[i]);
            list[i] = *tmp_list;
        }

        gettimeofday(&tvEnd, NULL);
        timeval_subtract(&tvEnd, &tvBegin);

        // compute
        gettimeofday(&tvBegin, NULL);
        compute_millers(temp2, list, ggg, m, pairing);
        gettimeofday(&tvEnd, NULL);
        timeval_subtract(&tvEnd, &tvBegin);

        gettimeofday(&tvBegin, NULL);
        element_prod_pairing(temp1, hh, gg, m);
        gettimeofday(&tvEnd, NULL);
        timeval_subtract(&tvEnd, &tvBegin);
    }

    /*change_NMAX(1);
    do_multi(10); do_multi(30);

    change_NMAX(2);
    do_multi(10); do_multi(30);

    change_NMAX(3);
    do_multi(10); do_multi(30);

    return 0;*/

    // n = 1
    change_NMAX(1);
    do_multi(2); do_multi(3); do_multi(4); do_multi(5);

    // n = 2
    change_NMAX(2);
    do_multi(2); do_multi(3); do_multi(4); do_multi(5);

    // n = 3
    change_NMAX(3);
    do_multi(2); do_multi(3); do_multi(4); do_multi(5);

    // n = 4
    change_NMAX(4);
    do_multi(2); do_multi(3); do_multi(4); do_multi(5);

	return 0;






	void do_count(){
        lpoly *list;

        // precomputation

        list = lpoly_init();
        precompute(list, pairing->r, sig, g);

        // DMAX
        printf("%i\n", list->MAXD);

        count_elem(list);
    }

    change_NMAX(1);
    do_count();

    change_NMAX(2);
    do_count();

    change_NMAX(3);
    do_count();

    change_NMAX(4);
    do_count();

	return 0;
}
