#include <stdio.h>
#include <stdlib.h>
#include "sylvan_edge_weights_qisq2.h"


/**********************<Some static utility functions>*************************/

static void
comp_cart_to_polar(fl_t r, fl_t i, fl_t *magnitude, fl_t *angle)
{
    // TODO
}

/**********************</Some static utility functions>************************/




/*****************<Implementation of edge_weights interface>*******************/

qisq2_t *
weight_qisq2_malloc()
{
    qisq2_t *res = malloc(sizeof(qisq2_t));
    return res;
}

void
_weight_qisq2_value(void *wgt_store, EVBDD_WGT a, qisq2_t *res)
{
    // TODO
}

EVBDD_WGT
_weight_qisq2_lookup_ptr(qisq2_t *a, void *wgt_store)
{
    // TODO
}

EVBDD_WGT
weight_qisq2_lookup(qisq2_t *a)
{
    return _weight_qisq2_lookup_ptr(a, wgt_storage);
}

void
init_qisq2_one_zero(void *wgt_store)
{
    qisq2_t a;
    a = qisq2_one();     EVBDD_ONE     = _weight_qisq2_lookup_ptr(&a, wgt_store);
    a = qisq2_zero();    EVBDD_ZERO    = _weight_qisq2_lookup_ptr(&a, wgt_store);
    a = qisq2_mone();    EVBDD_MIN_ONE = _weight_qisq2_lookup_ptr(&a, wgt_store);

}

void
weight_qisq2_abs(qisq2_t *a)
{
    qisq2_t *temp;
    init(temp);
    complexConjugate(temp,a);
    weight_qisq2_mul(a,temp);
    clear(temp);
}

void
weight_qisq2_neg(qisq2_t *num)
{
    mpq_neg(num->a, num->a);
    mpq_neg(num->b, num->b);
    mpq_neg(num->c, num->c);
    mpq_neg(num->d, num->d);
    reduce(num);
}

void
weight_qisq2_conj(qisq2_t *x)
{
    mpq_neg(x->c, x->c);
    mpq_neg(x->d, x->d);
    reduce(x);
}

void
weight_qisq2_sqr(qisq2_t *a)
{
    weight_qisq2_mul(a, a);
}

void
weight_qisq2_add(qisq2_t *x, qisq2_t *y)
{
    mpq_add(x->a, x->a, y->a);
    mpq_add(x->b, x->b, y->b);
    mpq_add(x->c, x->c, y->c);
    mpq_add(x->d, x->d, y->d);
    reduce(x);
}

void
weight_qisq2_sub(qisq2_t *x, qisq2_t *y)
{
    mpq_sub(x->a, x->a, y->a);
    mpq_sub(x->b, x->b, y->b);
    mpq_sub(x->c, x->c, y->c);
    mpq_sub(x->d, x->d, y->d);
    reduce(x);
}

void
weight_qisq2_mul(qisq2_t *x, qisq2_t *y)
{
    qisq2_t *result;
    init(result);
    mpq_t temp;
    mpq_init(temp);
    mpq_t two;
    mpq_init(two);
    mpq_set_ui(two, 2, 1);

    mpq_mul(result->a, x->a, y->a);
    mpq_mul(temp, x->b, y->b);
    mpq_mul(temp, temp, two);
    mpq_add(result->a, result->a, temp);

    mpq_mul(temp, x->c, y->c);
    mpq_sub(result->a, result->a, temp);

    mpq_mul(temp, x->d, y->d);
    mpq_mul(temp, temp, two);
    mpq_sub(result->a, result->a, temp);

    mpq_mul(result->b, x->a, y->b);

    mpq_mul(temp, x->c, y->d);
    mpq_sub(result->b, result->b, temp);

    mpq_mul(temp, x->b, y->a);
    mpq_add(result->b, result->b, temp);

    mpq_mul(temp, x->d, y->c);
    mpq_sub(result->b, result->b, temp);

    mpq_mul(result->c, x->a, y->c);

    mpq_mul(temp, x->b, y->d);
    mpq_mul(temp, temp, two);
    mpq_add(result->c, result->c, temp);

    mpq_mul(temp, x->c, y->a);
    mpq_add(result->c, result->c, temp);

    mpq_mul(temp, x->d, y->b);
    mpq_mul(temp, temp, two);
    mpq_add(result->c, result->c, temp);

    mpq_mul(result->d, x->a, y->d);

    mpq_mul(temp, x->d, y->a);
    mpq_add(result->d, result->d, temp);

    mpq_mul(temp, x->b, y->c);
    mpq_add(result->d, result->d, temp);

    mpq_mul(temp, x->c, y->b);
    mpq_add(result->d, result->d, temp);
    
    mpq_clear(temp);
    mpq_clear(two);

    reduce(result);
    *x = *result;

    mpq_clear(result);

}

void
complexConjugate(qisq2_t *result, qisq2_t *x) {    
    mpq_set(result->a, x->a);
    mpq_set(result->b, x->b);
    mpq_neg(result->c, x->c);
    mpq_neg(result->d, x->d);
    
    reduce(result);
}

void
sqrttwoConjugate(qisq2_t *result, qisq2_t *x) {
    
    mpq_set(result->a, x->a);
    mpq_neg(result->b, x->b);
    mpq_set(result->c, x->c);
    mpq_neg(result->d, x->d);

    reduce(result);
}

void
weight_qisq2_div(qisq2_t *x, qisq2_t *y)
{
    qisq2_t *temp, *cc, *sc;
    init(temp); 
    init(cc); 
    init(sc);

    *temp = *y;

    complexConjugate(cc, temp);
    weight_qisq2_mul(temp, cc);
    sqrttwoConjugate(sc, temp);
    weight_qisq2_mul(temp, sc);

    // y is real fraction now
    assert(mpq_cmp_ui(temp->b, 0, 1) == 0 && mpq_cmp_ui(temp->c, 0, 1) == 0 && mpq_cmp_ui(temp->d, 0, 1) == 0);

    // multiply numerator x with conjugates
    weight_qisq2_mul(x, cc);
    weight_qisq2_mul(x, sc);

    clear(cc);
    clear(sc);

    // calculate 1/y
    mpq_inv(temp->a, temp->a);

    weight_qisq2_mul(x, temp);
    clear(temp);
    
    reduce(x);
}

bool
weight_qisq2_eq(qisq2_t *x, qisq2_t *y)
{
    if (mpq_equal(x->a,y->a) != 0){
        return false;
    }
    if (mpq_equal(x->b,y->b) != 0){
        return false;
    }
    if (mpq_equal(x->c,y->c) != 0){
        return false;
    }
    if (mpq_equal(x->c,y->c) != 0){
        return false;
    }
    return true;

}

bool
weight_qisq2_eps_close(qisq2_t *x, qisq2_t *y, double eps)
{
    return weight_qisq2_eq(x, y);
}

bool
weight_qisq2_greater(qisq2_t *x, qisq2_t *y)
{
    if (mpq_cmp(x->a,y->a) > 0){
        return false;
    }
    if (mpq_cmp(x->b,y->b) > 0){
        return false;
    }
    if (mpq_cmp(x->c,y->c) > 0){
        return false;
    }
    if (mpq_cmp(x->c,y->c) > 0){
        return false;
    }
    return true;
}

EVBDD_WGT
wgt_qisq2_norm_L2(EVBDD_WGT *low, EVBDD_WGT *high)
{
    // TODO
}

EVBDD_WGT
wgt_qisq2_get_low_L2normed(EVBDD_WGT high)
{
    // TODO
}

void
weight_qisq2_fprint(FILE *stream, qisq2_t *x)
{
    fprintf(stream, "\t%f + %f sqrt(2) + %f i + %f i sqrt(2)\n",
               mpq_get_d(x->a), mpq_get_d(x->b), mpq_get_d(x->c), mpq_get_d(x->d));
}

/*****************</Implementation of edge_weights interface>******************/
