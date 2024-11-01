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
    // TODO
}

void
weight_qisq2_abs(qisq2_t *a)
{
    // TODO
}

void
weight_qisq2_neg(qisq2_t *a)
{
    // TODO
}

void
weight_qisq2_conj(qisq2_t *a)
{
    // TODO
}

void
weight_qisq2_sqr(qisq2_t *a)
{
    weight_qisq2_mul(a, a);
}

void
weight_qisq2_add(qisq2_t *a, qisq2_t *b)
{
    // TODO
}

void
weight_qisq2_sub(qisq2_t *a, qisq2_t *b)
{
    // TODO
}

void
weight_qisq2_mul(qisq2_t *a, qisq2_t *b)
{
    // TODO
}

void
weight_qisq2_div(qisq2_t *a, qisq2_t *b)
{
    // TODO
}

bool
weight_qisq2_eq(qisq2_t *a, qisq2_t *b)
{
    // TODO
}

bool
weight_qisq2_eps_close(qisq2_t *a, qisq2_t *b, double eps)
{
    // TODO
}

bool
weight_qisq2_greater(qisq2_t *a, qisq2_t *b)
{
    // TODO
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
weight_qisq2_fprint(FILE *stream, qisq2_t *a)
{
    // TODO
}

/*****************</Implementation of edge_weights interface>******************/
