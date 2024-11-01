#ifndef WGT_QISQ2_H
#define WGT_QISQ2_H

#include "sylvan_edge_weights.h"
#include "edge_weight_storage/flt.h"


/******************<Implementation of edge_weights interface>******************/

qisq2_t *weight_qisq2_malloc();
void _weight_qisq2_value(void *wgt_store, EVBDD_WGT a, qisq2_t *res);
EVBDD_WGT weight_qisq2_lookup(qisq2_t *a);
EVBDD_WGT _weight_qisq2_lookup_ptr(qisq2_t *a, void *wgt_store);

void init_qisq2_one_zero(void *wgt_store);

void weight_qisq2_abs(qisq2_t *a);
void weight_qisq2_neg(qisq2_t *a);
void weight_qisq2_conj(qisq2_t *a);
void weight_qisq2_sqr(qisq2_t *a);
void weight_qisq2_add(qisq2_t *a, qisq2_t *b);
void weight_qisq2_sub(qisq2_t *a, qisq2_t *b);
void weight_qisq2_mul(qisq2_t *a, qisq2_t *b);
void weight_qisq2_div(qisq2_t *a, qisq2_t *b);
bool weight_qisq2_eq(qisq2_t *a, qisq2_t *b);
bool weight_qisq2_eps_close(qisq2_t *a, qisq2_t *b, double eps);
bool weight_qisq2_greater(qisq2_t *a, qisq2_t *b);

EVBDD_WGT wgt_qisq2_norm_L2(EVBDD_WGT *low, EVBDD_WGT *high);
EVBDD_WGT wgt_qisq2_get_low_L2normed(EVBDD_WGT high);

void weight_qisq2_fprint(FILE *stream, qisq2_t *a);


static inline EVBDD_WGT
qisq2_lookup_angle(fl_t theta, fl_t mag)
{
	qisq2_t c; //= cmake_angle(theta, mag);
	return weight_lookup(&c);
}

static inline EVBDD_WGT
qisq2_lookup(fl_t r, fl_t i)
{
	qisq2_t c;// = cmake(r, i);
	return weight_lookup(&c);
}

/*****************</Implementation of edge_weights interface>******************/

#endif
