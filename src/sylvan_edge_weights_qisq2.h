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
void weight_qisq2_complexConjugate(qisq2_t *result, qisq2_t *x);
void weight_qisq2_sqrttwoConjugate(qisq2_t *result, qisq2_t *x);

void weight_qisq2_fprint(FILE *stream, qisq2_t *a);

void qisq2_init(qisq2_t *num);
void qisq2_clear(qisq2_t *num);
void qisq2_reduce(qisq2_t *num);

static inline qisq2_t qisq2_make(int64_t anum, int64_t aden, int64_t bnum, int64_t bden, int64_t cnum, int64_t cden, int64_t dnum, int64_t dden)
{
    qisq2_t res;
	qisq2_init(&res);
	mpq_set_si(res.a, anum, aden);
    mpq_set_si(res.b, bnum, bden);
    mpq_set_si(res.c, cnum, cden);
    mpq_set_si(res.d, dnum, dden);
	qisq2_reduce(&res);
    return res;
}

static inline EVBDD_WGT
qisq2_lookup(int64_t anum, int64_t aden, int64_t bnum, int64_t bden, int64_t cnum, int64_t cden, int64_t dnum, int64_t dden)
{
	qisq2_t c;
	c = qisq2_make(anum, aden, bnum, bden, cnum, cden, dnum, dden);
	//TODO: memory leakage is here, because c is initialized but not cleared
	return weight_lookup(&c);
}

static inline qisq2_t qisq2_zero() { return qisq2_make( 0, 1, 0, 1, 0, 1, 0, 1); }
static inline qisq2_t qisq2_one()  { return qisq2_make( 1, 1, 0, 1, 0, 1, 0, 1); }
static inline qisq2_t qisq2_mone() { return qisq2_make(-1, 1, 0, 1, 0, 1, 0, 1); }

/*****************</Implementation of edge_weights interface>******************/

#endif
