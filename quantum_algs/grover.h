#include <stdbool.h>

#include "sylvan.h"

/* Random bit array of lenght 'nbits' */
bool *qdd_grover_random_flag(BDDVAR nbits);

/**
 * Implementation of Grover where the gates are seen as functions applied to
 * the QDD.
 */
QDD qdd_grover(BDDVAR n, bool* flag);
#define qdd_grover_iteration(qdd,n,oracle) (CALL(qdd_grover_iteration,qdd,n,oracle));
TASK_DECL_3(QDD, qdd_grover_iteration, QDD, BDDVAR, bool*);

/**
 * Implementation of Grover where both the state vector and the gates are
 * represented as QDDs, and matrix-vector / matrix-matrix multiplication is
 * used to compute the result.
 */
QDD qdd_grover_matrix(BDDVAR n, bool *flag);
