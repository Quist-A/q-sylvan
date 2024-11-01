#include "wgt_storage_interface.h"


void * (*wgt_store_create)(uint64_t size, double tolerance);
void (*wgt_store_free)(void *wgt_storage);
int (*wgt_store_find_or_put)(const void *dbs, const void *v, uint64_t *ret);
void * (*wgt_store_get)(const void *dbs, const uint64_t ref);
uint64_t (*wgt_store_num_entries)(const void *dbs);
double (*wgt_store_get_tol)();


void init_wgt_storage_functions(wgt_storage_backend_t backend)
{

    switch (backend)
    {
    case COMP_HASHMAP:
        wgt_store_create      = &cmap_create;
        wgt_store_free        = &cmap_free;
        wgt_store_find_or_put = &cmap_find_or_put;
        wgt_store_get         = &cmap_get;
        wgt_store_num_entries = &cmap_count_entries;
        wgt_store_get_tol     = &cmap_get_tolerance;
        break;
    case QISQ2_MAP:
        wgt_store_create      = &qisq2_map_create;
        wgt_store_free        = &qisq2_map_free;
        wgt_store_find_or_put = &qisq2_map_find_or_put;
        wgt_store_get         = &qisq2_map_get;
        wgt_store_num_entries = &qisq2_map_count_entries;
        wgt_store_get_tol     = &qisq2_map_get_tolerance;
        break;
    default:
        break;
    }
}
