#include "qisq2_map.h"

#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "atomics.h"
#include "fast_hash.h"
#include "util.h"

#undef CACHE_LINE
#undef CACHE_LINE_SIZE

#define CACHE_LINE 8
#define CACHE_LINE_SIZE 256

// how many "blocks" of 64 bits for a single table entry
#define entry_size (4*sizeof(uint64_t)/8)

typedef union {
    qisq2_t         c;
    uint64_t        d[entry_size];
} bucket_t;

// float "equality" tolerance
static long double TOLERANCE = 1e-14l;
static const uint64_t EMPTY = 14738995463583502973ull;
static const uint64_t LOCK  = 14738995463583502974ull;
static const uint64_t CL_MASK = -(1ULL << CACHE_LINE);

/**
\typedef Lockless hastable database.
*/
typedef struct qisq2_map_s qisq2_map_t;
struct qisq2_map_s {
    size_t              size;
    size_t              mask;
    size_t              threshold;
    int                 seen_0;
    bucket_t  __attribute__(( __aligned__(32)))       *table;
    // Q: should this 32 change to 16 now that we use doubles instead of
    // long doubles for the real and imaginary components?
};

static void __attribute__((unused))
print_bucket_bits(bucket_t* b)
{
    printf("%016" PRIu64, b->d[0]);
    for (unsigned int k = 1; k < entry_size; k++) {
        printf(" %016" PRIu64, b->d[k]);
    }
    printf("\n");
}


int
qisq2_map_find_or_put(const void *dbs, const void *_v, uint64_t *ret)
{
    qisq2_t *v = (qisq2_t *)_v;
    qisq2_map_t *qisq2_map = (qisq2_map_t *) dbs;
    bucket_t *val  = (bucket_t *) v;

    // TODO: compute hash based on values (not on pointers)
    uint32_t hash = 0;
    //uint32_t hash  = SuperFastHash(&round_v, sizeof(complex_t), 0);
    uint32_t prime = odd_primes[hash & PRIME_MASK];

    assert (val->d[0] != LOCK);
    assert (val->d[0] != EMPTY);

    // Insert/lookup `v`
    for (unsigned int c = 0; c < qisq2_map->threshold; c++) {
        uint64_t            ref = hash & qisq2_map->mask;
        uint64_t            line_end = (ref & CL_MASK) + CACHE_LINE_SIZE;
        for (size_t i = 0; i < CACHE_LINE_SIZE; i++) {
            
            // 1. Get bucket
            bucket_t *bucket = &qisq2_map->table[ref];

            // 2. If bucket empty, insert new value here
            if (bucket->d[0] == EMPTY) {
                if (cas(&bucket->d[0], EMPTY, LOCK)) {
                    *ret = ref;
                    // write backwards (overwrite bucket->d[0] last)
                    for (int k = entry_size-1; k >= 0; k--) {
                        atomic_write (&bucket->d[k], val->d[k]);
                    }
                    return 0;
                }
            }

            // 3. Bucket not empty, wait for lock
            while (atomic_read(&bucket->d[0]) == LOCK) {}

            // 4. Bucket contains some complex value, check if close to `v`
            qisq2_t *in_table = (qisq2_t *)bucket;
            if (false) {
                // TODO: decide if values equal
                *ret = ref;
                return 1;
            }

            // If unsuccessful, try next
            ref += 1;
            ref = ref == line_end ? line_end - CACHE_LINE_SIZE : ref;
        }
        hash += prime << CACHE_LINE;
    }
    // amplitude table full, unable to add
    return -1;
}

void *
qisq2_map_get(const void *dbs, const uint64_t ref)
{
    qisq2_map_t *qisq2_map = (qisq2_map_t *) dbs;
    return &(qisq2_map->table[ref].c);
}

uint64_t
qisq2_map_count_entries(const void *dbs)
{
    qisq2_map_t *qisq2_map = (qisq2_map_t *) dbs;
    uint64_t entries = 0;
    for (unsigned int c = 0; c < qisq2_map->size; c++) {
        if (qisq2_map->table[c].d[0] != EMPTY)
            entries++;
    }
    return entries;
}

void
qisq2_map_print_bitvalues(const void *dbs, const uint64_t ref)
{
    qisq2_map_t *qisq2_map = (qisq2_map_t *) dbs;
    bucket_t *b = (bucket_t *) qisq2_map_get(qisq2_map, ref);
    printf("%016" PRIu64, b->d[0]);
    for (unsigned int k = 1; k < entry_size; k++) {
        printf(" %016" PRIu64, b->d[k]);
    }
}

void *
qisq2_map_create(uint64_t size, double tolerance)
{
    qisq2_map_t  *qisq2_map = calloc (1, sizeof(qisq2_map_t));
    qisq2_map->size = size;
    qisq2_map->mask = qisq2_map->size - 1;
    qisq2_map->table = calloc (qisq2_map->size, sizeof(bucket_t));
    for (unsigned int c = 0; c < qisq2_map->size; c++) {
        qisq2_map->table[c].d[0] = EMPTY;
    }
    qisq2_map->threshold = qisq2_map->size / 100;
    qisq2_map->threshold = min(qisq2_map->threshold, 1ULL << 16);
    qisq2_map->seen_0 = 0;
    return (void *) qisq2_map;
}

void
qisq2_map_free(void *dbs)
{
    qisq2_map_t * qisq2_map = (qisq2_map_t *) dbs;
    free (qisq2_map->table);
    free (qisq2_map);
}

double
qisq2_map_get_tolerance()
{
    return 0;
}
