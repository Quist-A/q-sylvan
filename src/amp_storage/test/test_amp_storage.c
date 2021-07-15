#include <stdio.h>

#include "test_assert.h"
#include "../amp_storage_interface.h"

bool VERBOSE = true;

static complex_t
comp_make(fl_t r, fl_t i)
{
    complex_t res;
    res.r = r;
    res.i = i;
    return res;
}

int test_cmap()
{
    void *ctable = cmap_create(1<<10, 1e-14);

    ref_t index1, index2;
    complex_t val1, val2, val3;
    int found;

    val1 = comp_make(3.5, 4.7);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    for(int k=0; k<10; k++){
        found = cmap_find_or_put(ctable, &val1, &index2);
        test_assert(found == 1);
        test_assert(index2 == index2);
    }

    val1 = comp_make(0.9, 2./3.);
    val2 = comp_make(0.9, 2./3.);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);

    val1 = comp_make(1.0/flt_sqrt(2.0),0); // 1/sqrt(2)
    val2 = comp_make(1.0/flt_sqrt(2.0),0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = cmap_get(ctable, index1);
    test_assert(flt_abs(val3.r - val1.r) < cmap_get_tolerance());
    test_assert(flt_abs(val3.i - val1.i) < cmap_get_tolerance());

    val1 = comp_make(2.99999999999999855, 0.0);
    val2 = comp_make(3.00000000000000123, 0.0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = cmap_get(ctable, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);

    val1 = comp_make(0.0005000000000012, 0.0);
    val2 = comp_make(0.0004999999999954, 0.0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = cmap_get(ctable, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);


    // test with tolerance = 0
    cmap_free(ctable);
    ctable = cmap_create(1<<10, 0.0);

    val1 = comp_make(0.9, 2./3.);
    val2 = comp_make(0.9, 2./3.);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);

    val1 = comp_make(2.99999999999999855, 0.0);
    val2 = comp_make(3.00000000000000123, 0.0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 0);
    test_assert(index1 != index2);
    val3 = cmap_get(ctable, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);

    cmap_free(ctable);
    if(VERBOSE) printf("cmap tests:               ok\n");
    return 0;
}

int test_rmap()
{
    void *rtable = rmap_create(1<<10, 1e-14);

    ref_t index1, index2;
    fl_t val1, val2, val3;
    int found;

    val1 = 3.5;
    found = rmap_find_or_put(rtable, &val1, &index1); test_assert(found == 0);
    for(int k=0; k<10; k++){
        found = rmap_find_or_put(rtable, &val1, &index2);
        test_assert(found == 1);
        test_assert(index2 == index2);
    }

    val1 = (2./3.);
    val2 = (2./3.);
    found = rmap_find_or_put(rtable, &val1, &index1); test_assert(found == 0);
    found = rmap_find_or_put(rtable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);


    val1 = 1.0/flt_sqrt(2.0);
    val2 = 1.0/flt_sqrt(2.0);
    found = rmap_find_or_put(rtable, &val1, &index1); test_assert(found == 0);
    found = rmap_find_or_put(rtable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *rmap_get(rtable, index1);
    test_assert((val3 - val1) < rmap_get_tolerance());

    val1 = 2.99999999999999855;
    val2 = 3.00000000000000123;
    found = rmap_find_or_put(rtable, &val1, &index1); test_assert(found == 0);
    found = rmap_find_or_put(rtable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *rmap_get(rtable, index1);
    test_assert(val3 == val1);

    val1 = 0.0005000000000012;
    val2 = 0.0004999999999954;
    found = rmap_find_or_put(rtable, &val1, &index1); test_assert(found == 0);
    found = rmap_find_or_put(rtable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *rmap_get(rtable, index1);
    test_assert(val3 == val1);

    // test with tolerance = 0
    rmap_free(rtable);
    rtable = rmap_create(1<<10, 0.0);

    val1 = (2./3.);
    val2 = (2./3.);
    found = rmap_find_or_put(rtable, &val1, &index1); test_assert(found == 0);
    found = rmap_find_or_put(rtable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);

    val1 = 2.99999999999999855;
    val2 = 3.00000000000000123;
    found = rmap_find_or_put(rtable, &val1, &index1); test_assert(found == 0);
    found = rmap_find_or_put(rtable, &val2, &index2); test_assert(found == 0);
    test_assert(index1 != index2);
    val3 = *rmap_get(rtable, index1);
    test_assert(val3 == val1);

    rmap_free(rtable);
    if(VERBOSE) printf("rmap tests:               ok\n");
    return 0;
}

int test_tree_map()
{
    void *tree_map = tree_map_create(1<<10, 1e-14);

    unsigned long index1, index2;
    fl_t val1, val2, val3;
    int found;

    val1 = 3.5;
    found = tree_map_find_or_put(tree_map, val1, &index1); test_assert(found == 0);
    for(int k=0; k<10; k++){
        found = tree_map_find_or_put(tree_map, val1, &index2);
        test_assert(found == 1);
        test_assert(index2 == index2);
    }

    val1 = (2./3.);
    val2 = (2./3.);
    found = tree_map_find_or_put(tree_map, val1, &index1); test_assert(found == 0);
    found = tree_map_find_or_put(tree_map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);

    val1 = 1.0/flt_sqrt(2.0);
    val2 = 1.0/flt_sqrt(2.0);
    found = tree_map_find_or_put(tree_map, val1, &index1); test_assert(found == 0);
    found = tree_map_find_or_put(tree_map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *tree_map_get(tree_map, index1);
    test_assert((val3 - val1) < tree_map_get_tolerance());

    val1 = 2.99999999999999855;
    val2 = 3.00000000000000123;
    found = tree_map_find_or_put(tree_map, val1, &index1); test_assert(found == 0);
    found = tree_map_find_or_put(tree_map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *tree_map_get(tree_map, index1);
    test_assert(val3 == val1);

    val1 = 0.0005000000000012;
    val2 = 0.0004999999999954;
    found = tree_map_find_or_put(tree_map, val1, &index1); test_assert(found == 0);
    found = tree_map_find_or_put(tree_map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *tree_map_get(tree_map, index1);
    test_assert(val3 == val1);

    val1 = 14.2;
    val2 = 14.25;
    found = tree_map_find_or_put(tree_map, val1, &index1); test_assert(found == 0);
    found = tree_map_find_or_put(tree_map, val2, &index2); test_assert(found == 0);
    test_assert(index1 != index2);

    // test with tolerance = 0
    tree_map_free(tree_map);
    tree_map = tree_map_create(1<<10, 0.0);

    val1 = (2./3.);
    val2 = (2./3.);
    found = tree_map_find_or_put(tree_map, val1, &index1); test_assert(found == 0);
    found = tree_map_find_or_put(tree_map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);

    val1 = 2.99999999999999855;
    val2 = 3.00000000000000123;
    found = tree_map_find_or_put(tree_map, val1, &index1); test_assert(found == 0);
    found = tree_map_find_or_put(tree_map, val2, &index2); test_assert(found == 0);
    test_assert(index1 != index2);
    val3 = *tree_map_get(tree_map, index1);
    test_assert(val3 == val1);

    tree_map_free(tree_map);
    if(VERBOSE) printf("tree map tests:           ok\n");
    return 0;
}


int runtests()
{
    if (test_cmap()) return 1;
    if (test_rmap()) return 1;
    if (test_tree_map()) return 1;
    return 0;
}

int main()
{
    return runtests();
}