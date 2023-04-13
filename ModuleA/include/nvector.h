#ifndef NVECTOR_H

#define NCOMP 2

typedef struct NVec {
     double comp[NCOMP];
} NVec;

//initialize
void init(NVec * restrict a, double b[NCOMP]);

//initialize to zeros
void zeros(NVec * restrict a);

//initialize to one in the 0 direction
void one(NVec *a);

//assignement
void equal(NVec * restrict a, NVec const * const restrict b);

//times equal
void timesequal(NVec *restrict a, double d);

//plusequal
void plusequal(NVec * restrict a, NVec const * const restrict b);

//minusequal
void minusequal(NVec * restrict a, NVec const * const restrict b);

//sum
void sum(NVec * restrict a, NVec const * const restrict b, NVec const * const restrict c);

//scalar product
double scalprod(NVec const * const restrict a, NVec const * const restrict b);

//rotate two components by angle
void rotate2(NVec * restrict a, int i, int j, double phi);

//normalize
void normalize(NVec * restrict a);

// norm
double norm(NVec * restrict a);

// random vector
void randvec(NVec * restrict a);

#endif
