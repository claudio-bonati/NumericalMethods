#ifndef NVECTOR_H

#define NCOMP 2

typedef struct NVec {
     double comp[NCOMP];
} NVec;

//initialize
void init(NVec *a, double b[NCOMP]);

//initialize to zeros
void zeros(NVec *a);

//initialize to one in the 0 direction
void one(NVec *a);

//assignement
void equal(NVec *a, NVec const * const b);

//times equal
void timesequal(NVec *a, double d);

//plusequal
void plusequal(NVec *a, NVec const * const b);

//minusequal
void minusequal(NVec *a, NVec const * const b);

//sum
void sum(NVec *a, NVec const * const b, NVec const * const c);

//scalar product
double scalprod(NVec const * const a, NVec const * const b);

//rotate two components by angle
void rotate2(NVec *a, int i, int j, double phi);

//normalize
void normalize(NVec *a);

// norm
double norm(NVec *a);

#endif
