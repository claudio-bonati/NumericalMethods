#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/nvector.h"

// all functions are defined as inline functions in nvector.h

//initialize
void init(NVec * restrict a, double b[NCOMP]);


//initialize to zeros
void zeros(NVec * restrict a);


//initialize to one
void one(NVec * restrict a);


//assignement
void equal(NVec * restrict a, NVec const * const restrict b);


//times equal
void timesequal(NVec * restrict a, double d);


//plusequal
void plusequal(NVec * restrict a, NVec const * const restrict b);


//minusequal
void minusequal(NVec * restrict a, NVec const * const restrict b);


//sum
void sum(NVec * restrict a, NVec const * const restrict b, NVec const * const restrict c);


//scalar product
double scalprod(NVec const * const a, NVec const * const b);


//rotate two components by angle
void rotate2(NVec * restrict a, int i, int j, double phi);


//normalize
void normalize(NVec * restrict a);


//norm
double norm(NVec * restrict a);


// random vector
void randvec(NVec * restrict a);

