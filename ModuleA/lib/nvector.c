#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/nvector.h"
#include"../include/random.h"

//initialize
void init(NVec * restrict a, double b[NCOMP])
  {
  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]=b[i];
     }
  }


//initialize to zeros
void zeros(NVec * restrict a)
  {
  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]=0.0;
     }
  }


//initialize to one
void one(NVec * restrict a)
  {
  int i;

  a->comp[0]=1.0;
  for(i=1; i<NCOMP; i++)
     {
     a->comp[i]=0.0;
     }
  }


//assignement
void equal(NVec * restrict a, NVec const * const restrict b)
  {
  #ifdef DEBUG
  if(a==b)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]=b->comp[i];
     }
  }


//times equal
void timesequal(NVec * restrict a, double d)
  {
  int i;

  for(i=0; i<NCOMP; i++)
     {
     (a->comp[i])*=d;
     }
  }


//plusequal
void plusequal(NVec * restrict a, NVec const * const restrict b)
  {
  #ifdef DEBUG
  if(a==b)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]+=b->comp[i];
     }
  }


//minusequal
void minusequal(NVec * restrict a, NVec const * const restrict b)
  {
  #ifdef DEBUG
  if(a==b)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]-=b->comp[i];
     }
  }


//sum
void sum(NVec * restrict a, NVec const * const restrict b, NVec const * const restrict c)
  {
  #ifdef DEBUG
  if(a==b || a==c || b==c)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]=b->comp[i]+c->comp[i];
     }
  }

//scalar product
double scalprod(NVec const * const restrict a, NVec const * const restrict b)
  {
  #ifdef DEBUG
  if(a==b)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i;
  double ris=0.0;

  for(i=0; i<NCOMP; i++)
     {
     ris+=(a->comp[i])*(b->comp[i]);
     }

  return ris;
  }


//rotate two components by angle
void rotate2(NVec * restrict a, int i, int j, double phi)
  {
  #ifdef DEBUG
  if(i==j)
    {
    fprintf(stderr, "error in rotate2 (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif 

  double v1, v2;

  v1=a->comp[i];
  v2=a->comp[j];

  a->comp[i]= v1*cos(phi)+v2*sin(phi);
  a->comp[j]=-v1*sin(phi)+v2*cos(phi);
  }


//normalize
void normalize(NVec * restrict a)
  {
  int i;
  double norm=0.0;
  
  for(i=0; i<NCOMP; i++)
     {
     norm+=(a->comp[i])*(a->comp[i]);
     }
  norm=1.0/sqrt(norm);
 
  for(i=0; i<NCOMP; i++)
     {
     (a->comp[i])*=norm;
     }
  }


//norm
double norm(NVec * restrict a)
  {
  int i;
  double norm=0.0;
  
  for(i=0; i<NCOMP; i++)
     {
     norm+=(a->comp[i])*(a->comp[i]);
     }
  return sqrt(norm);
  }


// random vector
void randvec(NVec * restrict a)
  {
  int i;
  double norma;

  do{
    for(i=0; i<NCOMP; i++)
       {
       a->comp[i]=1.0-2.0*myrand();
       }

    norma=norm(a);
    }
  while(norma>1);

  norma=1./norma;
  timesequal(a, norma);
  }


