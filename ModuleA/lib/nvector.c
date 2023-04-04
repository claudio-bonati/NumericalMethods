#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/nvector.h"

//initialize
void init(NVec *a, double b[NCOMP])
  {
  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]=b[i];
     }
  }


//initialize to zeros
void zeros(NVec *a)
  {
  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]=0.0;
     }
  }


//initialize to one
void one(NVec *a)
  {
  int i;

  a->comp[0]=1.0;
  for(i=1; i<NCOMP; i++)
     {
     a->comp[i]=0.0;
     }
  }


//assignement
void equal(NVec *a, NVec const * const b)
  {
  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]=b->comp[i];
     }
  }


//times equal
void timesequal(NVec *a, double d)
  {
  int i;

  for(i=0; i<NCOMP; i++)
     {
     (a->comp[i])*=d;
     }
  }


//plusequal
void plusequal(NVec *a, NVec const * const b)
  {
  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]+=b->comp[i];
     }
  }


//minusequal
void minusequal(NVec *a, NVec const * const b)
  {
  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]-=b->comp[i];
     }
  }


//sum
void sum(NVec *a, NVec const * const b, NVec const * const c)
  {
  int i;

  for(i=0; i<NCOMP; i++)
     {
     a->comp[i]=b->comp[i]+c->comp[i];
     }
  }

//scalar product
double scalprod(NVec const * const a, NVec const * const b)
  {
  int i;
  double ris=0.0;

  for(i=0; i<NCOMP; i++)
     {
     ris+=(a->comp[i])*(b->comp[i]);
     }

  return ris;
  }


//rotate two components by angle
void rotate2(NVec *a, int i, int j, double phi)
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
void normalize(NVec *a)
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
double norm(NVec *a)
  {
  int i;
  double norm=0.0;
  
  for(i=0; i<NCOMP; i++)
     {
     norm+=(a->comp[i])*(a->comp[i]);
     }
  return sqrt(norm);
  }


