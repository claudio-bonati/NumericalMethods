#include<math.h>

#include"../include/random.h"

// when using C99 M_PI is not defined in math.h header!
#ifndef M_PI
#  define M_PI  3.141592653589793238462643383279502884
#endif

// normal gaussian random number generator 
// polar form 
double gauss1(void)
  {
  double v1, v2, s, ris;

  do
    {
    v1=1.0-2.0*myrand();
    v2=1.0-2.0*myrand();
    s=v1*v1+v2*v2;
    }
  while(s > 1);

  ris=v1*sqrt(-2*log(s)/s);
  return ris;
  }


// normal gaussian random number generator
// polar form
void gauss2(double * restrict ris1, double * restrict ris2)
  {
  double v1, v2, s;

  do
    {
    v1=-1.0+2.0*myrand();
    v2=-1.0+2.0*myrand();
    s=v1*v1+v2*v2;
    }
  while(s > 1);

  *ris1=v1*sqrt(-2*log(s)/s);
  *ris2=v2*sqrt(-2*log(s)/s);
  }


// couple of independent normal random numbers
// basic form
void gauss2_basic(double * restrict ris1, double * restrict ris2)
  {
  double R, theta;

  R=sqrt(-2.0*log(myrand()));
  theta=2.0*M_PI*myrand();

  *ris1=R*cos(theta);
  *ris2=R*sin(theta);
  }

