#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include "int_lib.h"


double test(double x)
    {
    return x*x;
    }


int main(void)
    {
    int err;
    double ris;
    const double prec=1.0e-7;
    const double pi=3.141592653589793238462643383279;

    printf("Rectangle method\n");

    ris=int_rel_rec(0.0, 1.0, &test, prec, &err); // this is declared in int_lib.h and defined in int_lib.c
    if(err==0)
      {
      printf("%.14lf\n", ris);
      }

    ris=int_abs_rec(0.0, 1.0, &test, prec, &err);  // this is declared in int_lib.h and defined in int_lib.c
    if(err==0)
      {
      printf("%.14lf\n", ris);
      }

    ris=int_rel_rec(0.0, 2.0*pi, &sin, prec, &err);  // this is declared in int_lib.h and defined in int_lib.c
    if(err==0)
      {
      printf("%.14lf\n", ris);
      }

    ris=int_abs_rec(0.0, 2.0*pi, &sin, prec, &err); // this is declared in int_lib.h and defined in int_lib.c
    if(err==0)
      {
      printf("%.14lf\n", ris);
      }
    
    return EXIT_SUCCESS;
    }


