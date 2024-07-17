#include<math.h>
#include<stdio.h>
#include<stdlib.h>

// compute the integral \int_{xmin}^{xmax) f(x) dx 
// with relative accuracy relacc_target
// using the rectangle method
double int_rel_rec(double xmin, 
                   double xmax, 
                   double (*f)(double), // function pointer
                   double relacc_target,
                   int *err)
    {
    int iteration;
    const int maxiter=25;
    long int i, steps;
    double delta, ris, risold, absacc;

    iteration=1;
    steps=10;
    delta=(xmax-xmin)/(double) steps;
    
    ris=0.0;
    for(i=0; i<steps; i++)
       {
       ris+=f(xmin+(double)i*delta);
       }
    ris*=delta;

    absacc=relacc_target*ris+1;
    while(iteration<maxiter && absacc > relacc_target*ris)
         {
         iteration++;
         risold=ris;
  
         steps*=2;
         delta=(xmax-xmin)/(double)steps;

         ris=0.0;
         for(i=0; i<steps; i++)
            {
            ris+=f(xmin+(double)i*delta);
            }
         ris*=delta;

         absacc=fabs(ris-risold);
         }

    if(iteration==maxiter)
      {
      fprintf(stderr, "  Maximum iteration reached (%s, %d)\n", __FILE__, __LINE__);
      *err=1;
      }
    else
      {
      *err=0;
      }

    return ris;
    } 
   
 
// compute the integral \int_{xmin}^{xmax) f(x) dx 
// with  accuracy absacc_target
// using the rectangle method
double int_abs_rec(double xmin, 
                   double xmax, 
                   double (*f)(double), // function pointer
                   double absacc_target,
                   int *err)
    {
    int iteration;
    const int maxiter=25;
    long int i, steps;
    double delta, ris, risold, absacc;

    iteration=1;
    steps=10;
    delta=(xmax-xmin)/(double)steps;
    
    ris=0.0;
    for(i=0; i<steps; i++)
       {
       ris+=f(xmin+(double)i*delta);
       }
    ris*=delta;

    absacc=absacc_target+1;
    while(iteration<maxiter && absacc > absacc_target)
         {
         iteration++;
         risold=ris;
  
         steps*=2;
         delta=(xmax-xmin)/(double)steps;

         ris=0.0;
         for(i=0; i<steps; i++)
            {
            ris+=f(xmin+(double)i*delta);
            }
         ris*=delta;

         absacc=fabs(ris-risold);
         }

    if(iteration==maxiter)
      {
      fprintf(stderr, "  Maximum iteration reached (%s, %d)\n", __FILE__, __LINE__);
      *err=1;
      }
    else
      { 
      *err=0;
      }

    return ris;
    } 
 
