#include<math.h>

#include"../include/random.h"

// normal gaussian random number generator 
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
void gauss2(double *ris1, double *ris2)
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


