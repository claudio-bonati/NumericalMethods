#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include"../include/random.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

// main
int main(int argc, char **argv)
    {
    long int i, sample, counter;
    double x, y, ris, sigma;
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    if(argc != 2)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s sample\n\n", argv[0]);
      fprintf(stdout, "  sample = number of draws to be used\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  pi estimated my simple sampling MC using 'sample' draws\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input

      sample=atol(argv[1]);
      }

    if(sample<=0)
      {
      fprintf(stderr, "'sample' must be positive\n");
      return EXIT_FAILURE;
      }

    // initialize random number generator
    myrand_init(seed1, seed2);

    // count how many points fall inside the circle of radius 1
    counter=0;
    for(i=0; i<sample; i++)
       {
       x=myrand();
       y=myrand();

       // if the point is inside the circle
       if(x*x+y*y<1)
         {
         counter+=1;
         }
       }
   
    ris=(double)counter/(double) sample;

    // standard deviation of the mean a sequence of 0, 1
    sigma=sqrt(ris - ris*ris)/sqrt((double)sample-1.0);

    // the probability of falling inside the circle is pi/4 
    ris*=4;
    sigma*=4;

    printf("estimate-pi=%f ; ", ris-M_PI);
    printf("sigma=%f ; ", sigma);
    printf("(estimate-pi)/sigma=%f\n", (ris-M_PI)/sigma);

    return EXIT_SUCCESS;
    }

