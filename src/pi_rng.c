#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/random.h"

// main
int main (int argc, char **argv)
    {
    long int i, sample, counter;
    double x, y, ris, sigma;

    if(argc != 2)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s sample\n", argv[0]);
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  pi estimated my simple sampling MC using 'sample' trials\n");

      return EXIT_SUCCESS;
      }
    else
      {
      sample=atol(argv[1]);
      }

    if(sample<0)
      {
      fprintf(stderr, "'sample' mast be positive\n");
      return EXIT_FAILURE;
      }

    // initialize random number generator
    myrand_init(2302342, 2312311);

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
   
    // the probability of falling inside the circle is pi/4 
    ris=(double)counter/(double) sample;
    ris*=4;

    // standard deviation of a sequence of 0, 1
    sigma=sqrt((double)counter/(double) sample - ((double)counter/(double) sample)*((double)counter/(double) sample))/sqrt(sample);
    sigma*=4;

    printf("%lf %lf (accuracy: %lf)\n", ris, sigma, sigma/ris);

    return EXIT_SUCCESS;
    }

