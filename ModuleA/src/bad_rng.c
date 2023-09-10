#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/random.h"


// random number generator state 
unsigned long int rng_state;


// RANDU: random number generator in [0,1)
// x_{i+1}=65539*x_i mod 2^{31}
// seed must be an odd number
//
//
// "its very name RANDU is enough to bring dismay into the
// eyes and stomachs of many computer scientists!" 
// D. Knuth "The art of computer programming" vol 2, third edition, page.107
//
double randu()
  {
  const unsigned long int const1=65539; 
  const unsigned long int const2=2147483648; //2^31
  unsigned long int y;

  y=const1*rng_state;
  rng_state=y % const2;

  return (double)rng_state / (double) const2;
  }


#define STRING_LENGTH 50

// main 
int main(int argc, char **argv)
    {
    unsigned long int seed;
    int j, dim;
    long i, maxiter;
    double x, sigma_x, tmp;
    char datafile[STRING_LENGTH];
    FILE *fp;

    if(argc != 5)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s seed dim sample datafile\n\n", argv[0]);
      fprintf(stdout, "  seed = seed for the random number generator (0 = machine time)\n");
      fprintf(stdout, "  dim = random points in (0,1)^{dim} are generated\n");
      fprintf(stdout, "  sample = number of points to be generated\n");
      fprintf(stdout, "  datafile = name of the file in which to print three columns of 'sample' random data\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  some statistical test\n\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 

      seed=(unsigned long int)atoi(argv[1]);
      dim=atoi(argv[2]);
      maxiter=atol(argv[3]);

      if(strlen(argv[4]) >= STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
        }
      else
        {
        strcpy(datafile, argv[4]);
        }
      }

    // seed=0 is changed to machine time 
    if(seed==0)
      {
      seed=(unsigned long int) time(NULL);
      }

    // seed must be an odd number
    if(seed%2==0)
      {
      seed+=1;
      }

    // initialize the random number generator 
    rng_state=seed;

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize counters
    x=0.0;

    // loop on iterations
    for(i=0; i<maxiter; i++)
       {
       for(j=0; j<dim; j++)
          {
          tmp=randu();
          x+=tmp;
          fprintf(fp, "%f ", tmp);
          }
       fprintf(fp,"\n");
       }

    // close data file
    fclose(fp);

    // normalize
    x/=((double) (maxiter)*(double) dim);

    //<x> = 1/2
    sigma_x=sqrt(1./3. - 1./4.)/sqrt((double) maxiter); // theoretical std
    printf("<x>-exact=%f ; ", x-1./2);
    printf("th_sigma=%.6f ; ", sigma_x);
    printf("(<x>-exact)/th_sigma=%f\n", (x-1./2.)/sigma_x);

    return EXIT_SUCCESS;
    }

