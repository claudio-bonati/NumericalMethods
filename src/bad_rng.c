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
int main (int argc, char **argv)
    {
    unsigned long int seed;
    int i, maxiter;
    double x, x2, xy, xyz;
    double sigma_x, sigma_x2, sigma_xy, sigma_xyz;
    double tmp1, tmp2, tmp3;
    char datafile[STRING_LENGTH];
    FILE *fp;

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s seed sample datafile\n\n", argv[0]);
      fprintf(stdout, "  seed = seed for the random number generator (0 = machine time)\n");
      fprintf(stdout, "  sample = number of points to be generated\n");
      fprintf(stdout, "  datafile = name of the file in which to print three columns of 'sample' random data\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  some statistical test\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values

      seed=(unsigned long int)atoi(argv[1]);
      maxiter=atoi(argv[2]);

      if(strlen(argv[3]) >= STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
        }
      else
        {
        strcpy(datafile, argv[3]);
        }
      }

    // seed=0 is changed to machine time
    if(seed==0)
      {
      seed=(unsigned long int) time(NULL);
      }

    //initialize te random number generator
    rng_state=seed;

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // intialize average values
    x=0.0;
    x2=0.0;
    xy=0.0;
    xyz=0.0;

    // loop on iterations
    for(i=0; i<maxiter; i++)
       {
       tmp1=randu();
       tmp2=randu();
       tmp3=randu();

       fprintf(fp, "%lf %lf %lf\n", tmp1, tmp2, tmp3);

       x+=tmp1;
       x2+=tmp1*tmp1;
       xy+=tmp1*tmp2;
       xyz+=tmp1*tmp2*tmp3;
       }

    // close data file
    fclose(fp);

    // normalize
    x/=(double) (maxiter);
    x2/=(double) (maxiter);
    xy/=(double) (maxiter);
    xyz/=(double) (maxiter);

    //<x> = 1/2
    sigma_x=sqrt(1./3. - 1./4.)/sqrt(maxiter);
    printf("<x>-exact=%lf ; ", x-1./2);
    printf("th_sigma=%.6lf ; ", sigma_x);
    printf("(<x>-exact)/th_sigma=%lf\n", (x-1./2.)/sigma_x);

    //<x^2> = 1/3
    sigma_x2=sqrt(1./5. - 1./9.)/sqrt(maxiter);
    printf("<x2>-exact=%.6lf ; ", x2-1./3.);
    printf("th_sigma=%lf ; ", sigma_x2);
    printf("(<x2>-exact)/th_sigma=%lf\n", (x2-1./3.)/sigma_x2);

    //<xy> = 1/4
    sigma_xy=sqrt(1./9. - 1./16.)/sqrt(maxiter);
    printf("<xy>-exact=%.6lf ; ", xy-1./4.);
    printf("th_sigma=%lf ; ", sigma_xy);
    printf("(<xy>-exact)/th_sigma=%lf\n", (xy-1./4.)/sigma_xy);

    //<xyz> = 1/8
    sigma_xyz=sqrt(1./27. - 1./64.)/sqrt(maxiter);
    printf("<xyz>-exact=%.6lf ; ", xyz - 1./8.);
    printf("th_sigma=%lf ; ", sigma_xyz);
    printf("(<xyz>-exact)/th_sigma=%lf\n", (xyz-1./8.)/sigma_xyz);

    return EXIT_SUCCESS;
    }

