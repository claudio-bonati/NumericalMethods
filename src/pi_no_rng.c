#include<math.h>
#include<stdio.h>
#include<stdlib.h>

// main
int main (int argc, char **argv)
    {
    int i, nop;
    double delta, x;
    double ris_u, ris_l; // ris_u=upper estimate, ris_l=lower estimate

    if(argc != 2)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s nop\n\n", argv[0]);
      fprintf(stdout, "  nop = number of points to be used in the discretization of the [0,1] interval\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  pi estimated by the rectangle method using n points, upper and lower bound\n");

      return EXIT_SUCCESS;
      }
    else
      {
      nop=atoi(argv[1]);
      }

    if(nop<0)
      {
      fprintf(stderr, "'nop' mast be positive\n");
      return EXIT_FAILURE;
      }

    ris_u=0.0;
    ris_l=0.0; 

    delta=1./(double)nop;

    // compute the integral of sqrt(1-x^2)
    for(i=0; i<nop; i++)
       {
       x=(double)i*delta;
       ris_u+=sqrt(1-x*x);

       x=(double)(i+1)*delta;
       ris_l+=sqrt(1-x*x);
       }
  
     ris_u*=delta;
     ris_l*=delta;

    printf("%.12lf %.12lf (accuracy=%g)\n", 4*ris_u, 4*ris_l, (ris_u-ris_l)/(ris_u+ris_l)*2.0);

    return EXIT_SUCCESS;
    }

