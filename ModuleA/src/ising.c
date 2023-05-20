#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/geometry.h"

#define DIM 2  // dimensionality

// magnetization per site
double magn(int const * const restrict lattice, long int volume)
  {
  int r, sum;

  sum=0;
  for(r=0; r<volume; r++)
     {
     sum+=lattice[r];
     }

  return (double) sum / (double) volume;
  }


// energy per site
double energy(int const * const restrict lattice, 
              long int const * const restrict nnp, 
              long int volume)
  {
  long int r, sum;
  int i;

  sum=0;
  for(r=0; r<volume; r++)
     {
     for(i=0; i<DIM; i++)
        {
        sum+=-lattice[r]*lattice[nnp[dirgeo(r, i, volume)]];
        }
     }

  return (double) sum / (double) volume;
  }


// main
int main(int argc, char **argv)
    {
    int i, L, *lattice;
    long int r1, r2, raux, volume, twotovol; 
    long int *nnp, *nnm;
    double beta, locE, locM, E, E2, M, Mabs, M2, M4, Z;

    if(argc != 3)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s L beta\n\n", argv[0]);
      fprintf(stdout, "  L = linear size of the lattice (dimension defined by macro)\n");
      fprintf(stdout, "  beta = inverse temperature\n\n");
      fprintf(stdout, "Compiled for:\n");
      fprintf(stdout, "  dimensionality = %d\n\n", DIM);
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  beta, <E>, <E^2>-<E>^2, <M>, <|M|>, <M^2>-<|M|>^2,");
      fprintf(stdout, " <M^2>, <M^4>/<M^2>^2, Z (E=energy per site, M=magnetization per site)\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 

      L=atoi(argv[1]);
      beta=atof(argv[2]);
      }

    volume=1;
    for(i=0; i<DIM; i++)
       {
       volume*=L;
       }

    // allocate the lattice (lexicographic order)
    // and next neighbors: nnp[dirgeo(r, i, volume)]= next neighbor in positive "i" direction of site r 
    lattice=(int *)malloc((unsigned long int)(volume)*sizeof(int));
    if(lattice == NULL)
      {
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    nnp=(long int *)malloc((unsigned long int)(DIM*volume)*sizeof(long int));
    if(nnp == NULL){
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    nnm=(long int *)malloc((unsigned long int)(DIM*volume)*sizeof(long int));
    if(nnm == NULL){
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize nnp and nnm
    init_neighbors(nnp, nnm, L, DIM);

    // twotovol=2^{volume}
    twotovol=1;
    for(r1=0; r1<volume; r1++)
       {
       twotovol*=2;
       }

    // initialize values
    E=0.0;
    E2=0.0;
    M=0.0;
    Mabs=0.0;
    M2=0.0;
    M4=0.0;
    Z=0.0;

    // enumerate all lattices and compute expectation values 
    for(r1=0; r1<twotovol; r1++)
       {
       raux=r1;

       for(r2=0; r2<volume; r2++)
          {
          if(raux %2 ==0)
            {
            lattice[r2]=1;
//            printf("0");
            }
          else
            {
            lattice[r2]=-1;
//            printf("1");
            }

          raux=raux>>1;
          }
//       printf("\n"); 

       locE=energy(lattice, nnp, volume);
       locM=magn(lattice, volume);

       E+=locE*exp(-beta * locE * (double)volume);
       E2+=locE*locE*exp(-beta * locE * (double)volume);
       M+=locM*exp(-beta * locE * (double)volume);
       Mabs+=fabs(locM)*exp(-beta * locE * (double)volume);
       M2+=locM*locM*exp(-beta * locE * (double)volume);
       M4+=pow(locM, 4.0)*exp(-beta * locE * (double)volume);
       Z+=exp(-beta * locE * (double)volume);
       }

    // normalize
    E/=Z;
    E2/=Z;
    M/=Z;
    Mabs/=Z;
    M2/=Z;
    M4/=Z;

    printf("%f %f %f %f %f %f %f %f %f\n", beta, E, E2-E*E, M, Mabs, M2-Mabs*Mabs, M2, M4/(M2*M2), Z);

    free(lattice);
    free(nnp);
    free(nnm);

    return EXIT_SUCCESS;
    }


