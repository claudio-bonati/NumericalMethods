#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/geometry.h"
#include"../include/random.h"

#define DIM 2     // dimensionality
#define NSTATES 8 // number of states
//#define METROPOLIS // if this macro is defined Metropolis is used, otherwise heatbath

#define STRING_LENGTH 50

// when using C99 M_PI is not defined in math.h header!
#ifndef M_PI
#  define M_PI  3.141592653589793238462643383279502884
#endif


// magnetization per site (clock order parameter)
double magn(int const * const restrict lattice, long int volume)
  {
  long int r;
  double tmp1, tmp2, sum1, sum2;
  const double instates=1.0/(double) NSTATES;

  sum1=0.0;
  sum2=0.0;
  for(r=0; r<volume; r++)
     {
     tmp1=cos(lattice[r]*2.0*M_PI*instates);
     tmp2=sin(lattice[r]*2.0*M_PI*instates);

     sum1+=tmp1;
     sum2+=tmp2;
     }

  sum1/=(double) volume;
  sum2/=(double) volume;

  return sqrt(sum1*sum1+sum2*sum2);
  }

/* Possible alternative definition of magnetization
// magnetization per site (check state 0, since b.c. do not favor any state) 
double magn(int const * const restrict lattice, long int volume)
  {
  long int r, sum;

  sum=0;
  for(r=0; r<volume; r++)
     {
     if(lattice[r]==0)
       {
       sum+=1;
       }
     }

  return ((double) NSTATES * ((double) sum/(double) volume ) -1.0)/((double) NSTATES -1.0);
  }
*/

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
        if(lattice[r]==lattice[nnp[dirgeo(r, i, volume)]])
          {
          sum--;
          }
        }
     }

  return (double) sum / (double) volume;
  }

// heatbath update at site r
// return 1 (for consistency with metropolis update)
//
// remember that aux_prob[i]=exp(-beta*((double)i));
int heatbath(int * restrict lattice, 
            long int r, 
            long int const * const restrict nnp, 
            long int const * const restrict nnm, 
            long int volume, 
            double const * const restrict aux_prob)
   {
   int i, nn[NSTATES];
   double tmp, prob[NSTATES];

   //nn[i] = numbers of sites with state "i"
   //the energy will be -nn[lattice[r]]
   for(i=0; i<NSTATES; i++)
      {
      nn[i]=0;
      }
   for(i=0; i<DIM; i++)
      {
      nn[lattice[nnp[dirgeo(r, i, volume)]]]+=1;
      nn[lattice[nnm[dirgeo(r, i, volume)]]]+=1;
      }

   tmp=0.0;
   for(i=0; i<NSTATES; i++)
      { 
      tmp+=1.0/aux_prob[nn[i]];
      prob[i]=tmp;
      }
   for(i=0; i<NSTATES; i++)
      { 
      prob[i]/=tmp;
      }

   tmp=myrand();
   for(i=0; i<NSTATES; i++)
      {
      if(tmp<prob[i])
        {
        lattice[r]=i;
        i=NSTATES;
        }
      }

   return 1;
   }


// metropolis update at site r
// return 1 if accepted, else 0
//
// remember that aux_prob[i]=exp(-beta*((double)i));
int metropolis(int * restrict lattice, 
               long int r, 
               long int const * const restrict nnp, 
               long int const * const restrict nnm, 
               long int volume, 
               double const * const restrict aux_prob)
   {
   int i, trial, nn[NSTATES];
   double acc_prob;

   //nn[i] = numbers of sites with state "i"
   //the energy will be -nn[lattice[r]]
   for(i=0; i<NSTATES; i++)
      {
      nn[i]=0;
      }
   for(i=0; i<DIM; i++)
      {
      nn[lattice[nnp[dirgeo(r, i, volume)]]]+=1;
      nn[lattice[nnm[dirgeo(r, i, volume)]]]+=1;
      }

   trial=(int) (NSTATES*myrand());
   acc_prob=aux_prob[nn[lattice[r]]]/aux_prob[nn[trial]];

   if(acc_prob>1)
     {
     lattice[r]=trial;
     return 1;
     }
   else if(myrand()<acc_prob)
          {
          lattice[r]=trial;
          return 1;
          }

   return 0;
   }


// main
int main(int argc, char **argv)
    {
    int i, L, *lattice;
    long int r, volume, sample, iter, acc; 
    long int *nnp, *nnm;
    double beta, locE, locM;
    double aux_prob[2*DIM+1];
  
    char datafile[STRING_LENGTH];
    FILE *fp;

    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    #ifdef METROPOLIS
      long int raux;
    #endif


    if(argc != 5)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s L beta sample datafile\n\n", argv[0]);
      fprintf(stdout, "  L = linear size of the lattice (dimension defined by macro)\n");
      fprintf(stdout, "  beta = inverse temperature\n");
      fprintf(stdout, "  sample = number of drawn to be extracted\n");
      fprintf(stdout, "  datafile = name of the file on which to write the data\n\n");
      fprintf(stdout, "Compiled for:\n");
      fprintf(stdout, "  dimensionality = %d\n", DIM);
      fprintf(stdout, "  number of states = %d\n", NSTATES);
      #ifdef METROPOLIS
        fprintf(stdout, "  update with Metropolis\n\n");
      #else
        fprintf(stdout, "  update with heatbath\n\n");
      #endif
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  E, M (E=energy per site, M=magnetization per site), one line for each draw\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 
      L=atoi(argv[1]);
      beta=atof(argv[2]);
      sample=atol(argv[3]);

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

    if(L<=0)
      {
      fprintf(stderr, "'L' must be positive\n");
      return EXIT_FAILURE;
      }

    if(sample<=0)
      {
      fprintf(stderr, "'sample' must be positive\n");
      return EXIT_FAILURE;
      }

    // initialize random number generator
    myrand_init(seed1, seed2);

    // compute the volume
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

    // initialize lattice to ordered start
    for(r=0; r<volume; r++)
       {
       lattice[r]=0;
       }

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize auxilliary vector for acceptance probability
    for(i=0; i<2*DIM+1; i++)
       {
       aux_prob[i]=exp(-beta*((double)i));
       }

    acc=0;
    for(iter=0; iter<sample; iter++)
       {
       for(r=0; r<volume; r++)
          {
          #ifdef METROPOLIS
            // metropolis
            raux=(long int)((double)volume * myrand());
            acc+=metropolis(lattice, raux, nnp, nnm, volume, aux_prob);
          #else
            //heatbath
            acc+=heatbath(lattice, r, nnp, nnm, volume, aux_prob);
          #endif
          }

       locE=energy(lattice, nnp, volume);
       locM=magn(lattice, volume);

       fprintf(fp, "%.12f %.12f\n", locE, locM);
       }

    // close datafile
    fclose(fp);

    printf("Acceptance rate %f\n", (double)acc / (double)sample / (double) volume);

    free(lattice);
    free(nnp);
    free(nnm);

    return EXIT_SUCCESS;
    }


