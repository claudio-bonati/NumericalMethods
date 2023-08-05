#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/geometry.h"
#include"../include/nvector.h"
#include"../include/random.h"

#define DIM 2 // dimensionality
#define STRING_LENGTH 50


// magnetization per site
double magn(NVec const * const restrict lattice, long int volume)
  {
  long int r;
  NVec S;

  zeros(&S);
  for(r=0; r<volume; r++)
     {
     plusequal(&S, &(lattice[r]) );
     }
  timesequal(&S, 1./(double)volume);

  return  norm(&S);
  }


// energy per site
double energy(NVec const * const restrict lattice, 
              long int const * const restrict nnp, 
              long int volume)
  {
  long int r;
  int i;
  double sum;

  sum=0;
  for(r=0; r<volume; r++)
     {
     for(i=0; i<DIM; i++)
        {
        sum+=-scalprod(&(lattice[r]), &(lattice[nnp[dirgeo(r, i, volume)]]));
        }
     }

  return (double) sum / (double) volume;
  }


// metropolis update at site r
// return 1 if accepted, else 0
int metropolis(NVec * restrict lattice, 
               long int r, 
               long int const * const restrict nnp, 
               long int const * const restrict nnm, 
               long int volume, 
               double beta,
               double phimax) 
  {
  NVec sumnn, trial;
  int i, j, acc=0;
  double phi, deltaE;

  zeros(&sumnn);
  for(i=0; i<DIM; i++)
     {
     plusequal(&sumnn, &(lattice[nnp[dirgeo(r, i, volume)]]) );
     plusequal(&sumnn, &(lattice[nnm[dirgeo(r, i, volume)]]) );
     }

  // generate the trial vector
  equal(&trial, &(lattice[r]));
  phi=phimax*(1.0-2.0*myrand());
  i=(int)(NCOMP * myrand());
  j=i+1+(int)((NCOMP-1) * myrand()); // this is to ensure j!=i
  j=j % NCOMP; 
  rotate2(&trial, i, j, phi);  
  
  // deltaE = E(trial)-E(initial)
  deltaE=-scalprod(&trial, &sumnn);
  deltaE-=-scalprod(&(lattice[r]), &sumnn);

  // metropolis step
  if(deltaE<0)
    {
    equal(&(lattice[r]), &trial);
    acc=1;
    }
  else
    {
    if(myrand()<exp(-deltaE*beta))
      {
      equal(&(lattice[r]), &trial);
      acc=1;
      }
    }

  return acc;
  }


// microcanonic update at site r
void microcan(NVec * restrict lattice, 
              long int r, 
              long int const * const restrict nnp, 
              long int const * const restrict nnm, 
              long int volume) 
  {
  int i;
  NVec sumnn, trial;
  double norma;

  zeros(&sumnn);
  for(i=0; i<DIM; i++)
     {
     plusequal(&sumnn, &(lattice[nnp[dirgeo(r, i, volume)]]) );
     plusequal(&sumnn, &(lattice[nnm[dirgeo(r, i, volume)]]) );
     }
  norma=norm(&sumnn);

  if(norma>1.0e-15)
    {
    // generate the trial vector
    // trial = 2(sumnn,lattice)*sumnn/norm(sumnn)^2-lattice
     
    timesequal(&sumnn, 1.0/norma);
    equal(&trial, &sumnn);
    timesequal(&trial, 2.0*scalprod(&(lattice[r]), &sumnn) );
    minusequal(&trial, &(lattice[r]) );

    equal(&(lattice[r]), &trial);
    }
  }


// main
int main(int argc, char **argv)
    {
    NVec *lattice;
    int i, L;
    long int r, volume, sample, iter, acc, metro; 
    long int *nnp, *nnm;
    double beta, locE, locM;
    char datafile[STRING_LENGTH];
    FILE *fp;

    const int microupdates=5;
    const double phimax=3;
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

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
      fprintf(stdout, "  number of components = %d\n\n", NCOMP);
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
    lattice=(NVec *)malloc((unsigned long int)(volume)*sizeof(NVec));
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
       one(&(lattice[r]));
       }

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    acc=0;
    metro=0;
    for(iter=0; iter<sample; iter++)
       {
       if(myrand()<0.4)
         {
         // metropolis
         for(r=0; r<volume; r++)
            {
            acc+=metropolis(lattice, r, nnp, nnm, volume, beta, phimax);
            }
         metro+=1;
         }
       else
         {
         // microcanonical updates
         for(i=0; i<microupdates; i++) 
            {
            for(r=0; r<volume; r++)
               {
               microcan(lattice, r, nnp, nnm, volume);
               }
            }
         }

       // normalize the lattice
       for(r=0; r<volume; r++)
          {
          normalize(&(lattice[r]));
          }

       locE=energy(lattice, nnp, volume);
       locM=magn(lattice, volume);

       fprintf(fp, "%.12f %.12f\n", locE, locM);
       }

    // close datafile
    fclose(fp);

    printf("Acceptance rate %f\n", (double)acc / (double)metro / (double) volume);

    free(lattice);
    free(nnp);
    free(nnm);

    return EXIT_SUCCESS;
    }


