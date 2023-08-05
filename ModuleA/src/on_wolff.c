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


// non-recursive construction of the culter
void build_cluster_norec(NVec const * const restrict lattice, 
                         long int r, 
                         int * restrict occup, 
                         long int * restrict pointtoocc, 
                         long int * restrict clustersize,
                         long int const * const restrict nnp, 
                         long int const * const restrict nnm,
                         long int volume,
                         double beta,
                         NVec const * const restrict randv)
  {
  (void) r; // just to avoid warnings
  int i;
  double prob;
  long int index, r1;
  long int oldcs, oldcsnew;

  oldcs=0; // starting value, with *clustersize=1

  // if first neighbors have the same orientation and are not occupied
  // they are added to the cluster with probability prob

  while(*clustersize>oldcs) // this means that there are sites recently added, whose neighbors has not been checked yet, so we check them
       { 
       oldcsnew=*clustersize;

       for(index=oldcs; index<oldcsnew; index++)
          {
          r1=pointtoocc[index];

          for(i=0; i<DIM; i++)
             {
             // forward
             if(occup[nnp[dirgeo(r1, i, volume)]]==0)
               {
               prob=-2.0*beta*scalprod(randv, &(lattice[r1]))*scalprod(randv, &(lattice[nnp[dirgeo(r1, i, volume)]]));
               if(prob<0.0)
                 {
                 prob=1.0-exp(prob);
                 if(myrand()<prob)
                   {
                   occup[nnp[dirgeo(r1, i, volume)]]=1;
                   pointtoocc[*clustersize]=nnp[dirgeo(r1, i, volume)];
                   (*clustersize)++;
                   }
                 }
               }
        
             // backward
             if(occup[nnm[dirgeo(r1, i, volume)]]==0)
               {
               prob=-2.0*beta*scalprod(randv, &(lattice[r1]))*scalprod(randv, &(lattice[nnm[dirgeo(r1, i, volume)]]));
               if(prob<0.0)
                 {
                 prob=1.0-exp(prob);
                 if(myrand()<prob)
                   {
                   occup[nnm[dirgeo(r1, i, volume)]]=1;
                   pointtoocc[*clustersize]=nnm[dirgeo(r1, i, volume)];
                   (*clustersize)++;
                   }
                 }
               }
             } 
          }

       oldcs=oldcsnew;
       }
  }


// main
int main(int argc, char **argv)
    {
    NVec *lattice, randv, aux;
    int i, L, *occup;
    long int r, volume, sample, iter, clustersize; 
    long int *nnp, *nnm, *pointtoocc;
    double beta, locE, locM;
  
    char datafile[STRING_LENGTH];
    FILE *fp;

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
      fprintf(stdout, "  dimensionality = %d\n\n", DIM);
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

    // this structure will be used to keep trak of the occupied sites 
    // while building the cluster
    // 0 = free, 1=occupied
    occup=(int *)malloc((unsigned long int)(volume)*sizeof(int));
    if(occup== NULL)
      {
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // this structure will be used to keep trak of the occ 
    // while building the cluster
    // 0 = free, 1=occupied
    pointtoocc=(long int *)malloc((unsigned long int)(volume)*sizeof(long int));
    if(pointtoocc== NULL)
      {
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

    for(iter=0; iter<sample; iter++)
       {
       for(r=0; r<volume; r++) 
          {
          occup[r]=0;
          }
       clustersize=0;

       r=(int)((double)volume*myrand());
       occup[r]=1; // r is set as occupied
       pointtoocc[clustersize]=r; // a pointer to "r" is added in position "clustersize"
       clustersize++;

       // select random direction
       randvec(&randv);

       build_cluster_norec(lattice, r, occup, pointtoocc, &clustersize, nnp, nnm, volume, beta, &randv);

       // flip the cluster
       for(r=0; r<clustersize; r++)
          {
          equal(&aux, &randv);
          timesequal(&aux, -2.0*scalprod(&(lattice[pointtoocc[r]]), &randv) );
          plusequal(&aux, &(lattice[pointtoocc[r]]));

          equal(&(lattice[pointtoocc[r]]), &aux);
          }

       locE=energy(lattice, nnp, volume);
       locM=magn(lattice, volume);

       fprintf(fp, "%.12f %.12f\n", locE, locM);

       if(iter%10==0)
         {
         for(r=0; r<volume; r++)
            {
            normalize(&(lattice[r]));
            }
         }
       }

    // close datafile
    fclose(fp);

    free(lattice);
    free(occup);
    free(pointtoocc);
    free(nnp);
    free(nnm);

    return EXIT_SUCCESS;
    }
