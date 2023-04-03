#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/geometry.h"
#include"../include/random.h"

#define DIM 2
#define STRING_LENGTH 50

// magnetization per site
double magn(int *lattice, long int volume)
  {
  long int r, sum;

  sum=0;
  for(r=0; r<volume; r++)
     {
     sum+=lattice[r];
     }

  return (double) sum / (double) volume;
  }


// energy per site
double energy(int *lattice, long int *nnp, long int volume)
  {
  long int r, sum;
  int i;

  sum=0;
  for(r=0; r<volume; r++)
     {
     for(i=0; i<DIM; i++)
        {
        sum+=-lattice[r]*lattice[nnp[i*volume + r]];
        }
     }

  return (double) sum / (double) volume;
  }


// main
int main(int argc, char **argv)
    {
    int i, L, *lattice, sumnn;
    long int r, volume, sample, iter, acc; 
    long int *nnp, *nnm;
    double beta, locE, locM;
    double acc_prob[2*DIM+1];
  
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
    // and next neighbors: nnp[i*volume+r]= next neighbor in positive "i" direction of site r 
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
       lattice[r]=1;
       }

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize acceptance probability
    for(i=0; i<2*DIM+1; i++)
       {
       acc_prob[i]=exp(-2.0*beta*((double)i));
       }

    acc=0;
    for(iter=0; iter<sample; iter++)
       {
       for(r=0; r<volume; r++)
          {
          // the relevant part of the energy will be -lattice[r]*sum
          sumnn=0;
          for(i=0; i<DIM; i++)
             {
             sumnn+=lattice[nnp[i*volume+r]];
             sumnn+=lattice[nnm[i*volume+r]];
             }
          sumnn*=lattice[r];

          // metropolis step
          if(sumnn<0)
            {
            lattice[r]=-lattice[r];
            acc++;
            }
          else
            {
            if(myrand()<acc_prob[sumnn])
              {
              lattice[r]=-lattice[r];
              acc++;
              }
            }
          }

       locE=energy(lattice, nnp, volume);
       locM=magn(lattice, volume);

       fprintf(fp, "%f %f\n", locE, locM);
       }

    // close datafile
    fclose(fp);

    printf("Acceptance rate %f\n", (double)acc / (double)sample / (double) volume);

    free(lattice);
    free(nnp);
    free(nnm);

    return EXIT_SUCCESS;
    }


