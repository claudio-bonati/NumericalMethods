#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/random.h"

#define SIZE 4
#define STRING_LENGTH 50

// magnetization per site
double magn(int lattice[SIZE][SIZE])
  {
  long int rx, ry, sum;

  sum=0;
  for(rx=0; rx<SIZE; rx++)
     {
     for(ry=0; ry<SIZE; ry++)
        {
        sum+=lattice[rx][ry];
        }
     }

  return (double) sum / (double) (SIZE*SIZE);
  }


// energy per site
double energy(int lattice[SIZE][SIZE]) 
  {
  long int rx, ry, tmp, sum;

  sum=0;
  for(rx=0; rx<SIZE; rx++)
     {
     for(ry=0; ry<SIZE; ry++)
        {
        tmp=(rx+1)%SIZE;
        sum+=-lattice[rx][ry]*lattice[tmp][ry];
       
        tmp=(ry+1)%SIZE;
        sum+=-lattice[rx][ry]*lattice[rx][tmp];
        }
     }

  return (double) sum / (double) (SIZE*SIZE);
  }


// metropolis update at site r
// return 1 if accepted, else 0
int metropolis(int lattice[SIZE][SIZE], 
               long int rx,
               long int ry,
               double const * const restrict acc_prob)
  {
  long int tmp;
  int sumnn, acc=0;

  // the relevant part of the energy will be -lattice[r]*sumnn
  sumnn=0;
  
  tmp=(rx+1)%SIZE;
  sumnn+=lattice[tmp][ry];

  tmp=(rx-1);
  if(tmp<0)
    {
    tmp=SIZE-1;
    }
  sumnn+=lattice[tmp][ry];

  tmp=(ry+1)%SIZE;
  sumnn+=lattice[rx][tmp];

  tmp=(ry-1);
  if(tmp<0)
    {
    tmp=SIZE-1;
    }
  sumnn+=lattice[rx][tmp];

  sumnn*=lattice[rx][ry];

  // metropolis step
  if(sumnn<0)
    {
    lattice[rx][ry]=-lattice[rx][ry];
    acc=1;
    }
  else
    {
    if(myrand()<acc_prob[sumnn])  // remember that acc_prob[i]=exp(-2.0*beta*((double)i));
      {
      lattice[rx][ry]=-lattice[rx][ry];
      acc=1;
      }
    }

  return acc;
  }


// main
int main(int argc, char **argv)
    {
    int i, lattice[SIZE][SIZE];
    long int rx, ry, rxaux, ryaux, sample, iter, acc; 
    double beta, locE, locM;
    double acc_prob[5];
  
    char datafile[STRING_LENGTH];
    FILE *fp;

    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s beta sample datafile\n\n", argv[0]);
      fprintf(stdout, "  beta = inverse temperature\n");
      fprintf(stdout, "  sample = number of drawn to be extracted\n");
      fprintf(stdout, "  datafile = name of the file on which to write the data\n");
      fprintf(stdout, "  size of the lattice define by macro: now it is %d\n\n", SIZE);
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  E, M (E=energy per site, M=magnetization per site), one line for each draw\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 
      beta=atof(argv[1]);
      sample=atol(argv[2]);

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

    if(sample<=0)
      {
      fprintf(stderr, "'sample' must be positive\n");
      return EXIT_FAILURE;
      }

    // initialize random number generator
    myrand_init(seed1, seed2);

    // initialize lattice to ordered start
    for(rx=0; rx<SIZE; rx++)
       {
       for(ry=0; ry<SIZE; ry++)
          {
          lattice[rx][ry]=1;
          }
       }

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize acceptance probability
    for(i=0; i<5; i++)
       {
       acc_prob[i]=exp(-2.0*beta*((double)i));
       }

    acc=0;
    for(iter=0; iter<sample; iter++)
       {
       for(rx=0; rx<SIZE; rx++)
          {
          for(ry=0; ry<SIZE; ry++)
             {
             rxaux=(long int)((double)SIZE * myrand());
             ryaux=(long int)((double)SIZE * myrand());

             acc+=metropolis(lattice, rxaux, ryaux, acc_prob);
             }
          }

       locE=energy(lattice);
       locM=magn(lattice);

       fprintf(fp, "%.12f %.12f\n", locE, locM);
       }

    // close datafile
    fclose(fp);

    printf("Acceptance rate %f\n", (double)acc / (double)sample / (double) (SIZE*SIZE));

    return EXIT_SUCCESS;
    }


