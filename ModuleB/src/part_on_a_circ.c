#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/random.h"

#define STRING_LENGTH 50


// stucture for the configuration 
typedef struct Conf {
  double *lattice;  

  double simbeta;
  long int Nt;
  double eta;
  } Conf;


// initialize the configuration
void init_conf(Conf *config, double simbeta, long int Nt)
  {
  long int r;

  config->lattice=(double *)malloc((unsigned long int)(Nt)*sizeof(double));
  if(config->lattice == NULL)
    {
    fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  for(r=0; r<Nt; r++)
     {
     config->lattice[r]=0.5;
     }

  config->simbeta=simbeta;
  config->Nt=Nt;
  config->eta=simbeta/(double) Nt;
  }


// free the configuration
void free_conf(Conf *config)
  {
  free(config->lattice);
  }


// oriented distance between x, y in [0, 1)
double dist(double x, double y)
  {
  if(fabs(x-y)<= 0.5)
    {
    return x-y;
    }
  else if(x-y>0.5)
         {
         return x-y-1.0;
         }
  else 
    {
    return x-y+1.0;
    }
  }


// Metropolis update, return the number of accepted steps
long int metropolis(Conf * restrict config, 
                    long int const * const restrict nnp, 
                    long int const * const restrict nnm)
  {
  const double delta=0.5; //10.0*sqrt(config->eta);
  double trial, Eold, Enew, aux;
  long int r, acc;

  acc=0; 

  for(r=0; r<config->Nt; r++)
     {
     aux=dist(config->lattice[nnp[r]], config->lattice[r]);
     Eold=1.0/(2.0*config->eta)*aux*aux;
     aux=dist(config->lattice[r], config->lattice[nnm[r]]);
     Eold+=1.0/(2.0*config->eta)*aux*aux;
   
     trial=config->lattice[r] + delta*(1.0-2.0*myrand());

     while(trial >=1.0)
       {
       trial-=1.0;
       }
     while(trial<0.0)
       {
       trial+=1.0;
       }

     aux=dist(config->lattice[nnp[r]], trial);
     Enew=1.0/(2.0*config->eta)*aux*aux;
     aux=dist(trial, config->lattice[nnm[r]]);
     Enew+=1.0/(2.0*config->eta)*aux*aux;

     if(Enew<Eold)
       {
       config->lattice[r]=trial;
       acc+=1;
       }
     else if(myrand()<exp(-(Enew-Eold)) )
            {
            config->lattice[r]=trial;
            acc+=1;
            }
     } 

  return acc;
  }


// winding number
long int calc_Q(Conf const * const restrict config,
                long int const * const restrict nnp)
  {
  double risf;
  long int r, Q; 

  risf=0.0;
  for(r=0; r<config->Nt; r++)
     {
     risf+=dist(config->lattice[nnp[r]], config->lattice[r]);
     }
 
  Q=lround(risf); 

  return Q;
  }



// main
int main(int argc, char **argv)
    {
    Conf config;
    long int Nt, r, sample, iter, acc; 
    long int *nnp, *nnm;
    double simbeta;
    long int Q; 
  
    char datafile[STRING_LENGTH];
    FILE *fp;

    const int measevery=10;

    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    if(argc != 5)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s simbeta Nt sample datafile\n\n", argv[0]);
      fprintf(stdout, "  simbeta = hbar^2/(4*pi^2*m*R^2*k_B*T)\n");
      fprintf(stdout, "  Nt = number of temporal steps\n");
      fprintf(stdout, "  sample = number of drawn to be extracted\n");
      fprintf(stdout, "  datafile = name of the file on which to write the data\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  winding number, one line for each configuration\n\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 
      simbeta=atof(argv[1]);
      Nt=atol(argv[2]);
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

    if(Nt<=1)
      {
      fprintf(stderr, "'Nt' must be larger than 1\n");
      return EXIT_FAILURE;
      }

    if(sample<=0)
      {
      fprintf(stderr, "'sample' must be positive\n");
      return EXIT_FAILURE;
      }

    // initialize random number generator
    myrand_init(seed1, seed2);

    // allocate next neighbors
    nnp=(long int *)malloc((unsigned long int)(Nt)*sizeof(long int));
    if(nnp == NULL){
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    nnm=(long int *)malloc((unsigned long int)(Nt)*sizeof(long int));
    if(nnm == NULL){
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize nnp and nnm for periodic b.c.
    for(r=0; r<Nt; r++)
       {
       nnp[r]=r+1;
       nnm[r]=r-1;
       }
    nnp[Nt-1]=0;
    nnm[0]=Nt-1;

    // allocate and initialize the configuration
    init_conf(&config, simbeta, Nt); 

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    acc=0;
    for(iter=0; iter<sample; iter++)
       {
       // complete update of the configuration
       acc+=metropolis(&config, nnp, nnm);

       if(iter%measevery==0)
         {
         Q=calc_Q(&config, nnp);
 
         fprintf(fp, "%ld\n", Q);
         }
       }
    printf("Acceptance rate %f\n", (double)acc / (double)sample / (double) Nt);

    // close datafile
    fclose(fp);

    // free memory
    free_conf(&config);
    free(nnp);
    free(nnm);

    return EXIT_SUCCESS;
    }


