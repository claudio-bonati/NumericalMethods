#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/boxmuller.h"
#include"../include/random.h"

#define STRING_LENGTH 50

//#define METRO // if this macro is defined Metropolis is uded, otherwise heatbath

// average position square
double calc_x2(double const * const restrict lattice, long int Nt)
  {
  long int r;
  double ris;

  ris=0.0;
  for(r=0; r<Nt; r++)
     {
     ris+=lattice[r]*lattice[r];
     }
  
  return ris/(double)Nt;
  }

// correlators at temporal distance deltat
//  ris[0]: x x correlator
//  ris[1]: x^2 x^2 correlator
//  ris[2]: x^3 x^3 correlator
//  ris[3]: x^3-(3/2)x x^3-(3/2)x correlator  
void correlators(double const * const restrict lattice, 
                 long Nt,
                 long deltat,
                 double ris[4])
  {
  long int r, raux;
  int i;

  for(i=0; i<4; i++)
     {
     ris[i]=0.0;
     }

  for(r=0; r<Nt; r++)
     {
     raux=(r+deltat)%Nt;

     ris[0]+=lattice[r]*lattice[raux];
     ris[1]+=pow(lattice[r],2)*pow(lattice[raux],2);
     ris[2]+=pow(lattice[r],3)*pow(lattice[raux],3);
     ris[3]+=(pow(lattice[r],3)-1.5*lattice[r])*(pow(lattice[raux],3)-1.5*lattice[raux]);
     }
 
  for(i=0; i<4; i++)
     {
     ris[i]/=(double)Nt;     
     }
  }


// Metropolis update, return 1 if accepted
int metropolis(double * restrict lattice, long int r, double nnsum, double eta)
  {
  const double delta=10.0*sqrt(eta);
  double trial, Eold, Enew;

  Eold=lattice[r]*lattice[r]*(eta/2.0+1./eta)-lattice[r]*nnsum/eta;
  trial=lattice[r]+delta*(1.0-2.0*myrand());
  Enew=trial*trial*(eta/2.0+1./eta)-trial*nnsum/eta;

  if(Enew<Eold)
    {
    lattice[r]=trial;
    return 1;
    }
  else if(myrand()<exp(-(Enew-Eold)) )
         {
         lattice[r]=trial;
         return 1;
         }

  return 0; 
  }


// Heatbath update, return 1
int heatbath(double * restrict lattice, long int r, double nnsum, double eta)
  {
  const double std=1.0/sqrt(eta+2.0/eta);
  const double average=nnsum/eta/(eta+2.0/eta);
 
  lattice[r]=average+std*gauss1();

  return 1; 
  }


// overrelaxation
void overrelaxation(double * restrict lattice, long int r, double nnsum, double eta)
  {
  const double average=nnsum/eta/(eta+2.0/eta);
  double ris=2.0*average-lattice[r];
  lattice[r]=ris;
  }


// main
int main(int argc, char **argv)
    {
    double *lattice;
    long int Nt, r, sample, iter, acc, deltat; 
    long int *nnp, *nnm;
    double simbeta, eta, nnsum;
    double x2, corr[4];  // 4 different correlators are measured, see the function "correlators"
    int j;
  
    char datafile[STRING_LENGTH];
    FILE *fp;

    const int measevery=10;
    const int overrelaxsteps=5;

    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    if(argc != 5)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s simbeta Nt sample datafile\n\n", argv[0]);
      fprintf(stdout, "  simbeta = hbar*omega/(k_B T)\n");
      fprintf(stdout, "  Nt = number of temporal steps\n");
      fprintf(stdout, "  sample = number of drawn to be extracted\n");
      fprintf(stdout, "  datafile = name of the file on which to write the data\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  Each line starts with <x^2(0)> and then there is a sequence of\n");
      fprintf(stdout, "  <x(0)x(t)>, <x^2(0)x^2(t)>, <x^3(0)x^3(t)>, <A(0)A(t)> with A=x^3-(3/2)x\n");
      fprintf(stdout, "  starting from t=0 up to Nt/4. Different lines correspond to different configurations\n\n");

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

    // definition of eta
    eta=simbeta/(double)Nt;

    // initialize random number generator
    myrand_init(seed1, seed2);

    // allocate the lattice and next neighbors
    lattice=(double *)malloc((unsigned long int)(Nt)*sizeof(double));
    if(lattice == NULL)
      {
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
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

    // initialize lattice 
    for(r=0; r<Nt; r++)
       {
       lattice[r]=0.0;
       }

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
       if(myrand()<0.5)
         {
         for(r=0; r<Nt; r++)
            {
            nnsum=lattice[nnp[r]]+lattice[nnm[r]];

            #ifdef METRO
              acc += metropolis(lattice, r, nnsum, eta);
            #else
              acc += heatbath(lattice, r, nnsum, eta);
            #endif
            }
         }
       else
         {
         // overrelaxation
         for(j=0; j<overrelaxsteps; j++)
            {
            for(r=0; r<Nt; r++)
               {
               nnsum=lattice[nnp[r]]+lattice[nnm[r]];

               overrelaxation(lattice, r, nnsum, eta);
               }
            }
         }

       if(iter%measevery==0)
         {
         x2=calc_x2(lattice, Nt); 
         fprintf(fp, "%.12f ", x2);

         for(deltat=0; deltat<Nt/4; deltat++)
            {
            correlators(lattice, Nt, deltat, corr);
            
            for(j=0; j<4; j++)
               {
               fprintf(fp, "%.12f ", corr[j]);
               }
            }
         fprintf(fp, "\n");
         }
       }

    printf("Acceptance rate %f\n", (double)acc / (double)sample / (double) Nt);

    // close datafile
    fclose(fp);

    free(lattice);
    free(nnp);
    free(nnm);

    return EXIT_SUCCESS;
    }


