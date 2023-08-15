#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/boxmuller.h"
#include"../include/random.h"

#define STRING_LENGTH 50

//#define METRO // if defined Metropolis update is used

// stucture for the configuration 
typedef struct Node {
  double value;  

  struct Node *nnp; // point to the next node of the configuration
  struct Node *nnm; // point to the previous node of the configuration

  } Node;


// intialize the configuration with the value x0
void init_conf(Node * restrict part, long int Nt, double x0)
  {
  long int r;

  for(r=0; r<Nt; r++)
     {
     part[r].value=x0;
     }

  for(r=0; r<Nt-1; r++)
     {
     part[r].nnp=&(part[r+1]);
     }
  part[Nt-1].nnp=&(part[0]);
 
  for(r=1; r<Nt; r++)
     {
     part[r].nnm=&(part[r-1]);
     }
  part[0].nnm=&(part[Nt-1]);
  }


// Metropolis update, return 1 if accepted
int metropolis(Node * restrict node, double eta)
  {
  const double delta=10.0*sqrt(eta);
  double trial, Eold, Enew, nsum;

  double nnsum = (node->nnp)->value + (node->nnm)->value;

  Eold=(node->value)*(node->value)*(eta/2.0+1./eta)-(node->value)*nnsum/eta;
  trial=node->value + delta*(1.0-2.0*myrand());
  Enew=trial*trial*(eta/2.0+1./eta)-trial*nnsum/eta;

  if(Enew<Eold)
    {
    node->value=trial;
    return 1;
    }
  else if(myrand()<exp(-(Enew-Eold)) )
         {
         node->value=trial;
         return 1;
         }

  return 0; 
  }


// Heatbath update, return 1
int heatbath(Node * restrict node, double eta)
  {
  double nnsum = (node->nnp)->value + (node->nnm)->value;

  const double std=1.0/sqrt(eta+2.0/eta);
  const double average=nnsum/eta/(eta+2.0/eta);
 
  (node->value)=average+std*gauss1();

  return 1; 
  }


// compute x and x2
void compute_x_and_x2(Node const * const restrict part, 
                      long int Nt, 
                      double * restrict x,  
                      double * restrict x2)
  {
  long int r;
  double aux;

  *x=0.0;
  *x2=0.0;

  for(r=0; r<Nt; r++)
     {
     aux=part[r].value;
   
     *x+=aux;
     *x2+=aux*aux;
     }
  
  *x/=(double) Nt;
  *x2/=(double) Nt;
  }
  

// compute x12 and x12_2
void compute_x12_and_x12_2(Node const * const restrict part1, 
                           Node const * const restrict part2,
                           long int Nt, 
                           double * restrict x12,  
                           double * restrict x12_2)
  {
  long int r;
  double aux;

  *x12=0.0;
  *x12_2=0.0;

  for(r=0; r<Nt; r++)
     {
     aux=part1[r].value-part2[r].value;
   
     *x12+=aux;
     *x12_2+=aux*aux;
     }
  
  *x12/=(double) Nt;
  *x12_2/=(double) Nt;
  }


// main
int main(int argc, char **argv)
    {
    Node *part1, *part2;
    long int Nt, r, sample, iter, acc, accswap; 
    double simbeta, eta;
    double x1, x2, x12, x1_2, x2_2, x12_2;
    int twisted=0; // if =0 the two configuration are not twisted together
  
    char datafile[STRING_LENGTH];
    FILE *fp;

    const int measevery=10;
    const int swapevery=10;

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
      fprintf(stdout, "  TO BE MODIFIED\n");

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

    // allocate the configuration
    part1=(Node *)malloc((unsigned long int)(Nt)*sizeof(Node));
    if(part1 == NULL){
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    part2=(Node *)malloc((unsigned long int)(Nt)*sizeof(Node));
    if(part1 == NULL){
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // intialize the configuration
    init_conf(part1, Nt, 0.1);
    init_conf(part2, Nt, 0.2);

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    acc=0;
    accswap=0;
    for(iter=0; iter<sample; iter++)
       {
       #ifndef METRO
         for(r=0; r<Nt; r++)
            { 
            acc+=heatbath(&(part1[r]), eta);
            }  
         for(r=0; r<Nt; r++)
            { 
            acc+=heatbath(&(part2[r]), eta);
            }  
       #else
         for(r=0; r<Nt; r++)
            { 
            acc+=metropolis(&(part1[r]), eta);
            }  
         for(r=0; r<Nt; r++)
            { 
            acc+=metropolis(&(part2[r]), eta);
            }  
       #endif

       if(iter%measevery==0)
         {
         compute_x_and_x2(part1, Nt, &x1, &x1_2);
         compute_x_and_x2(part2, Nt, &x2, &x2_2);
         compute_x12_and_x12_2(part1, part2, Nt, &x12, &x12_2);

         fprintf(fp, "%f %f %f %f %f %f %d\n", x1, x1_2, x2, x2_2, x12, x12_2, twisted);
         }
       }
    printf("Acceptance rate %f\n", (double)acc / (double)sample / (double) Nt /2.0);

    // close datafile
    fclose(fp);

    // free memory
    free(part1);
    free(part2);

    return EXIT_SUCCESS;
    }


