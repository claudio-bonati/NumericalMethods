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
 
  part[0].nnm=&(part[Nt-1]);
  for(r=1; r<Nt; r++)
     {
     part[r].nnm=&(part[r-1]);
     }
  }


// Heatbath update
void heatbath(Node * restrict node, long int Nt, double eta)
  {
  double nnsum, std, average; 
  long int r; 

  for(r=0; r<Nt; r++)
     {
     nnsum = (node[r].nnp)->value + (node[r].nnm)->value;
     std=1.0/sqrt(eta+2.0/eta);
     average=nnsum/eta/(eta+2.0/eta);
 
     (node[r].value)=average+std*gauss1();
     }
  }


// overrelaxation update
void overrelaxation(Node * restrict node, long int Nt, double eta)
  {
  double nnsum, average, ris;
  long int r; 

  for(r=0; r<Nt; r++)
     {
     nnsum = (node[r].nnp)->value + (node[r].nnm)->value;
     average=nnsum/eta/(eta+2.0/eta);
     ris=2.0*average-(node[r].value);

     node[r].value=ris;
     }
  }


// random swap: return 1 if accepted 
int rand_swap(Node * restrict part1,
              Node * restrict part2, 
              long int Nt,
              double eta,
              int * restrict twisted)
  {
  Node *nodetmp;
  double Enew, Eold, aux;
  long int r;

  r=(long int)((double) Nt*myrand());

  aux=part1[r].value-(part1[r].nnp)->value;
  Eold=aux*aux/2.0/eta;
  aux=part2[r].value-(part2[r].nnp)->value;
  Eold+=aux*aux/2.0/eta;

  aux=part2[r].value-(part1[r].nnp)->value;
  Enew=aux*aux/2.0/eta;
  aux=part1[r].value-(part2[r].nnp)->value;
  Enew+=aux*aux/2.0/eta;
 
  if(Enew<Eold)
    {
    nodetmp=part1[r].nnp;
    part1[r].nnp=part2[r].nnp;
    (part1[r].nnp)->nnm=&(part1[r]);
    part2[r].nnp=nodetmp;
    (part2[r].nnp)->nnm=&(part2[r]);

    *twisted=(*twisted+1)%2;

    return 1;
    }
  else if(myrand()<exp(-(Enew-Eold)) )
         {
         nodetmp=part1[r].nnp;
         part1[r].nnp=part2[r].nnp;
         (part1[r].nnp)->nnm=&(part1[r]);
         part2[r].nnp=nodetmp;
         (part2[r].nnp)->nnm=&(part2[r]);
     
         *twisted=(*twisted+1)%2;

         return 1;
         }

  return 0; 
  }


// optimized swap: return 1 if accepted 
int opt_rand_swap(Node * restrict part1,
                  Node * restrict part2, 
                  long int Nt,
                  double eta,
                  int * restrict twisted)
  {
  Node *nodetmp;
  double Enew, Eold, aux, min;
  long int r, r1;

  // r1 is the value at which the two trajectories are closer to each other
  aux=part1[0].value-part2[0].value;
  min=aux*aux;
  r=0; 
  for(r1=1; r1<Nt; r1++)
     {
     aux=part1[r1].value-part2[r1].value;
     if(aux*aux<min)
       {
       min=aux*aux;
       r=r1;
       }
     }
     
  aux=part1[r].value-(part1[r].nnp)->value;
  Eold=aux*aux/2.0/eta;
  aux=part2[r].value-(part2[r].nnp)->value;
  Eold+=aux*aux/2.0/eta;

  aux=part2[r].value-(part1[r].nnp)->value;
  Enew=aux*aux/2.0/eta;
  aux=part1[r].value-(part2[r].nnp)->value;
  Enew+=aux*aux/2.0/eta;
 
  if(Enew<Eold)
    {
    nodetmp=part1[r].nnp;
    part1[r].nnp=part2[r].nnp;
    (part1[r].nnp)->nnm=&(part1[r]);
    part2[r].nnp=nodetmp;
    (part2[r].nnp)->nnm=&(part2[r]);

    *twisted=(*twisted+1)%2;

    return 1;
    }
  else if(myrand()<exp(-(Enew-Eold)) )
         {
         nodetmp=part1[r].nnp;
         part1[r].nnp=part2[r].nnp;
         (part1[r].nnp)->nnm=&(part1[r]);
         part2[r].nnp=nodetmp;
         (part2[r].nnp)->nnm=&(part2[r]);
     
         *twisted=(*twisted+1)%2;

         return 1;
         }

  return 0; 
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
  

// main
int main(int argc, char **argv)
    {
    Node *part1, *part2;
    long int Nt, sample, iter, accswap; 
    double simbeta, eta;
    double x1, x2, x1_2, x2_2;
    int j, twisted=0; // if =0 the two configuration are not twisted together
  
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
      fprintf(stdout, "  x1, x1^2, x2, x2^2, twisted for every configuration, where\n");
      fprintf(stdout, "  twisted=0 if the two trajectories are independent, =1 if a cycle is present.\n\n");

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

    accswap=0;
    for(iter=0; iter<sample; iter++)
       {
       if(myrand()<0.5)
         {
         // heatbath
         heatbath(part1, Nt, eta);
         heatbath(part2, Nt, eta);
         }
       else
         {
         // overrelaxation
         for(j=0; j<overrelaxsteps; j++)
            {
            overrelaxation(part1, Nt, eta);
            overrelaxation(part2, Nt, eta);
            }
         }

       // swap
       accswap+=rand_swap(part1, part2, Nt, eta, &twisted);

       if(iter%measevery==0)
         {
         compute_x_and_x2(part1, Nt, &x1, &x1_2);
         compute_x_and_x2(part2, Nt, &x2, &x2_2);

         fprintf(fp, "%.12f %.12f %.12f %.12f %d\n", x1, x1_2, x2, x2_2, twisted);
         }
       }
    printf("Acceptance rate for swap %f\n", (double)accswap / (double)sample);

    // close datafile
    fclose(fp);

    // free memory
    free(part1);
    free(part2);

    return EXIT_SUCCESS;
    }


