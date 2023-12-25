#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/boxmuller.h"
#include"../include/geometry_st.h"
#include"../include/random.h"

// when using C99 M_PI is not defined in math.h header!
#ifndef M_PI
#  define M_PI  3.141592653589793238462643383279502884
#endif

#define STDIM 2  // space-time dimensionality
#define STRING_LENGTH 50

void phi_at_fixed_momentum(double const * const restrict lattice,
                           double const p[STDIM-1],
                           int Nt, 
                           int Ns,
                           double complex * restrict O)
  {
  long int r, stvol;
  int i, coord[STDIM];
  double sum;

  for(i=0; i<Nt; i++)
     {
     O[i]=0.0+I*0.0;
     }

  stvol=Nt;
  for(i=1; i<STDIM; i++)
     {
     stvol*=Ns;
     }

  for(r=0; r<stvol; r++)
     {
     lex_to_cart_st(coord, r, Nt, Ns, STDIM);

     sum=0.0;
     for(i=1; i<STDIM; i++)
        {
        sum+=p[i-1]*coord[i];
        }
     
     O[coord[0]]+=lattice[r]*exp(I*sum);
     } 

  for(i=0; i<Nt; i++)
     {
     O[i]/=sqrt((double)stvol);
     }
  }


void time_corr_meas(double const * const restrict lattice,
                    int Nt, 
                    int Ns,
                    double *corr_p0)
   {
   int i, j, k;
   double p[STDIM-1];
   double complex *O;  

   for(i=0; i<Nt/4; i++)
      {
      corr_p0[i]=0.0;
      }

   // allocate O
   O=(double complex *)malloc((unsigned long int)(Nt)*sizeof(double complex));
   if(O == NULL)
     {
     fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   // zero momentum
   for(i=0; i<STDIM-1; i++)
      {
      p[i]=0.0;
      }

   phi_at_fixed_momentum(lattice, p, Nt, Ns, O);
   for(i=0; i<Nt/4; i++)
      {
      for(j=0; j<Nt; j++)
         {
         k=(j+i)%Nt;
    
         corr_p0[i]+=creal(conj(O[j])*O[k]);
         }
      corr_p0[i]/=Nt;
      }

   free(O);
   }


// heatbath update at site r
void heatbath(double * restrict lattice, 
              long int r, 
              long int const * const restrict nnp,  
              long int const * const restrict nnm,
              long int stvol,
              double hatm)
   {
   int i;
   double nnsum, mean, std;

   nnsum=0.0;
   for(i=0; i<STDIM; i++)
      {
      nnsum+=lattice[nnp[dirgeo(r, i, stvol)]];
      nnsum+=lattice[nnm[dirgeo(r, i, stvol)]];
      }
   
   mean=nnsum/(hatm*hatm+2.0*(double)STDIM);
   std=1.0/sqrt(hatm*hatm+2.0*(double)STDIM);
   
   lattice[r]=mean+std*gauss1();
   }


// microcanonical update
void overrelaxation(double * restrict lattice, 
                    long int r, 
                    long int const * const restrict nnp,  
                    long int const * const restrict nnm, 
                    long int stvol, 
                    double hatm)
   {
   int i;
   double nnsum, mean, old;

   nnsum=0.0;
   for(i=0; i<STDIM; i++)
      {
      nnsum+=lattice[nnp[dirgeo(r, i, stvol)]];
      nnsum+=lattice[nnm[dirgeo(r, i, stvol)]];
      }
   
   mean=nnsum/(hatm*hatm+2.0*(double)STDIM);
   old=lattice[r];
   
   lattice[r]=2*mean-old;
   }


// main
int main(int argc, char **argv)
   {
   int i, Nt, Ns;
   double hatm, rand, *lattice, *corr_p0;
   long int *nnp, *nnm;
   long int iter, sample, r, stvolume;

   char datafile[STRING_LENGTH];
   FILE *fp;

   const int overrelax=5;
   const int measevery=1;
   
   const unsigned long int seed1=(unsigned long int) time(NULL);
   const unsigned long int seed2=seed1+127;

   if(argc != 6)
     {
     fprintf(stdout, "How to use this program:\n");
     fprintf(stdout, "  %s Nt Ns hatm sample datafile\n\n", argv[0]);
     fprintf(stdout, "  Nt = temporal size of the lattice\n");
     fprintf(stdout, "  Ns = spatial size of the lattice (space-time dimension defined by macro STDIM)\n");
     fprintf(stdout, "  hatm = a*m\n");
     fprintf(stdout, "  sample = number of drawn to be extracted\n");
     fprintf(stdout, "  datafile = name of the file on which to write the data\n\n");
     fprintf(stdout, "Compiled for:\n");
     fprintf(stdout, "  dimensionality = %d\n\n", STDIM);
     fprintf(stdout, "Output:\n");
     fprintf(stdout, "  each line is made up of corr_0[i] for 0<=i<Nt/4 where corr_0[i] are the\n"); 
     fprintf(stdout, "  correlators of the scalar field computed at zero momentum\n");

     return EXIT_SUCCESS;
     }
   else
     {  
     // read input values 
     Nt=atoi(argv[1]);
     Ns=atoi(argv[2]);
     hatm=atof(argv[3]);
     sample=atol(argv[4]);

     if(strlen(argv[5]) >= STRING_LENGTH)
       {
       fprintf(stderr, "File name too long. Increse STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
       return EXIT_FAILURE;
       }
     else
       {
       strcpy(datafile, argv[5]);
       }
     }

   // initialize random number generator
   myrand_init(seed1, seed2);

   // compute the spacetime volume
   stvolume=Nt;
   for(i=1; i<STDIM; i++)
      {
      stvolume*=Ns;
      }

   // allocate the lattice (lexicographic order)
   // and next neighbors: nnp[dirgeo(r, i, volume)]= next neighbor in positive "i" direction of site r 
   lattice=(double *)malloc((unsigned long int)(stvolume)*sizeof(double));
   if(lattice == NULL)
     {
     fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
     return EXIT_FAILURE;
     }
   nnp=(long int *)malloc((unsigned long int)(STDIM*stvolume)*sizeof(long int));
   if(nnp == NULL){
     fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
     return EXIT_FAILURE;
     }
   nnm=(long int *)malloc((unsigned long int)(STDIM*stvolume)*sizeof(long int));
   if(nnm == NULL){
     fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
     return EXIT_FAILURE;
     }

   // allocate correlators
   corr_p0=(double *)malloc((unsigned long int)(Nt/4)*sizeof(double));
   if(corr_p0 == NULL)
     {
     fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
     return EXIT_FAILURE;
     }

   // initialize nnp and nnm
   init_neighbors_st(nnp, nnm, Nt, Ns, STDIM);

   // initialize lattice to ordered start
   for(r=0; r<stvolume; r++)
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

   for(iter=0; iter<sample; iter++)
      {
      rand=myrand(); 

      if(rand<1.0/(double)overrelax)
        {
        // heatbath       
        for(r=0; r<stvolume; r++)
           {
           heatbath(lattice, r, nnp, nnm, stvolume, hatm);
           }
         }
      else
        {
        // overrelaxation
        for(r=0; r<stvolume; r++)
           {
           overrelaxation(lattice, r, nnp, nnm, stvolume, hatm);
           }
         }
      
      if(iter%measevery==0)
        {
        // perform measures
        time_corr_meas(lattice, Nt, Ns, corr_p0);

        for(i=0; i<Nt/4; i++)
           {
           fprintf(fp, "%.12f ", corr_p0[i]);
           }
        fprintf(fp, "\n");
        }
      }

   // close datafile
   fclose(fp);

   // free memory
   free(lattice);
   free(nnp);
   free(nnm);
   free(corr_p0);


   return EXIT_SUCCESS;
   }


