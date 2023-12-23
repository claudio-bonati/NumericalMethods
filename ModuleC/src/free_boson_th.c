#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/boxmuller.h"
#include"../include/geometry_st.h"
#include"../include/random.h"

#define STDIM 2  // space-time dimensionality
#define STRING_LENGTH 50

// O1 = (1/stvol) sum_r (hatm^2 phi_r^2)
double obsO1(double const * const restrict lattice, 
             long int stvol, 
             double hatm)
   {
   long int r;
   double ris;

   ris=0.0;
   for(r=0; r<stvol; r++)
      {
      ris+=hatm*hatm*lattice[r]*lattice[r];
      }
   ris/=(double) stvol;

   return ris;
   }


// O2 = (1/stvol) sum_r sum_{mu>0} (phi_(r+mu)-phi_r)^2
double obsO2(double const * const restrict lattice, 
             long int const * const restrict nnp, 
             long int stvol) 
   {
   int i;
   long int r;
   double ris, aux;

   ris=0.0;
   for(r=0; r<stvol; r++)
      {
      for(i=1; i<STDIM; i++)
         {
         aux=lattice[nnp[dirgeo(r, i, stvol)]]-lattice[r];
         ris+=aux*aux;
         }
      }
   ris/=(double) stvol;

   return ris;
   }


// O3 = (1/stvol) sum_r (phi_(r+0)-phi_r)^2
double obsO3(double const * const restrict lattice, 
             long int const * const restrict nnp, 
             long int stvol) 
   {
   long int r;
   double ris, aux;

   ris=0.0;
   for(r=0; r<stvol; r++)
      {
      aux=lattice[nnp[dirgeo(r, 0, stvol)]]-lattice[r];
      ris+=aux*aux;
      }
   ris/=(double) stvol;

   return ris;
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
   double hatm, rand, *lattice, O1, O2, O3;
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
     fprintf(stdout, "  O1 O2 O3\n");
     fprintf(stdout, "  O1 = 1/(Nt*Ns^{STDIM-1}) sum_r(hatm*phi_r)^2\n");
     fprintf(stdout, "  O2 = 1/(Nt*Ns^{STDIM-1}) sum_r sum_{mu>0} (phi_{r+mu}-phi_r)^2\n");
     fprintf(stdout, "  O3 = 1/(Nt*Ns^{STDIM-1}) sum_r (phi_{r+0}-phi_r)^2\n");

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
        O1=obsO1(lattice, stvolume, hatm);
        O2=obsO2(lattice, nnp, stvolume);
        O3=obsO3(lattice, nnp, stvolume);

        fprintf(fp, "%.12f %.12f %.12f\n", O1, O2, O3);
        }
      }

   // close datafile
   fclose(fp);

   // free memory
   free(lattice);
   free(nnp);
   free(nnm);

   return EXIT_SUCCESS;
   }


