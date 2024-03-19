#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/geometry_st.h"
#include"../include/random.h"

// when using C99 M_PI is not defined in math.h header!
#ifndef M_PI
  #define M_PI  3.141592653589793238462643383279502884
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))  // min of two numbers

#define STDIM 2  // space-time dimensionality
#define STRING_LENGTH 50

void staple(double complex ** restrict lattice, 
            long int const * const restrict nnp,
            long int const * const restrict nnm,
            long int r,
            int dir,
            long int stvol,
            double complex *ris)
  {
//               ^ dir
//         (4)   |   (1)
//     +----<----+---->----+
//     |         |         |
//  (5)|         |         |
//     V         ^         V (2)
//     |         |         |
//     +---->----+----<----+---> i
//     r1  (6)   r   (3)
//
  
  int i, k;
  long int r1;
  double complex aux;
  
  *ris=0.0+I*0.0;
  
  for(k=1; k<STDIM; k++)
     {
     i=(dir+k)%STDIM;
 
     // forward
     aux=lattice[nnp[dirgeo(r,dir,stvol)]][i];    // 1
     aux*=conj(lattice[nnp[dirgeo(r,i,stvol)]][dir]);  // 2
     aux*=conj(lattice[r][i]);  //3
     *ris+=aux;

     // backward
     r1=nnm[dirgeo(r,i,stvol)];
     aux=conj(lattice[nnp[dirgeo(r1,dir,stvol)]][i]);  // 4
     aux*=conj(lattice[r1][dir]);
     aux*=lattice[r1][i];
     *ris+=aux; 
     }
  }


// return 1 if the update is accepted, 0 otherwise
int metropolis(double complex ** restrict lattice, 
               long int const * const restrict nnp,
               long int const * const restrict nnm,
               long int r,
               int dir,
               double epsilon,
               double beta,
               long int stvol)
  {
  int acc=0;
  double oldS, newS, deltaS, aux;
  double complex stap, change;

  aux=epsilon*(2.0*myrand()-1.0);
  change=1.0/sqrt(1.0+aux*aux) + I*aux/sqrt(1.0+aux*aux);  // complex phase that becomes 1 if epsilon=0
  
  // compute staple
  staple(lattice, nnp, nnm, r, dir, stvol, &stap);

  oldS=-beta*creal(lattice[r][dir]*stap);
  newS=-beta*creal(lattice[r][dir]*stap*change);
  deltaS=newS-oldS;

  if(deltaS<0) 
    {
    lattice[r][dir]*=change;
    acc=1;    
    }
  else
    {
    if(myrand()<exp(-deltaS))
      {
      lattice[r][dir]*=change;
      acc=1;    
      }
    }

  return acc;
  }


void overrelaxation(double complex ** restrict lattice, 
                    long int const * const restrict nnp,
                    long int const * const restrict nnm,
                    long int r,
                    int dir,
                    long int stvol)
  {
  double complex stap, aux;

  // compute staple
  staple(lattice, nnp, nnm, r, dir, stvol, &stap);

  if(cabs(stap)>1.0e-10)
    {
    aux=conj(lattice[r][dir]*stap*stap);
    aux/=(cabs(stap)*cabs(stap));

    lattice[r][dir]=aux;
    }
  }


// multihit 
double complex multihit(double complex ** restrict lattice, 
                        long int const * const restrict nnp,
                        long int const * const restrict nnm,
                        long int r,
                        int dir,
                        double epsilon,
                        double beta,
                        int nhits,
                        long int stvol)
  {
  int i, count;
  double oldS, newS, deltaS, aux;
  double complex stap, link, change, ris=0.0+0.0*I;

  link=lattice[r][dir];

  // compute staple
  staple(lattice, nnp, nnm, r, dir, stvol, &stap);

  count=0;
  for(i=0; i<nhits; i++)
     {
     // metropolis
     aux=epsilon*(2.0*myrand()-1.0);
     change=1.0/sqrt(1.0+aux*aux) + I*aux/sqrt(1.0+aux*aux);  // complex phase that becomes 1 if epsilon=0

     oldS=-beta*creal(link*stap);
     newS=-beta*creal(link*stap*change);
     deltaS=newS-oldS;

     if(deltaS<0) 
       {
       link*=change;
       }
     else
       {
       if(myrand()<exp(-deltaS))
         {
         link*=change;
         }
       }
     ris+=link;
     count++;
     
     // overrelaxation    
     if(cabs(stap)>1.0e-10)
       {
       change=conj(link*stap*stap);
       change/=(cabs(stap)*cabs(stap));
   
       link=change;
   
       ris+=link;
       count++;
       }
     }
 
  ris/=(double complex) (count);

  return ris;
  }


double Wilsonloop(double complex ** restrict lattice, 
                  long int const * const restrict nnp,
                  long int const * const restrict nnm,
                  int Wt,
                  int Ws,
                  double epsilon,
                  double beta,
                  long int stvol)
  {
  int i, dir;
  const int nhits=10;
  long int r;
  double ris=0.0;
  double complex aux;

  for(r=0; r<stvol; r++)
     {
     for(dir=1; dir<STDIM; dir++) 
        {
//      ^ 0
//      | 
//   r1 +--------+r2       
//      |        |
//      |Wt      |
//      |        |
//      |        |
//      +--------+--> dir   
//     r    Ws   r3

        aux=1.0+0.0*I;
 
        if(Wt==1 || Ws==1)
          {
          for(i=0; i<Wt; i++)
             {
             aux*=lattice[r][0];
             r=nnp[dirgeo(r,0,stvol)];
             }
          // now we are in r1
          
          for(i=0; i<Ws; i++)
             {
             aux*=lattice[r][dir];
             r=nnp[dirgeo(r,dir,stvol)];
             }
          // now we are in r2
          
          for(i=0; i<Wt; i++)
             {
             r=nnm[dirgeo(r,0,stvol)];
             aux*=conj(lattice[r][0]);
             }
          // now we are in r3

          for(i=0; i<Ws; i++)
             {
             r=nnm[dirgeo(r,dir,stvol)];
             aux*=conj(lattice[r][dir]);
             }
          }
        else
          {
          for(i=0; i<Wt; i++)
             {
             if(i==0)
               {
               aux*=lattice[r][0];
               }
             else
               {
               aux*=multihit(lattice, nnp, nnm, r, 0, epsilon, beta, nhits, stvol);
               }
             r=nnp[dirgeo(r,0,stvol)];
             }
          // now we are in r1
          
          for(i=0; i<Ws; i++)
             {
             if(i==0)
               {
               aux*=lattice[r][dir];
               }
             else
               {
               aux*=multihit(lattice, nnp, nnm, r, dir, epsilon, beta, nhits, stvol);
               }
             r=nnp[dirgeo(r,dir,stvol)];
             }
          // now we are in r2
          
          for(i=0; i<Wt; i++)
             {
             r=nnm[dirgeo(r,0,stvol)];
             if(i==0)
               {
               aux*=conj(lattice[r][0]);
               }
             else
               {
               aux*=conj(multihit(lattice, nnp, nnm, r, 0, epsilon, beta, nhits, stvol));
               }
             }
          // now we are in r3

          for(i=0; i<Ws; i++)
             {
             r=nnm[dirgeo(r,dir,stvol)];
             if(i==0)
               {
               aux*=conj(lattice[r][dir]);
               }
             else
               {
               aux*=conj(multihit(lattice, nnp, nnm, r, dir, epsilon, beta, nhits, stvol));
               }
             }
          } 
        ris+=creal(aux);
        }
     }

  ris/=(double)stvol;
  ris/=(double)(STDIM-1);

  return ris;
  }


// main
int main(int argc, char **argv)
   {
   int i, Nt, Ns, Wt, Ws;
   double beta, rand, W;
   double complex **lattice;
   long int *nnp, *nnm;
   long int iter, sample, r, stvolume, acc, count;

   char datafile[STRING_LENGTH];
   FILE *fp;

   const int overrelax=5;
   const int measevery=50;
   const int unitarizeevery=10;
   const double epsilon=1;
   
   const unsigned long int seed1=(unsigned long int) time(NULL);
   const unsigned long int seed2=seed1+127;

   if(argc != 6)
     {
     fprintf(stdout, "How to use this program:\n");
     fprintf(stdout, "  %s Nt Ns beta sample datafile\n\n", argv[0]);
     fprintf(stdout, "  Nt = temporal size of the lattice\n");
     fprintf(stdout, "  Ns = spatial size of the lattice (space-time dimension defined by macro STDIM)\n");
     fprintf(stdout, "  beta = coupling\n");
     fprintf(stdout, "  sample = number of drawn to be extracted\n");
     fprintf(stdout, "  datafile = name of the file on which to write the data\n\n");
     fprintf(stdout, "Compiled for:\n");
     fprintf(stdout, "  dimensionality = %d\n\n", STDIM);
     fprintf(stdout, "Output:\n");
     fprintf(stdout, "  Wilson loops (Wt, Ws) using the following order\n");
     fprintf(stdout, "  for(Ws=1; Ws<=Ns/4; Ws++){\n");
     fprintf(stdout, "      for(Wt=1; Wt<=MIN(Nt/4,8); Wt++){\n");
     fprintf(stdout, "         Ws, Wt Wilson loop }}\n");

     return EXIT_SUCCESS;
     }
   else
     {  
     // read input values 
     Nt=atoi(argv[1]);
     Ns=atoi(argv[2]);
     beta=atof(argv[3]);
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

   // allocate the lattice
   // and next neighbors: nnp[dirgeo(r, i, volume)]= next neighbor in positive "i" direction of site r 
   lattice=(double complex **)malloc((unsigned long int)(stvolume)*sizeof(double complex *));
   if(lattice == NULL)
     {
     fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
     return EXIT_FAILURE;
     }
   for(r=0; r<stvolume; r++)
      {
      lattice[r]=(double complex *)malloc((unsigned long int)(STDIM)*sizeof(double complex));
      if(lattice[r] == NULL)
        {
        fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
        }
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
      for(i=0; i<STDIM; i++)
         {
         lattice[r][i]=1.0+I*0.0;
         }
      }

   // open data file
   fp=fopen(datafile, "w");
   if(fp==NULL)
     {
     fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
     return EXIT_FAILURE;
     }

   acc=0.0;
   count=0;
   for(iter=0; iter<sample; iter++)
      {
      rand=myrand(); 

      if(rand<1.0/(double)overrelax)
        {
        count++;
        // metropolis
        for(r=0; r<stvolume; r++)
           {
           for(i=0; i<STDIM; i++)
              {
              acc += metropolis(lattice, nnp, nnm, r, i, epsilon, beta, stvolume);
              }
           }
         }
      else
        {
        // overrelaxation
        for(r=0; r<stvolume; r++)
           {
           for(i=0; i<STDIM; i++)
              {
              overrelaxation(lattice, nnp, nnm, r, i, stvolume);
              }
           }
        } 

      // reunitarize
      if(iter%unitarizeevery==0)
        {
        for(r=0; r<stvolume; r++)
           {
           for(i=0; i<STDIM; i++)
              {
              lattice[r][i]/=sqrt(cabs(lattice[r][i]));
              }
            }
        }

      // perform measures
      if(iter%measevery==0)
        {
        // perform measures
        for(Ws=1; Ws<=Ns/4; Ws++)
           {
           for(Wt=1; Wt<=MIN(Nt/4, 8); Wt++)
              {
              W=Wilsonloop(lattice, nnp, nnm, Wt, Ws, epsilon, beta, stvolume);
              fprintf(fp, "%.12f ", W);
              }
           }
        fprintf(fp,"\n");
        }
      }

   printf("Acceptance rate %f\n", (double)acc/(double)count/(double)stvolume/(double)STDIM);

   // close datafile
   fclose(fp);

   // free memory
   for(r=0; r<stvolume; r++)
      {
      free(lattice[r]);
      }
   free(lattice);
   free(nnp);
   free(nnm);

   return EXIT_SUCCESS;
   }


