#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/random.h"

// try to generate a sald-avoiding walk of given length
// acc=1 if success, else =0
// d2=|start-end|^2 
void saw_generate(int **lattice, 
                 int length,
                 int *d2, 
                 int *acc)
  {
  int i, j, tmp, fail;
  int x0, y0, x, y;
  int step;

  // initialize the lattice with zero (unoccupied points)
  for(i=0; i<2*length+2; i++)
     {
     for(j=0; j<2*length+2; j++)
        {
        lattice[i][j]=0;
        }
     }

   // starting position
   x0=length+1;
   y0=length+1;

   // set the starting position as occupied site
   x=x0;
   y=y0;
   step=0;
   lattice[x][y]=1;

   // first step
   tmp=(int) (4*myrand());
   switch(tmp)
         {
         case 0:
           x++;
           break;
         case 1:
           x--;
           break;
         case 2:
           y++;
           break;
         default:
           y--;
         }
   
   // set as occupied
   lattice[x][y]=1;
   step++;

   fail=0;
   // build the walk
   while(step<length && fail==0)
        {
        // select the next move (avoiding the previous point)
        switch(tmp)
              {
              case 0:
                tmp=2+(int)(3*myrand());
                break;
              case 1:
                tmp=1+(int)(3*myrand());
                break;
              case 2:
                tmp=0+(int)(3*myrand());
                break;
              default:
                tmp=3+(int)(3*myrand());
              }
        tmp=tmp%4;

        // check if the move is possible
        switch(tmp)
              {
              case 0:
                 if(lattice[x+1][y]==0){x++;}
                 else {fail=1;}
                 break;
               case 1:
                 if(lattice[x-1][y]==0){x--;}
                 else {fail=1;}
                 break;
               case 2:
                 if(lattice[x][y+1]==0){y++;}
                 else {fail=1;}
                 break;
               default:
                 if(lattice[x][y-1]==0){y--;}
                 else {fail=1;}
               }

        // set as occupied
        step++;
        lattice[x][y]=1;
        }

   if(fail==0) // i.e. not aborted
     {
     *d2=(x-x0)*(x-x0) + (y-y0)*(y-y0);  // square distance from the start
     *acc=1;
     }
   else
     {
     *acc=0;
     *d2=0;
     }
   }


// main
int main(int argc, char **argv)
    {
    int **lattice;
    int length, d2_loc, acc, i;
    double d2f, d4f;
    long int sample, counter, accepted;
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    if(argc != 3)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s length sample\n\n", argv[0]);
      fprintf(stdout, "  length = length of the self-avoiding random walk\n");
      fprintf(stdout, "  sample = number of walk to be generate\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  length, <r^2>, its error, and the success rate in the generation\n");

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 

      length=atoi(argv[1]);
      sample=atol(argv[2]);
      }

    if(length<=0)
      {
      fprintf(stderr, "'length' must be positive\n");
      return EXIT_FAILURE;
      }

    if(sample<=0)
      {
      fprintf(stderr, "'sample' must be positive\n");
      return EXIT_FAILURE;
      }

    // initialize the random number generator
    myrand_init(seed1, seed2);

    // allocate the lattice 
    // (this is not strictly needed, but makes the code faster
    // at the expence of a larger used memory)
    //
    // lattice=0 free site, lattice=1 occupied site
    lattice=(int **)malloc((long unsigned int)(2*length+2)*sizeof(int*));
    if(lattice == NULL)
      {
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    for(i=0; i<2*length+2; i++)
       {
       lattice[i]=(int *)malloc((long unsigned int)(2*length+2)*sizeof(int));
       if(lattice[i] == NULL)
         {
         fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
         return EXIT_FAILURE;
         }
       }

    // initialize averages and counters
    d2f=0.0;
    d4f=0.0;
    counter=0;
    accepted=0;

    // loop over the saw sample
    while(accepted<sample)
         {
         saw_generate(lattice, length, &d2_loc, &acc);
   
         accepted+=acc;
         counter++;

         if(acc==1)
           {
           d2f+=(double)d2_loc;
           d4f+=(double)d2_loc*(double)d2_loc;
           }
         }

    d2f/=(double) sample;
    d4f/=(double) sample;

    // free the lattice
    for(i=0; i<2*length+2; i++)
       {
       free(lattice[i]);
       }
    free(lattice);

    //print results
    printf("%d ", length);
    printf("%lf %lf ", d2f, sqrt(d4f-d2f*d2f)/sqrt((double) sample));
    printf("%lf\n", (double) accepted/(double) counter);

    return EXIT_SUCCESS;
    }
