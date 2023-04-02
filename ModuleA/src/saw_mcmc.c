#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/random.h"

#define STRING_LENGTH 50

// a trial is generate by selecting a pivot point
// and rotating/reflecting the path starting from the pivot point
void generate_trial(int **path_trial, int * const * const path, int L)
  {
  int i, tmp;
  const int pivot=(int)(L*myrand());

  for(i=0; i<=pivot; i++)
     {
     path_trial[i][0]=path[i][0];
     path_trial[i][1]=path[i][1];
     }

  tmp=(int)(7*myrand());
  
  switch(tmp)
    {
    case 0: // pi/2 rotation
      for(i=pivot+1; i<L; i++)
         {
         path_trial[i][0]=path[pivot][0] -(path[i][1]-path[pivot][1]);
         path_trial[i][1]=path[pivot][1] +(path[i][0]-path[pivot][0]);
         }
      break;
    
    case 1: // -pi/2 rotation
      for(i=pivot+1; i<L; i++)
         {
         path_trial[i][0]=path[pivot][0] +(path[i][1]-path[pivot][1]);
         path_trial[i][1]=path[pivot][1] -(path[i][0]-path[pivot][0]);
         }
      break;

    case 2: // pi rotation
      for(i=pivot+1; i<L; i++)
         {
         path_trial[i][0]=path[pivot][0] -(path[i][0]-path[pivot][0]);
         path_trial[i][1]=path[pivot][1] -(path[i][1]-path[pivot][1]);
         }
      break;

    case 3: // reflection along x=const 
      for(i=pivot+1; i<L; i++)
         {
         path_trial[i][0]=path[i][0];
         path_trial[i][1]=path[pivot][1] -(path[i][1]-path[pivot][1]);
         }
      break;

    case 4: // reflection along y=const
      for(i=pivot+1; i<L; i++)
         {
         path_trial[i][0]=path[pivot][0] -(path[i][0]-path[pivot][0]);
         path_trial[i][1]=path[i][1];
         }
      break;

    case 5: // reflection along y=x
      for(i=pivot+1; i<L; i++)
         {
         path_trial[i][0]=path[pivot][0] +(path[i][1]-path[pivot][1]);
         path_trial[i][1]=path[pivot][1] +(path[i][0]-path[pivot][0]);
         }
      break;

    case 6: // reflection along y=-x
      for(i=pivot+1; i<L; i++)
         {
         path_trial[i][0]=path[pivot][0] -(path[i][1]-path[pivot][1]);
         path_trial[i][1]=path[pivot][1] -(path[i][0]-path[pivot][0]);
         }
      break;
    }
  }


// if no intersection is found 
// path <- path_trial and return 1
// else return 0
int check_and_update(int **path, int * const * const path_trial, int L)
  {
  int i, j, k, success;

  success=1;

  // check for intersections
  for(i=0; i<L && success==1; i++)
     {
     for(j=i+1; j<i+L; j++)
        {
        k = j % L; // in this way runs on all values between 0<= <L different from i

        if(path_trial[i][0]==path_trial[k][0] && path_trial[i][1]==path_trial[k][1])
          {
          success=0;
          }
        }
     }

  if(success==1)
    {
    for(i=0; i<L; i++)
       {
       path[i][0]=path_trial[i][0];
       path[i][1]=path_trial[i][1];
       }
    }

  return success;
  }


// main
int main(int argc, char **argv)
    {
    int **path, **path_trial;
    long int sample, iteration, accepted;
    int i, d2, length, L;
    char datafile[STRING_LENGTH];
    FILE *fp;
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s length sample datafile\n\n", argv[0]);
      fprintf(stdout, "  length = length of the 2d self-avoiding random walk\n");
      fprintf(stdout, "  sample = number of walk to be generated\n");
      fprintf(stdout, "  datafile = file where to write data\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  r^2 for every sample (single column)\n");

      /* 
      * The algorithm used is the pivot algorithm, see 
      *
      * Neal Madras, Alan D. Sokal 
      * The pivot algorithm: A highly efficient Monte Carlo method for the self-avoiding walk
      * Journal of Statistical Physics vol. 50, pages 109-186 (1988)
      */

      return EXIT_SUCCESS;
      }
    else
      {  
      // read input values 

      length=atoi(argv[1]);
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

    // initialize the random number generator
    myrand_init(seed1, seed2);

    L=length+1; // L is the number of points of the path

    // allocate the paths 
    // path[i][j]= j-th component of the i-th point
    // j=0,1 corresponds to x,y
    // all paths start from (0,0)
    path=(int **)malloc((unsigned long int)L*sizeof(int*));
    if(path == NULL)
      {
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    for(i=0; i<L; i++)
       {
       path[i]=(int *)malloc(2*sizeof(int));
       if(path[i] == NULL)
         {
         fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
         return EXIT_FAILURE;
         }
       }

    path_trial=(int **)malloc((unsigned long int)L*sizeof(int*));
    if(path_trial == NULL){
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    for(i=0; i<L; i++)
       {
       path_trial[i]=(int *)malloc(2*sizeof(int));
       if(path_trial[i] == NULL)
         {
         fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
         return EXIT_FAILURE;
         }
       }

    // initialize the path with a straight line
    for(i=0; i<L; i++)
       {
       path[i][0]=i;
       path[i][1]=0;
       }

    // open datafile for writing
    fp=fopen(datafile, "w");

    accepted=0;
    // loop on iterations
    for(iteration=0; iteration<sample; iteration++) 
       {
       generate_trial(path_trial, path, L);
       accepted+=check_and_update(path, path_trial, L);

       d2=path[L-1][0]*path[L-1][0] + path[L-1][1]*path[L-1][1];
       fprintf(fp, "%d\n", d2);   
       }

    // close datafile
    fclose(fp);

    // free the path
    for(i=0; i<L; i++)
       {
       free(path[i]);
       free(path_trial[i]);
       }
    free(path);
    free(path_trial);


    // print simulation details
    printf("Acceptance rate: %f\n", (double)accepted / (double)sample);


    return EXIT_SUCCESS;
    }
