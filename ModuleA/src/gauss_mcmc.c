#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/random.h"

#define STRING_LENGTH 50

// main
int main(int argc, char **argv)
    {
    long int i, sample, acc;
    double step, state, trial, deltaenergy;
    double x, x2, x4;
    char datafile[STRING_LENGTH];
    FILE *fp;
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s step sample datafile\n\n", argv[0]);
      fprintf(stdout, "  step = size of the random move\n");
      fprintf(stdout, "  sample = number of draws\n");
      fprintf(stdout, "  datafile = name of the file in which to print the data\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  a single column of 'sample' number with normal gaussian distribution in 'datafile'\n");
      fprintf(stdout, "  the acceptance rate\n");
      fprintf(stdout, "  the [naive!] values of <x> and <x^2> with their errors (sample/10 thermalization steps)\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input
      step=atof(argv[1]);
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

    if(step<=0)
      {
      fprintf(stderr, "'step' must be positive\n");
      return EXIT_FAILURE;
      }

    if(sample<=0)
      {
      fprintf(stderr, "'sample' must be positive\n");
      return EXIT_FAILURE;
      }

    // initialize random number generator
    myrand_init(seed1, seed2);

    // open data file
    fp=fopen(datafile, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
      return EXIT_FAILURE;
      }
    
    // intial value = 1
    state=1;

    // initialize the acceptance and the averages
    acc=0;
    x=0.0;
    x2=0.0;
    x4=0.0;

    // loop on iterations
    for(i=0; i<sample; i++)
       {
       if(i>sample/10)
	 {
         x+=state;
         x2+=pow(state,2.0);
         x4+=pow(state,4.0);
	 }

       fprintf(fp, "%f\n", state);
       trial=state+step*(1.0-2.0*myrand());

       // the target distribution is exp[-x^2/2]/Z
       deltaenergy=pow(trial,2.0)/2.0 - pow(state,2.0)/2.0;

       if(deltaenergy<0) // the trial is accepted
         {
         state=trial;
         acc++;
         }
       else if(myrand()<exp(-deltaenergy)) // the trial is accepted
              {
              state=trial;
              acc++;
              } 
       }

    // close datadile
    fclose(fp);

    // normalize averages
    x/=((double)sample * 9.0/10.0);
    x2/=((double)sample * 9.0/10.0);
    x4/=((double)sample * 9.0/10.0);

    printf("Acceptance rate=%f\n", (double) acc/(double) sample);
    printf("<x>[naive!]=%f %f\n", x, sqrt((x2-x*x)/((double) sample) * 9.0/10.0) );
    printf("<x^2>[naive!]=%f %f\n", x2, sqrt((x4-x2*x2)/((double) sample) * 9.0/10.0) );

    return EXIT_SUCCESS;
    }

