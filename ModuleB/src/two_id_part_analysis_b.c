#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/read_data.h"

#define STRING_LENGTH 50

// REMEMBER: the data file is a 5 columns file, each raw corresponds x1, x1^2, x2, x2^2, twisted(=0,1)

// compute the jacknife samples 
void computejack(double * restrict datajack, 
                 double const * const restrict data, 
                 long int numberofbins, 
                 int binsize)
  {
  long int i, r;
  int j;
  const long int sampleeff=numberofbins*(long int) binsize;
  double x1tot, x2tot, r2tot, signtot;
  double x1, x2, r2, sign;

  x1tot=0.0;
  x2tot=0.0;
  r2tot=0.0;
  signtot=0.0;

  for(i=0; i<sampleeff; i++)
     {
     x1tot+=data[5*i+0];
     x2tot+=data[5*i+2];
     r2tot+=(data[5*i+1]+data[5*i+3]);
     signtot+=(1.0-2.0*data[5*i+4]); 
     }

  for(i=0; i<numberofbins; i++)
     {
     x1=x1tot;
     x2=x2tot;
     r2=r2tot;
     sign=signtot;

     for(j=0; j<binsize; j++)
        {
        r=i*binsize+j;

        x1-=data[5*r+0];
        x2-=data[5*r+2];
        r2-=(data[5*r+1]+data[5*r+3]);
        sign-=(1.0-2.0*data[5*r+4]);
        }

     x1/=(double)((numberofbins-1)*binsize); 
     x2/=(double)((numberofbins-1)*binsize);
     r2/=(double)((numberofbins-1)*binsize);
     sign/=(double)((numberofbins-1)*binsize);
 
     datajack[4*i+0]=x1;
     datajack[4*i+1]=x2;
     datajack[4*i+2]=r2;
     datajack[4*i+3]=sign;
     }
  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize, j;
    long int sample, numberofbins, sampleeff, i;
    double *data, *datajack, ris[4], err[4];   // 4 becuse there are 4 observables
    char datafile[STRING_LENGTH];

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize datafile\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed\n\n");
      fprintf(stdout, "Output [Bosonic case]:\n");
      fprintf(stdout, "  <x1>, err, <x2>, err, <x1^2+x2^2>, err, <(-1)^{perm}>, err\n");
      fprintf(stdout, "  computed using binning and jackknife\n\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input
      therm=atoi(argv[1]);
      binsize=atoi(argv[2]);

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

    if(binsize<=0)
      {
      fprintf(stderr, "'binsize' must be positive\n");
      return EXIT_FAILURE;
      }

    // determine the length of the file
    sample=linecounter_mc(datafile, 5);

    // initialize numberofbins and sampleeff
    numberofbins=(sample-therm)/binsize;
    sampleeff=numberofbins*binsize;

    // allocate data arrays
    data=(double *)malloc((unsigned long int)(5*sampleeff)*sizeof(double));
    if(data==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // allocate jackknife samples
    datajack=(double *)malloc((unsigned long int)(4*numberofbins)*sizeof(double)); // 4 because there are 4 observables
    if(datajack==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize data
    readdata_mc(datafile, therm, sampleeff, data, 5);

    // compute jackknife resamplings
    computejack(datajack, data, numberofbins, binsize);

    // compute average
    for(j=0; j<4; j++)
       {
       ris[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          ris[j]+=datajack[4*i+j];
          }
       ris[j]/=(double)numberofbins;
       }

    // compute error
    for(j=0; j<4; j++)
       {
       err[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          err[j]+=pow(ris[j]-datajack[4*i+j], 2.0);
          }
       // this corrects for a factor that is irrelevant but we leave it just for clarity
       err[j]*=(double)(numberofbins-1);
       err[j]/=(double)numberofbins;
       err[j]=sqrt(err[j]);
       }

    // free data arrays
    free(data);
    free(datajack);

    for(j=0; j<4; j++)
       {
       printf("%.12f %.12f ", ris[j], err[j]);
       }
    printf("\n");

    return EXIT_SUCCESS;
    }

