#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/read_data.h"

#define STRING_LENGTH 50

// REMEMBER: the data file is a single column file, each raw corresponds to the measure of Q taken on one configuraton

// compute the jacknife samples of <Q>, <Q^2>, -(<Q^4>-3<Q^2>^2)/(12<Q^2>)
void computejack(double * restrict datajack, 
                 double const * const restrict data, 
                 long int numberofbins, 
                 int binsize)
  {
  long int i, r;
  int j;
  const long int sampleeff=numberofbins*(long int) binsize;
  double Qtot, Q2tot, Q4tot;
  double Q, Q2, Q4;

  Qtot=0.0;
  Q2tot=0.0;
  Q4tot=0.0;

  for(i=0; i<sampleeff; i++)
     {
     Qtot+=data[i];
     Q2tot+=pow(data[i], 2.0);
     Q4tot+=pow(data[i], 4.0);
     }

  for(i=0; i<numberofbins; i++)
     {
     Q=Qtot;
     Q2=Q2tot;
     Q4=Q4tot;

     for(j=0; j<binsize; j++)
        {
        r=i*binsize+j;

        Q-=data[r];
        Q2-=pow(data[r], 2.0);
        Q4-=pow(data[r], 4.0);
        }

     Q/=(double)((numberofbins-1)*binsize); 
     Q2/=(double)((numberofbins-1)*binsize);
     Q4/=(double)((numberofbins-1)*binsize);
  
     datajack[3*i+0]=Q;
     datajack[3*i+1]=Q2;
     datajack[3*i+2]=-(Q4-3.0*Q2*Q2)/(12.0*Q2);
     }
  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize, j;
    long int sample, numberofbins, sampleeff, i;
    double *data, *datajack, ris[3], err[3];   // 3 becuse there are 3 observables
    char datafile[STRING_LENGTH];

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize datafile\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed (1 column, Q)\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  <Q>, err, <Q^2>, err, -(<Q^4>-3<Q^2>^2)/(12<Q^2>), err\n");
      fprintf(stdout, "  computed using binning and jackknife\n");

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
    sample=linecounter_sc(datafile);

    // initialize numberofbins and sampleeff
    numberofbins=(sample-therm)/binsize;
    sampleeff=numberofbins*binsize;

    // allocate data arrays
    data=(double *)malloc((unsigned long int)(sampleeff)*sizeof(double));
    if(data==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // allocate jackknife samples
    datajack=(double *)malloc((unsigned long int)(3*numberofbins)*sizeof(double)); // 3 because there are 3 observables
    if(datajack==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize data
    readdata_sc(datafile, therm, sampleeff, data);

    // compute jackknife resamplings
    computejack(datajack, data, numberofbins, binsize);

    // compute average
    for(j=0; j<3; j++)
       {
       ris[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          ris[j]+=datajack[3*i+j];
          }
       ris[j]/=(double)numberofbins;
       }

    // compute error
    for(j=0; j<3; j++)
       {
       err[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          err[j]+=pow(ris[j]-datajack[3*i+j], 2.0);
          }
       // this corrects for a factor that is irrelevant but we leave it just for clarity
       err[j]*=(double)(numberofbins-1);
       err[j]/=(double)numberofbins;
       err[j]=sqrt(err[j]);
       }

    // free data arrays
    free(data);
    free(datajack);

    for(j=0; j<3; j++)
       {
       printf("%f %f ", ris[j], err[j]);
       }
    printf("\n");

    return EXIT_SUCCESS;
    }

