#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/read_data.h"

#define STRING_LENGTH 50

// REMEMBER: the data file is a two column file, each raw corresponding to measures taken on the same configuraton
// first column = energy per site
// second column = magnetization per site

// compute the jacknife samples of <E>, <E^2>-<E>^2, <M>, <|M|>, <M^2>-<|M|^2>, <M2>, <M^4>/<M^2>^2
void computejack(double * restrict datajack, 
                 double const * const restrict data, 
                 long int numberofbins, 
                 int binsize)
  {
  long int i, r;
  int j;
  const long int sampleeff=numberofbins*(long int) binsize;
  double Etot, E2tot, Mtot, Mabstot, M2tot, M4tot;
  double E, E2, M, Mabs, M2, M4;

  Etot=0.0;
  E2tot=0.0;
  Mtot=0.0;
  Mabstot=0.0;
  M2tot=0.0;
  M4tot=0.0;

  for(i=0; i<sampleeff; i++)
     {
     Etot+=data[2*i+0];
     E2tot+=pow(data[2*i+0], 2.0);
     Mtot+=data[2*i+1];
     Mabstot+=fabs(data[2*i+1]);
     M2tot+=pow(data[2*i+1],2.0);
     M4tot+=pow(data[2*i+1],4.0);
     }

  for(i=0; i<numberofbins; i++)
     {
     E=Etot;
     E2=E2tot;
     M=Mtot;
     Mabs=Mabstot;
     M2=M2tot;
     M4=M4tot;

     for(j=0; j<binsize; j++)
        {
        r=i*binsize+j;

        E-=data[2*r+0];
        E2-=pow(data[2*r+0], 2.0);
        M-=data[2*r+1];
        Mabs-=fabs(data[2*r+1]);
        M2-=pow(data[2*r+1],2.0);
        M4-=pow(data[2*r+1],4.0);
        }

     E/=(double)((numberofbins-1)*binsize); 
     E2/=(double)((numberofbins-1)*binsize);
     M/=(double)((numberofbins-1)*binsize);
     Mabs/=(double)((numberofbins-1)*binsize);
     M2/=(double)((numberofbins-1)*binsize);
     M4/=(double)((numberofbins-1)*binsize);
  
     datajack[7*i+0]=E;
     datajack[7*i+1]=E2-E*E;
     datajack[7*i+2]=M;
     datajack[7*i+3]=Mabs;
     datajack[7*i+4]=M2-Mabs*Mabs;
     datajack[7*i+5]=M2;
     datajack[7*i+6]=M4/(M2*M2);
     }
  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize, j;
    long int sample, numberofbins, sampleeff, i;
    double *data, *datajack, ris[7], err[7];   // 7 becuse there are 7 observables
    char datafile[STRING_LENGTH];

    if(argc != 4)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize datafile\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed (2 columns: energy and magnetization per site)\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  <E>, err, <E^2>-<E>^2, err, <M>, err, <|M|>, err,");
      fprintf(stdout, " <M^2>-<|M|>^2, err, <M^2>, err, <M^4>/<M^2>^2, err\n");
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
    sample=linecounter_mc(datafile, 2);

    // initialize numberofbins and sampleeff
    numberofbins=(sample-therm)/binsize;
    sampleeff=numberofbins*binsize;

    // allocate data arrays
    data=(double *)malloc((unsigned long int)(2*sampleeff)*sizeof(double));
    if(data==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // allocate jackknife samples
    datajack=(double *)malloc((unsigned long int)(7*numberofbins)*sizeof(double)); // 7 because there are 7 observables
    if(datajack==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize data
    readdata_mc(datafile, therm, sampleeff, data, 2);

    // compute jackknife resamplings
    computejack(datajack, data, numberofbins, binsize);

    // compute average
    for(j=0; j<7; j++)
       {
       ris[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          ris[j]+=datajack[7*i+j];
          }
       ris[j]/=(double)numberofbins;
       }

    // compute error
    for(j=0; j<7; j++)
       {
       err[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          err[j]+=pow(ris[j]-datajack[7*i+j], 2.0);
          }
       // this corrects for a factor that is irrelevant but we leave it just for clarity
       err[j]*=(double)(numberofbins-1);
       err[j]/=(double)numberofbins;
       err[j]=sqrt(err[j]);
       }

    // free data arrays
    free(data);
    free(datajack);

    for(j=0; j<7; j++)
       {
       printf("%f %f ", ris[j], err[j]);
       }
    printf("\n");

    return EXIT_SUCCESS;
    }

