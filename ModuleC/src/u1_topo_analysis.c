#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/read_data.h"

#define STRING_LENGTH 50

// REMEMBER: the data file is a 2 column file, each raw corresponding to measures taken on the same configuraton
// first column = plaq
// second column = top. charge
// compute the jacknife samples of <plaq>  <top. charge>  top.susc.
void computejack(double * restrict datajack, 
                 double const * const restrict data, 
                 long int numberofbins, 
                 int binsize,
                 long int stvolume) 
  {
  long int i, r;
  int j;
  const long int sampleeff=numberofbins*(long int) binsize;
  double ris1tot, ris2tot, ris3tot;
  double ris1, ris2, ris3;

  ris1tot=0.0;
  ris2tot=0.0;
  ris3tot=0.0;


  for(i=0; i<sampleeff; i++)
     {
     ris1tot+=data[2*i+0]; 
     ris2tot+=data[2*i+1]; 
     ris3tot+=data[2*i+1]*data[2*i+1]; 
     }

  for(i=0; i<numberofbins; i++)
     {
     ris1=ris1tot;
     ris2=ris2tot;
     ris3=ris3tot;


     for(j=0; j<binsize; j++)
        {
        r=i*binsize+j;

        ris1-=data[2*r+0];
        ris2-=data[2*r+1];
        ris3-=data[2*r+1]*data[2*r+1];
        }

     ris1/=(double)((numberofbins-1)*binsize); 
     ris2/=(double)((numberofbins-1)*binsize); 
     ris3/=(double)((numberofbins-1)*binsize); 

     datajack[3*i+0]=ris1;
     datajack[3*i+1]=ris2;
     datajack[3*i+2]=ris3/(double) stvolume;
     }
  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize, j;
    long int sample, numberofbins, sampleeff, i, stvolume;
    double *data, *datajack, ris[3], err[3];  // 3 observables
    char datafile[STRING_LENGTH];

    if(argc != 5)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize stvolume datafile\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  stvolume = space-time volume\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  <plaq>  <(top. charge)>  susc.top.\n");
      fprintf(stdout, "  computed using binning and jackknife\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input
      therm=atoi(argv[1]);
      binsize=atoi(argv[2]);
      stvolume=atol(argv[3]);

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

    if(binsize<=0)
      {
      fprintf(stderr, "'binsize' must be positive\n");
      return EXIT_FAILURE;
      }

    // determine the length of the file
    sample=linecounter_mc(datafile, 2); // 2 columns

    // initialize numberofbins and sampleeff
    numberofbins=(sample-therm)/binsize;
    sampleeff=numberofbins*binsize;

    // allocate data arrays
    data=(double *)malloc((unsigned long int)(2*sampleeff)*sizeof(double));  // 2 columns
    if(data==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // allocate jackknife samples
    datajack=(double *)malloc((unsigned long int)(3*numberofbins)*sizeof(double)); // 3 because there are 3 observables to be computed
    if(datajack==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize data
    readdata_mc(datafile, therm, sampleeff, data, 2);

    // compute jackknife resamplings
    computejack(datajack, data, numberofbins, binsize, stvolume);

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
       printf("%.12f %.12f ", ris[j], err[j]);
       }
    printf("\n");

    return EXIT_SUCCESS;
    }

