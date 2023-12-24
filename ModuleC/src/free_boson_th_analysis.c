#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/read_data.h"

#define STRING_LENGTH 50

// REMEMBER: the data file is a 3 column file, each raw corresponding to measures taken on the same configuraton
// first column = O1
// second column = O2
// third column = O3

// compute the jacknife samples of <O1+O2-O3>/(2 Nt^{stdim})
void computejack(double * restrict datajack, 
                 double const * const restrict data, 
                 long int numberofbins, 
                 int binsize) 
  {
  long int i, r;
  int j;
  const long int sampleeff=numberofbins*(long int) binsize;
  double ris1tot, ris2tot;
  double ris1, ris2;

  ris1tot=0.0;
  ris2tot=0.0;

  for(i=0; i<sampleeff; i++)
     {
     ris1tot+=data[3*i+0]+data[3*i+1]-data[3*i+2]; // O1+O2-O3
     ris2tot+=data[3*i+0]; // O1
     }

  for(i=0; i<numberofbins; i++)
     {
     ris1=ris1tot;
     ris2=ris2tot;

     for(j=0; j<binsize; j++)
        {
        r=i*binsize+j;

        ris1-=(data[3*r+0]+data[3*r+1]-data[3*r+2]); // O1+O2-O3
        ris2-=data[3*r+0]; // O1
        }

     ris1/=(double)((numberofbins-1)*binsize); 
     ris2/=(double)((numberofbins-1)*binsize); 

     ris1/=2.0;

     datajack[2*i+0]=ris1;
     datajack[2*i+1]=ris2;
     }
  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize, j, Nt;
    long int sample, numberofbins, sampleeff, i;
    double *data, *datajack, ris[2], err[2], hatm;  // two observables
    char datafile[STRING_LENGTH];

    if(argc != 6)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize Nt hatm datafile\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  Nt = temporal extent of the lattice\n");
      fprintf(stdout, "  hatm = a*m\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  Nt hatm <O1+O2-O3>/2 err <O1> err\n");
      fprintf(stdout, "  computed using binning and jackknife, where the two observables\n");
      fprintf(stdout, "  once renormalized and multiplied by Nt^{stdim} are related to\n");
      fprintf(stdout, "  the energy density and to the trace of the energy momentum tensor.\n");

      return EXIT_SUCCESS;
      }
    else
      {
      // read input
      therm=atoi(argv[1]);
      binsize=atoi(argv[2]);
      Nt=atoi(argv[3]);
      hatm=atof(argv[4]);

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

    if(binsize<=0)
      {
      fprintf(stderr, "'binsize' must be positive\n");
      return EXIT_FAILURE;
      }

    // determine the length of the file
    sample=linecounter_mc(datafile, 3); // 3 columns

    // initialize numberofbins and sampleeff
    numberofbins=(sample-therm)/binsize;
    sampleeff=numberofbins*binsize;

    // allocate data arrays
    data=(double *)malloc((unsigned long int)(3*sampleeff)*sizeof(double));  // 3 columns
    if(data==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // allocate jackknife samples
    datajack=(double *)malloc((unsigned long int)(2*numberofbins)*sizeof(double)); // 2 because there are 2 observables to be computed
    if(datajack==NULL)
      {
      fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
      }

    // initialize data
    readdata_mc(datafile, therm, sampleeff, data, 3);

    // compute jackknife resamplings
    computejack(datajack, data, numberofbins, binsize);

    // compute average
    for(j=0; j<2; j++)
       {
       ris[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          ris[j]+=datajack[2*i+j];
          }
       ris[j]/=(double)numberofbins;
       }

    // compute error
    for(j=0; j<2; j++)
       {
       err[j]=0.0;
       for(i=0; i<numberofbins; i++)
          {
          err[j]+=pow(ris[j]-datajack[2*i+j], 2.0);
          }
       // this corrects for a factor that is irrelevant but we leave it just for clarity
       err[j]*=(double)(numberofbins-1);
       err[j]/=(double)numberofbins;
       err[j]=sqrt(err[j]);
       }

    // free data arrays
    free(data);
    free(datajack);

    printf("%d %.12f ", Nt, hatm);
    for(j=0; j<2; j++)
       {
       printf("%.12f %.12f ", ris[j], err[j]);
       }
    printf("\n");

    return EXIT_SUCCESS;
    }

