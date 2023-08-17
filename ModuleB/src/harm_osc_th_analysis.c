#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/read_data.h"

#define STRING_LENGTH 50

// REMEMBER: the data file is a two column file, each raw corresponding to measures taken on the same configuraton
// first column = x (averaged on Nt)
// second column = x^2 (averaged on Nt)
// thir columns = K (averaged on Nt)

// compute the jacknife samples of <x>, <x^2>, <Knaive>, <Hnaive>, <H>
void computejack(double * restrict datajack, 
                 double const * const restrict data, 
                 long int numberofbins, 
                 int binsize,
                 double eta)
  {
  long int i, r;
  int j;
  const long int sampleeff=numberofbins*(long int) binsize;
  double xtot, x2tot, Ktot;
  double x, x2, K;

  xtot=0.0;
  x2tot=0.0;
  Ktot=0.0;

  for(i=0; i<sampleeff; i++)
     {
     xtot +=data[3*i+0];
     x2tot+=data[3*i+1];
     Ktot +=data[3*i+2];
     }

  for(i=0; i<numberofbins; i++)
     {
     x=xtot;
     x2=x2tot;
     K=Ktot;

     for(j=0; j<binsize; j++)
        {
        r=i*binsize+j;

        x -=data[3*r+0];
        x2-=data[3*r+1];
        K -=data[3*r+2];
        }

     x/=(double)((numberofbins-1)*binsize); 
     x2/=(double)((numberofbins-1)*binsize);
     K/=(double)((numberofbins-1)*binsize);
  
     datajack[4*i+0]=x;
     datajack[4*i+1]=x2;
     datajack[4*i+2]=K;
     datajack[4*i+3]=0.5*x2 -K +0.5/eta;
     }
  }


// main
int main(int argc, char **argv)
    {
    int therm, binsize, j;
    long int sample, numberofbins, sampleeff, i, Nt;
    double simbeta, eta;
    double *data, *datajack, ris[4], err[4];   // 4 becuse there are 4 observables
    char datafile[STRING_LENGTH];

    if(argc != 6)
      {
      fprintf(stdout, "How to use this program:\n");
      fprintf(stdout, "  %s therm binsize datafile simbeta Nt\n\n", argv[0]);
      fprintf(stdout, "  therm = number of lines to be discarded as thermalization\n");
      fprintf(stdout, "  binsize = size of the bin to be used in binning/blocking\n");
      fprintf(stdout, "  datafile = name of the data file to be analyzed\n");
      fprintf(stdout, "  simbeta = hbar*omega/(k_B T)\n");
      fprintf(stdout, "  Nt = number of temporal steps\n\n");
      fprintf(stdout, "Output:\n");
      fprintf(stdout, "  <x>, err, <x^2>, err, <Knaive>, err, <H>, err\n");
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
    simbeta=atof(argv[4]);
    Nt=atol(argv[5]);

    if(Nt<=1)
      {
      fprintf(stderr, "'Nt' must be larger than 1\n");
      return EXIT_FAILURE;
      }

    // definition of eta
    eta=simbeta/(double) Nt;

    if(binsize<=0)
      {
      fprintf(stderr, "'binsize' must be positive\n");
      return EXIT_FAILURE;
      }

    // determine the length of the file
    sample=linecounter_mc(datafile, 3);

    // initialize numberofbins and sampleeff
    numberofbins=(sample-therm)/binsize;
    sampleeff=numberofbins*binsize;

    // allocate data arrays
    data=(double *)malloc((unsigned long int)(3*sampleeff)*sizeof(double)); // 3 columns!
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
    readdata_mc(datafile, therm, sampleeff, data, 3);

    // compute jackknife resamplings
    computejack(datajack, data, numberofbins, binsize, eta);

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

