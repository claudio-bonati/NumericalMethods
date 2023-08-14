#ifndef READ_DATA_H 
#define READ_DATA_H 

// determine the length of the file (single column!)
long int linecounter_sc(char const * const filename);

// initialize data (single column!)
void readdata_sc(char const * const filename, int therm, long int sampleeff, double * restrict data);

// determine the length of the file with 'col' columns
long int linecounter_mc(char const * const filename, int col);

// initialize data from file with 'col' columns
// *data is structured as data[i*col+j]= j-th column of the i-th raw 
void readdata_mc(char const * const filename, int therm, long int sampleeff, double * restrict data, int col);

#endif 
