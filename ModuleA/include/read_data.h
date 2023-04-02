#ifndef READ_DATA_H 

// determine the length of the file (single column!)
long int linecounter_sc(char const * const filename);

// initialize data (single column!)
void readdata_sc(char const * const filename, int therm, long int sampleeff, double *data);

#endif 
