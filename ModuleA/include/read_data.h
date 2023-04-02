#ifndef READ_DATA_H 

// determine the length of the file (single column!)
long int linecounter_sc(char const * const filename);

// initialize data (single column!)
void readdata_sc(char const * const filename, int therm, long int sampleeff, double *data);

// determine the length of the file with 'col' columns
long int linecounter_mc(char const * const filename, int col);

// initialize data from file with 'col' columns
// **data has indices ordered as data[raw][col]
void readdata_mc(char const * const filename, int therm, long int sampleeff, double **data, int col);

#endif 
