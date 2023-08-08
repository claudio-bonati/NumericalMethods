CC=gcc
CFLAGS=-O2 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99

binaries=harm_osc_th harm_osc_th_analysis harm_osc_spectrum

all: $(binaries)

boxmuller.o: ../lib/boxmuller.c ../include/boxmuller.h
	$(CC) $(CFLAGS) -c ../lib/boxmuller.c

random.o: ../lib/random.c ../include/pcg32min.h ../include/random.h
	$(CC) $(CFLAGS) -c ../lib/random.c

pcg32min.o: ../lib/pcg32min.c ../include/pcg32min.h 
	$(CC) $(CFLAGS) -c ../lib/pcg32min.c

read_data.o: ../lib/read_data.c ../include/read_data.h
	$(CC) $(CFLAGS) -c ../lib/read_data.c

harm_osc_th: ../src/harm_osc_th.c boxmuller.o pcg32min.o random.o
	$(CC) $(CFLAGS) -c ../src/harm_osc_th.c 
	$(CC) $(CFLAGS) harm_osc_th.o boxmuller.o pcg32min.o random.o -o $@ -lm

harm_osc_th_analysis: ../src/harm_osc_th_analysis.c read_data.o
	$(CC) $(CFLAGS) -c ../src/harm_osc_th_analysis.c 
	$(CC) $(CFLAGS) harm_osc_th_analysis.o read_data.o -o $@ -lm

harm_osc_spectrum: ../src/harm_osc_spectrum.c boxmuller.o pcg32min.o random.o
	$(CC) $(CFLAGS) -c ../src/harm_osc_spectrum.c 
	$(CC) $(CFLAGS) harm_osc_spectrum.o boxmuller.o pcg32min.o random.o -o $@ -lm


.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o
