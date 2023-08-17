CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99


binaries=harm_osc_th harm_osc_th_analysis harm_osc_spectrum harm_osc_spectrum_analysis part_on_a_circ part_on_a_circ_analysis part_on_a_circ_par_temp two_id_part two_id_part_analysis_b two_id_part_analysis_f

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

harm_osc_spectrum_analysis: ../src/harm_osc_spectrum_analysis.c read_data.o
	$(CC) $(CFLAGS) -c ../src/harm_osc_spectrum_analysis.c 
	$(CC) $(CFLAGS) harm_osc_spectrum_analysis.o read_data.o -o $@ -lm

part_on_a_circ: ../src/part_on_a_circ.c pcg32min.o random.o
	$(CC) $(CFLAGS) -c ../src/part_on_a_circ.c 
	$(CC) $(CFLAGS) part_on_a_circ.o pcg32min.o random.o -o $@ -lm

part_on_a_circ_analysis: ../src/part_on_a_circ_analysis.c read_data.o
	$(CC) $(CFLAGS) -c ../src/part_on_a_circ_analysis.c 
	$(CC) $(CFLAGS) part_on_a_circ_analysis.o read_data.o -o $@ -lm

part_on_a_circ_par_temp: ../src/part_on_a_circ_par_temp.c pcg32min.o random.o
	$(CC) $(CFLAGS) -c ../src/part_on_a_circ_par_temp.c 
	$(CC) $(CFLAGS) part_on_a_circ_par_temp.o pcg32min.o random.o -o $@ -lm

two_id_part: ../src/two_id_part.c boxmuller.o pcg32min.o random.o
	$(CC) $(CFLAGS) -c ../src/two_id_part.c 
	$(CC) $(CFLAGS) two_id_part.o boxmuller.o pcg32min.o random.o -o $@ -lm

two_id_part_analysis_b: ../src/two_id_part_analysis_b.c read_data.o
	$(CC) $(CFLAGS) -c ../src/two_id_part_analysis_b.c 
	$(CC) $(CFLAGS) two_id_part_analysis_b.o read_data.o -o $@ -lm

two_id_part_analysis_f: ../src/two_id_part_analysis_f.c read_data.o
	$(CC) $(CFLAGS) -c ../src/two_id_part_analysis_f.c 
	$(CC) $(CFLAGS) two_id_part_analysis_f.o read_data.o -o $@ -lm

.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o
