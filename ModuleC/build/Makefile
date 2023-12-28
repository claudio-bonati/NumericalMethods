CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99 

binaries = free_boson_th free_boson_th_analysis free_boson_th_doublers free_boson_spectrum free_boson_spectrum_analysis u1_topo u1_topo_analysis u1_conf u1_conf_analysis

all: $(binaries)

boxmuller.o: ../lib/boxmuller.c ../include/boxmuller.h
	$(CC) $(CFLAGS) -c ../lib/boxmuller.c

geometry_st.o: ../lib/geometry_st.c ../include/geometry_st.h 
	$(CC) $(CFLAGS) -c ../lib/geometry_st.c  

pcg32min.o: ../lib/pcg32min.c ../include/pcg32min.h 
	$(CC) $(CFLAGS) -c ../lib/pcg32min.c

random.o: ../lib/random.c ../include/pcg32min.h ../include/random.h
	$(CC) $(CFLAGS) -c ../lib/random.c

read_data.o: ../lib/read_data.c ../include/read_data.h
	$(CC) $(CFLAGS) -c ../lib/read_data.c

free_boson_th: ../src/free_boson_th.c boxmuller.o geometry_st.o random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/free_boson_th.c  
	$(CC) $(CFLAGS) free_boson_th.o boxmuller.o geometry_st.o random.o pcg32min.o -o $@ -lm

free_boson_th_analysis: ../src/free_boson_th_analysis.c read_data.o
	$(CC) $(CFLAGS) -c ../src/free_boson_th_analysis.c 
	$(CC) $(CFLAGS) free_boson_th_analysis.o read_data.o -o $@ -lm

free_boson_th_doublers: ../src/free_boson_th_doublers.c boxmuller.o geometry_st.o random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/free_boson_th_doublers.c  
	$(CC) $(CFLAGS) free_boson_th_doublers.o boxmuller.o geometry_st.o random.o pcg32min.o -o $@ -lm

free_boson_spectrum: ../src/free_boson_spectrum.c boxmuller.o geometry_st.o random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/free_boson_spectrum.c  
	$(CC) $(CFLAGS) free_boson_spectrum.o boxmuller.o geometry_st.o random.o pcg32min.o -o $@ -lm

free_boson_spectrum_analysis: ../src/free_boson_spectrum_analysis.c read_data.o
	$(CC) $(CFLAGS) -c ../src/free_boson_spectrum_analysis.c 
	$(CC) $(CFLAGS) free_boson_spectrum_analysis.o read_data.o -o $@ -lm

u1_topo: ../src/u1_topo.c geometry_st.o random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/u1_topo.c  
	$(CC) $(CFLAGS) u1_topo.o geometry_st.o random.o pcg32min.o -o $@ -lm

u1_topo_analysis: ../src/u1_topo_analysis.c read_data.o
	$(CC) $(CFLAGS) -c ../src/u1_topo_analysis.c 
	$(CC) $(CFLAGS) u1_topo_analysis.o read_data.o -o $@ -lm

u1_conf: ../src/u1_conf.c geometry_st.o random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/u1_conf.c  
	$(CC) $(CFLAGS) u1_conf.o geometry_st.o random.o pcg32min.o -o $@ -lm

u1_conf_analysis: ../src/u1_conf_analysis.c read_data.o
	$(CC) $(CFLAGS) -c ../src/u1_conf_analysis.c 
	$(CC) $(CFLAGS) u1_conf_analysis.o read_data.o -o $@ -lm


.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o
