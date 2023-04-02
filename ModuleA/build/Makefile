CC=gcc
CFLAGS=-O2 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99

all: bad_rng pi_no_rng pi_ss gauss_distr rew_problem gauss_mcmc mcmc_analysis1 mcmc_analysis2 mcmc_analysis3 saw_ss saw_mcmc

bad_rng: ../src/bad_rng.c
	$(CC) $(CFLAGS) ../src/bad_rng.c -o $@ -lm

pi_no_rng: ../src/pi_no_rng.c
	$(CC) $(CFLAGS) ../src/pi_no_rng.c -o $@ -lm

random.o: ../lib/random.c ../include/pcg32min.h ../include/random.h
	$(CC) $(CFLAGS) -c ../lib/random.c

pcg32min.o: ../lib/pcg32min.c ../include/pcg32min.h 
	$(CC) $(CFLAGS) -c ../lib/pcg32min.c

pi_ss: ../src/pi_ss.c random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/pi_ss.c  
	$(CC) $(CFLAGS) pcg32min.o random.o pi_ss.o -o $@ -lm

boxmuller.o: ../lib/boxmuller.c ../include/boxmuller.h
	$(CC) $(CFLAGS) -c ../lib/boxmuller.c

gauss_distr: ../src/gauss_distr.c random.o pcg32min.o boxmuller.o
	$(CC) $(CFLAGS) -c ../src/gauss_distr.c  
	$(CC) $(CFLAGS) boxmuller.o gauss_distr.o pcg32min.o random.o -o $@ -lm

rew_problem: ../src/rew_problem.c random.o pcg32min.o boxmuller.o
	$(CC) $(CFLAGS) -c ../src/rew_problem.c  
	$(CC) $(CFLAGS) boxmuller.o rew_problem.o pcg32min.o random.o -o $@ -lm

gauss_mcmc: ../src/gauss_mcmc.c random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/gauss_mcmc.c  
	$(CC) $(CFLAGS) gauss_mcmc.o pcg32min.o random.o -o $@ -lm

read_data.o: ../lib/read_data.c ../include/read_data.h
	$(CC) $(CFLAGS) -c ../lib/read_data.c

mcmc_analysis1: ../src/mcmc_analysis1.c read_data.o
	$(CC) $(CFLAGS) -c ../src/mcmc_analysis1.c 
	$(CC) $(CFLAGS) mcmc_analysis1.o pcg32min.o random.o read_data.o -o $@ -lm

mcmc_analysis2: ../src/mcmc_analysis2.c pcg32min.o random.o read_data.o
	$(CC) $(CFLAGS) -c ../src/mcmc_analysis2.c  
	$(CC) $(CFLAGS) mcmc_analysis2.o pcg32min.o random.o read_data.o -o $@ -lm

mcmc_analysis3: ../src/mcmc_analysis3.c pcg32min.o random.o read_data.o
	$(CC) $(CFLAGS) -c ../src/mcmc_analysis3.c  
	$(CC) $(CFLAGS) mcmc_analysis3.o pcg32min.o random.o read_data.o -o $@ -lm

saw_ss: ../src/saw_ss.c random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/saw_ss.c  
	$(CC) $(CFLAGS) pcg32min.o random.o saw_ss.o -o $@ -lm

saw_mcmc: ../src/saw_mcmc.c random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/saw_mcmc.c  
	$(CC) $(CFLAGS) pcg32min.o random.o saw_mcmc.o -o $@ -lm

clean:
	rm -f bad_rng pi_no_rng pi_ss gauss_distr rew_problem gauss_mcmc mcmc_analysis1 mcmc_analysis2 mcmc_analysis3 saw_ss saw_mcmc *.o

cleanobj:
	rm -f *.o
