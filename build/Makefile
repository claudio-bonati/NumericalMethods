CC=gcc
CFLAGS=-O2 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99

all: bad_rng pi_no_rng pi_ss gauss_distr

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
	$(CC) $(CFLAGS) boxmuller.o pcg32min.o random.o gauss_distr.o -o $@ -lm

clean:
	rm -f bad_rng pi_no_rng pi_ss gauss_distr *.o

cleanobj:
	rm -f *.o
