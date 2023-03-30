CC=gcc
CFLAGS=-O2

all: bad_rng pi_no_rng pi_rng gauss_rng

bad_rng: ../src/bad_rng.c
	$(CC) $(CFLAGS) ../src/bad_rng.c -o $@ -lm

pi_no_rng: ../src/pi_no_rng.c
	$(CC) $(CFLAGS) ../src/pi_no_rng.c -o $@ -lm

random.o: ../lib/random.c ../include/pcg32min.h ../include/random.h
	$(CC) $(CFLAGS) -c ../lib/random.c

pcg32min.o: ../lib/pcg32min.c ../include/pcg32min.h 
	$(CC) $(CFLAGS) -c ../lib/pcg32min.c

pi_rng: ../src/pi_rng.c random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/pi_rng.c  
	$(CC) $(CFLAGS) pcg32min.o random.o pi_rng.o -o $@ -lm

boxmuller.o: ../lib/boxmuller.c ../include/boxmuller.h
	$(CC) $(CFLAGS) -c ../lib/boxmuller.c

gauss_rng: ../src/gauss_rng.c random.o pcg32min.o boxmuller.o
	$(CC) $(CFLAGS) -c ../src/gauss_rng.c  
	$(CC) $(CFLAGS) boxmuller.o pcg32min.o random.o gauss_rng.o -o $@ -lm

clean:
	rm -f bad_rng pi_no_rng pi_rng gauss_rng *.o

cleanobj:
	rm -f *.o
