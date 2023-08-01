CC=gcc
CFLAGS=-O2 -Wall # -Wextra -Werror -pedantic -Wconversion -std=c99

binaries=harm_osc

all: $(binaries)

harm_osc: ../src/harm_osc.c random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/harm_osc.c 
	$(CC) $(CFLAGS) harm_osc.o pcg32min.o random.o -o $@ -lm

boxmuller.o: ../lib/boxmuller.c ../include/boxmuller.h
	$(CC) $(CFLAGS) -c ../lib/boxmuller.c

random.o: ../lib/random.c ../include/pcg32min.h ../include/random.h
	$(CC) $(CFLAGS) -c ../lib/random.c

pcg32min.o: ../lib/pcg32min.c ../include/pcg32min.h 
	$(CC) $(CFLAGS) -c ../lib/pcg32min.c


.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o
