all: bad_rng pi_no_rng

bad_rng: ./src/bad_rng.c
	gcc -O2 ./src/bad_rng.c -o bad_rng -lm

pi_no_rng: ./src/pi_no_rng.c
	gcc -O2 ./src/pi_no_rng.c -o pi_no_rng -lm

clean:
	rm -f bad_rng pi_no_rng
