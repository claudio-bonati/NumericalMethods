all: bad_rng

bad_rng: ./src/bad_rng.c
	gcc -O2 ./src/bad_rng.c -o bad_rng -lm

clean:
	rm -f bad_rng
