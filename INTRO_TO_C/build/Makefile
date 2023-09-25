CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion

binaries=hello_world command_line1 command_line2 functions rng dynamic_memory

all: $(binaries)

hello_world: ../src/hello_world.c
	$(CC) $(CFLAGS) ../src/hello_world.c -o $@ 

command_line1: ../src/command_line1.c
	$(CC) $(CFLAGS) ../src/command_line1.c -o $@ 

command_line2: ../src/command_line2.c
	$(CC) $(CFLAGS) ../src/command_line2.c -o $@ -lm

functions: ../src/functions.c
	$(CC) $(CFLAGS) ../src/functions.c -o $@ 

rng: ../src/rng.c
	$(CC) $(CFLAGS) ../src/rng.c -o $@ -lm 

dynamic_memory: ../src/dynamic_memory.c
	$(CC) $(CFLAGS) ../src/dynamic_memory.c -o $@ 


.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o
