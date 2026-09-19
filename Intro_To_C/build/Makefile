CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion

binaries=arrays flow_control hello_world command_line1 command_line2 dynamic_memory functions integrator rng struct types types_complex write_read integrator2
all: $(binaries)

arrays: ../src/arrays.c
	$(CC) $(CFLAGS) ../src/arrays.c -o $@ 

hello_world: ../src/hello_world.c
	$(CC) $(CFLAGS) ../src/hello_world.c -o $@ 

command_line1: ../src/command_line1.c
	$(CC) $(CFLAGS) ../src/command_line1.c -o $@ 

command_line2: ../src/command_line2.c
	$(CC) $(CFLAGS) ../src/command_line2.c -o $@ -lm

dynamic_memory: ../src/dynamic_memory.c
	$(CC) $(CFLAGS) ../src/dynamic_memory.c -o $@ 

flow_control: ../src/flow_control.c
	$(CC) $(CFLAGS) ../src/flow_control.c -o $@ 

functions: ../src/functions.c
	$(CC) $(CFLAGS) ../src/functions.c -o $@ -lm 

integrator: ../src/integrator.c
	$(CC) $(CFLAGS) ../src/integrator.c -o $@ -lm

int_lib.o: ../src/int_lib.c ../src/int_lib.h
	$(CC) $(CFLAGS) -c ../src/int_lib.c 

integrator2: ../src/integrator2.c int_lib.o
	$(CC) $(CFLAGS) ../src/integrator2.c int_lib.o -o $@ -lm

rng: ../src/rng.c
	$(CC) $(CFLAGS) ../src/rng.c -o $@ -lm 

struct: ../src/struct.c
	$(CC) $(CFLAGS) ../src/struct.c -o $@  

types: ../src/types.c
	$(CC) $(CFLAGS) ../src/types.c -o $@  

types_complex: ../src/types_complex.c
	$(CC) $(CFLAGS) ../src/types_complex.c -o $@ -std=c99

write_read: ../src/write_read.c
	$(CC) $(CFLAGS) ../src/write_read.c -o $@  


.PHONY: clean
clean:
	rm -f $(binaries) *.o *.dat

cleanobj:
	rm -f *.o
