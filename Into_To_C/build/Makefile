CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion

binaries=arrays flow_control hello_world command_line1 command_line2 dynamic_memory functions integrator rng types write_read
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
	$(CC) $(CFLAGS) ../src/functions.c -o $@ 

integrator: ../src/integrator.c
	$(CC) $(CFLAGS) ../src/integrator.c -o $@ -lm

rng: ../src/rng.c
	$(CC) $(CFLAGS) ../src/rng.c -o $@ -lm 

types: ../src/types.c
	$(CC) $(CFLAGS) ../src/types.c -o $@  

write_read: ../src/write_read.c
	$(CC) $(CFLAGS) ../src/write_read.c -o $@  


.PHONY: clean
clean:
	rm -f $(binaries) *.o *.dat

cleanobj:
	rm -f *.o
