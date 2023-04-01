all: 
	cd ModuleA/build && make

clean:
	cd ModuleA/build && make clean

cleanobj:
	cd ModuleA/build && make cleanobj

dist: 
	make clean
	mkdir MNcodes
	cp -r Makefile README Module* MNcodes
	tar -czvf MNcodes.tar.gz MNcodes
	rm -rf MNcodes
