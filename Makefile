all: 
	cd build && make

clean:
	cd build && make clean

cleanobj:
	cd build && make cleanobj

dist: 
	make clean
	mkdir MNcodes
	cp -r Makefile README build include lib src MNcodes
	tar -czvf MNcodes.tar.gz MNcodes
	rm -rf MNcodes
