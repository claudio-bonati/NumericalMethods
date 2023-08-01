all: 
	cd build && make

clean:
	cd build && make clean

cleanobj:
	cd build && make cleanobj

dist: 
	make clean
	mkdir MNcodesB
	cp -r Makefile README build include lib src MNcodesB
	tar -czvf MNcodesB.tar.gz MNcodesB
	rm -rf MNcodesB
