all: 
	cd build && make

clean:
	cd build && make clean

cleanobj:
	cd build && make cleanobj

dist: 
	make clean
	mkdir MNcodesA
	cp -r Makefile README build include lib src MNcodesA
	tar -czvf MNcodesA.tar.gz MNcodesA
	rm -rf MNcodesA
