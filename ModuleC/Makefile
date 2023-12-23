all: 
	cd build && make

clean:
	cd build && make clean

cleanobj:
	cd build && make cleanobj

dist: 
	make clean
	mkdir MNcodesC
	cp -r Makefile README build include lib src MNcodesC
	tar -czvf MNcodesC.tar.gz MNcodesC
	rm -rf MNcodesC
