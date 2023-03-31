buid:
	cd build && make

clean:
	cd build && make clean

cleanobj:
	cd build && make cleanobj

dist:
	mkdir MN
	cp -r Makefile README build include lib src MN
	tar -czvf MN.tar.gz MN
	rm -rf MN
