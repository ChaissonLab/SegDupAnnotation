all: comask bemask


htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --disable-s3 --disable-lzma --disable-bz2 --prefix=$(PWD)/htslib/ && make -j 4 && make install

comask: CombineMask.cpp htslib/lib/libhts.a
	g++ -O2  $< -o $@  -I htslib/include -Lhtslib/lib -lhts -Wl,-rpath,$(PWD)/htslib/lib  -lhts -lz -lpthread


bemask: MaskBed.cpp htslib/lib/libhts.a
	g++ -O2  $< -o $@  -I htslib/include -Lhtslib/lib -lhts -Wl,-rpath,$(PWD)/htslib/lib  -lhts -lz -lpthread
