all: comask bemask mask2bed nl toupper sedef/sedef


htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --disable-s3 --disable-lzma --disable-bz2 --prefix=$(PWD)/htslib/ && make -j 4 && make install

sedef/sedef:
	cd sedef && make -j 4


toupper: ToUpper.cpp htslib/lib/libhts.a
	g++ -O2  $< -o $@  -I htslib/include -Lhtslib/lib -lhts -Wl,-rpath,$(PWD)/htslib/lib  -lhts -lz -lpthread


comask: CombineMask.cpp htslib/lib/libhts.a
	g++ -O2  $< -o $@  -I htslib/include -Lhtslib/lib -lhts -Wl,-rpath,$(PWD)/htslib/lib  -lhts -lz -lpthread

bemask: MaskBed.cpp htslib/lib/libhts.a
	g++ -O2  $< -o $@  -I htslib/include -Lhtslib/lib -lhts -Wl,-rpath,$(PWD)/htslib/lib  -lhts -lz -lpthread

mask2bed: MaskedToBed.cpp htslib/lib/libhts.a
	g++ -O2  $< -o $@  -I htslib/include -Lhtslib/lib -lhts -Wl,-rpath,$(PWD)/htslib/lib  -lhts -lz -lpthread

countn: CountN.cpp
	g++ -O2 CountN.cpp -o countn

nl: CountRep.cpp
	g++ -O2 CountRep.cpp -o nl
