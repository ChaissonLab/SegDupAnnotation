all: comask bemask mask2bed nl toupper sedef/sedef hmcnc/src/hmmcnc bamToFreq


sedef/sedef:
	cd sedef && make -j 4

hmcnc/src/hmmcnc: hmcnc/src/hmmcnc.cpp
	cd hmcnc/src && make

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

bamToFreq: BamToFreq.cpp $(CONDA_PREFIX)/lib/libhts.so
	g++ -O2 $< -o $@ -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib -lhts -lpthread -lz -Wl,-rpath,$(CONDA_PREFIX)/lib
