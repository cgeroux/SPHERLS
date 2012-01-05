bin_PROGRAMS = SPHERLSanal
SPHERLSanal_SOURCES = main.cpp main.h userguide.h ../eos.h ../eos.cpp ../exception2.cpp ../exception2.h 
SPHERLSanal_LDADD = -L$(FFTW_LIB_PATH) -L$(HDF_LIB_PATH) -lfftw3 -lm -ldf -lmfhdf -ljpeg -lz
SPHERLSanal_CPPFLAGS = -I$(FFTW_INCLUDE_PATH) -I$(HDF_INCLUDE_PATH) -I../

docs: ${SPHERLSanal_SOURCES} Doxyfile
	doxygen
	(cd ./docs/latex; make pdf;)
	(cd ./docs/latex;mv refman.pdf ../../../../docs/manuals/SPHERLSanal_manual.pdf)
.PHONY:docs

cleandocs:
	rm -r -f ./docs