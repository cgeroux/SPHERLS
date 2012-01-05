bin_PROGRAMS = SPHERLSgen
SPHERLSgen_SOURCES = main.cpp main.h userguide.h ../exception2.cpp ../exception2.h ../xmlFunctions.cpp ../xmlFunctions.h ../eos.h ../eos.cpp ../xmlParser.h ../xmlParser.cpp
SPHERLSgen_CPPFLAGS = -I../

docs: ${SPHERLSgen_SOURCES} Doxyfile
	doxygen
	(cd ./docs/latex; make pdf;)
	(cd ./docs/latex;mv refman.pdf ../../../../docs/manuals/SPHERLSgen_manual.pdf)
.PHONY:docs

cleandocs:
	rm -r -f ./docs