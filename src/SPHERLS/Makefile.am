CXX=mpic++
bin_PROGRAMS = SPHERLS
SPHERLS_SOURCES = dataManipulation.cpp dataMonitoring.cpp global.cpp main.cpp physEquations.cpp watchzone.cpp dataManipulation.h dataMonitoring.h global.h main.h physEquations.h userguide.h watchzone.h ../eos.cpp ../eos.h ../exception2.cpp ../exception2.h ../xmlFunctions.cpp ../xmlFunctions.h ../xmlParser.cpp ../xmlParser.h
SPHERLS_LDADD = -L$(PETSC_LIB_PATH) -lpetsc
SPHERLS_CPPFLAGS = -I$(PETSC_INCLUDE_PATH) -I../

docs: ${SPHERLS_SOURCES} Doxyfile
	doxygen
	(cd ./docs/latex; make pdf;)
	(cd ./docs/latex;mv refman.pdf ../../../../docs/manuals/SPHERLS_manual.pdf)
.PHONY:docs

cleandocs:
	rm -r -f ./docs
