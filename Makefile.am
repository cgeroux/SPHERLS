CXX=mpic++

bin_PROGRAMS= SPHERLS SPHERLSanal SPHERLSgen

SPHERLS_SOURCES = \
  src/SPHERLS/dataManipulation.cpp \
  src/SPHERLS/dataMonitoring.cpp \
  src/SPHERLS/global.cpp \
  src/SPHERLS/main.cpp \
  src/SPHERLS/physEquations.cpp \
  src/SPHERLS/watchzone.cpp \
  src/SPHERLS/dataManipulation.h \
  src/SPHERLS/dataMonitoring.h \
  src/SPHERLS/global.h \
  src/SPHERLS/main.h \
  src/SPHERLS/physEquations.h \
  src/SPHERLS/userguide.h \
  src/SPHERLS/watchzone.h \
  src/eos.cpp src/eos.h \
  src/exception2.cpp \
  src/exception2.h \
  src/xmlFunctions.cpp \
  src/xmlFunctions.h \
  src/xmlParser.cpp \
  src/xmlParser.h
SPHERLS_LDADD = -L$(PETSC_LIB_PATH) -lpetsc
SPHERLS_CPPFLAGS = -I$(PETSC_INCLUDE_PATH) -Isrc/

SPHERLSanal_SOURCES = \
  src/SPHERLSanal/main.cpp \
  src/SPHERLSanal/main.h \
  src/SPHERLSanal/userguide.h \
  src/eos.h \
  src/eos.cpp \
  src/exception2.cpp \
  src/exception2.h 
SPHERLSanal_LDADD = -L$(FFTW_LIB_PATH) -L$(HDF_LIB_PATH) -lfftw3 -lm -ldf -lmfhdf -ljpeg -lz
SPHERLSanal_CPPFLAGS = -I$(FFTW_INCLUDE_PATH) -I$(HDF_INCLUDE_PATH) -Isrc/

SPHERLSgen_SOURCES = \
  src/SPHERLSgen/main.cpp \
  src/SPHERLSgen/main.h \
  src/SPHERLSgen/userguide.h \
  src/exception2.cpp \
  src/exception2.h \
  src/xmlFunctions.cpp \
  src/xmlFunctions.h \
  src/eos.h \
  src/eos.cpp \
  src/xmlParser.h \
  src/xmlParser.cpp
SPHERLSgen_CPPFLAGS = -Isrc/

.PHONY:docs
docs: docs/manuals/SPHERLS_manual.pdf docs/manuals/SPHERLSgen_manual.pdf \
  docs/manuals/SPHERLSanal_manual.pdf

docs/manuals/SPHERLS_manual.pdf: ${SPHERLS_SOURCES} src/SPHERLS/Doxyfile
	(cd src/SPHERLS; doxygen)
	(cd src/SPHERLS/docs/latex; make pdf;)
	(cd src/SPHERLS/docs/latex;mv refman.pdf ../../../../docs/manuals/SPHERLS_manual.pdf)

docs/manuals/SPHERLSgen_manual.pdf: ${SPHERLSgen_SOURCES} src/SPHERLSgen/Doxyfile
	(cd src/SPHERLSgen; doxygen)
	(cd src/SPHERLSgen/docs/latex; make pdf;)
	(cd src/SPHERLSgen/docs/latex;mv refman.pdf ../../../../docs/manuals/SPHERLSgen_manual.pdf)

docs/manuals/SPHERLSanal_manual.pdf: ${SPHERLSanal_SOURCES} src/SPHERLSanal/Doxyfile
	(cd src/SPHERLSanal; doxygen)
	(cd src/SPHERLSanal/docs/latex; make pdf;)
	(cd src/SPHERLSanal/docs/latex;mv refman.pdf ../../../../docs/manuals/SPHERLSanal_manual.pdf)


