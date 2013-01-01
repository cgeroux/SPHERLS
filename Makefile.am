#CXX=mpic++, this should really be the case after running configure

###############################################################################
# EXTRA_DIST ##################################################################
###############################################################################

#energy profiles
EXTRA_DIST=data/energy_pro/Epro.dat
EXTRA_DIST+=data/energy_pro/Epro_mod.dat
EXTRA_DIST+=data/energy_pro/Epro_points.dat

#OPAL/Alexander Furguson EOS and opacity files
EXTRA_DIST+=data/OPAL_EOS_2005/EOS5_0*
EXTRA_DIST+=data/OPAL_OPACITY_GN93/GN93hz
EXTRA_DIST+=data/AF_OPACITY_GN93_1994/g*.*.tron

#eos files
EXTRA_DIST+=data/eos/convert.cpp
EXTRA_DIST+=data/eos/Makefile
EXTRA_DIST+=data/eos/eosCTOPAC
EXTRA_DIST+=data/eos/eosY240Z002
EXTRA_DIST+=data/eos/eosY300Z002
EXTRA_DIST+=data/eos/eosNewY240Z002_wider_finer
EXTRA_DIST+=data/eos/eosNewY300Z002_wider_finer

#ref calcs, will include .svn files :(
EXTRA_DIST+=data/ref_calcs/1DNA/1*
EXTRA_DIST+=data/ref_calcs/1DNA/*.xml
EXTRA_DIST+=data/ref_calcs/1DNA/*.txt
EXTRA_DIST+=data/ref_calcs/2DNA/2*
EXTRA_DIST+=data/ref_calcs/2DNA/*.xml
EXTRA_DIST+=data/ref_calcs/2DNA/*.txt
EXTRA_DIST+=data/ref_calcs/3DNA/3*
EXTRA_DIST+=data/ref_calcs/3DNA/*.xml
EXTRA_DIST+=data/ref_calcs/3DNA/*.txt

#velocity profiles
EXTRA_DIST+=data/velocity_pro/*.dat

#scripts
EXTRA_DIST+=scripts/*.py

#docs
EXTRA_DIST+=docs/manuals/*.pdf
EXTRA_DIST+=docs/templateXML/*.xml

#python c++ extensions
EXTRA_DIST+=src/pythonextensions/hdfwrapper/hdfmodule.cpp

#supporting libraries
EXTRA_DIST+=libs/fftw-3.3.3.tar.gz
EXTRA_DIST+=libs/hdf-4.2.8.tar.gz
EXTRA_DIST+=libs/jpegsrc.v8d.tar.gz
EXTRA_DIST+=libs/openmpi-1.6.3.tar.gz
EXTRA_DIST+=libs/petsc-3.1-p8.tar.gz

###############################################################################
# INSTALL DATA ################################################################
###############################################################################

install-data-local:
	$(MKDIR_P) $(DESTDIR)$(datadir)/$(PACKAGE)/data
	for f in a b; do \
	  $(INSTALL_DATA) $(srcdir)/data/$$f \
	  $(DESTDIR)$(datadir)/$(PACKAGE)/data; done
uninstall-local:
	rm $(DESTDIR)$(datadir)/$(PACKAGE)/data/[ab]
	-rmdir $(DESTDIR)$(datadir)/$(PACKAGE)/data

###############################################################################
# DOCS ########################################################################
###############################################################################

DOC_FILES=./src/userguide.txt
  
###############################################################################
# BINS + SOURCE ###############################################################
###############################################################################

bin_PROGRAMS=	SPHERLS	SPHERLSanal	SPHERLSgen
#dist_pkgdata_DATA=./data/energy_pro/*	./data/eos/*	./data/ref_cals/*	./data/velocity_pro/*
SPHERLS_SOURCES	=	\
	src/SPHERLS/dataManipulation.cpp	\
	src/SPHERLS/dataMonitoring.cpp	\
	src/SPHERLS/dataMonitoring.h	\
	src/SPHERLS/global.cpp	\
	src/SPHERLS/global.h	\
	src/SPHERLS/main.cpp	\
	src/SPHERLS/main.h	\
	src/SPHERLS/physEquations.cpp	\
	src/SPHERLS/physEquations.h	\
	src/SPHERLS/watchzone.cpp	\
	src/SPHERLS/dataManipulation.h	\
	src/SPHERLS/time.h	\
	src/SPHERLS/time.cpp	\
	src/SPHERLS/procTop.h	\
	src/SPHERLS/procTop.cpp	\
	src/SPHERLS/watchzone.h	\
	src/SPHERLS/profileData.cpp	\
	src/SPHERLS/profileData.h	\
	src/SPHERLS/fileExists.h	\
	src/SPHERLS/fileExists.cpp	\
	src/eos.cpp	\
	src/eos.h	\
	src/exception2.cpp	\
	src/exception2.h	\
	src/xmlFunctions.cpp	\
	src/xmlFunctions.h	\
	src/xmlParser.cpp	\
	src/xmlParser.h	
SPHERLS_LDADD	=	-L$(PETSC_LIB_PATH)	-lpetsc
SPHERLS_CPPFLAGS	=	-I$(PETSC_INCLUDE_PATH)	-Isrc/

SPHERLSanal_SOURCES	=	\
	src/SPHERLSanal/main.cpp	\
	src/SPHERLSanal/main.h	\
	src/eos.h	\
	src/eos.cpp	\
	src/exception2.cpp	\
	src/exception2.h	

if FFTW_ENABLE
FFTW_LDADDS=-L$(FFTW_LIB_PATH) -lfftw3
FFTW_CPPFLAGS=-I$(FFTW_INCLUDE_PATH)
else
FFTW_LDADDS= 
FFTW_CPPFLAGS= 
endif
if HDF_ENABLE
HDF_LDADDS=-L$(HDF_LIB_PATH) -lm -lmfhdf -ldf -ljpeg -lz
HDF_CPPFLAGS=-I$(HDF_INCLUDE_PATH)
BUILT_SOURCES=hdf.so
else
HDF_LDADDS= 
HDF_CPPFLAGS= 
endif

SPHERLSanal_LDADD=$(FFTW_LDADDS) $(HDF_LDADDS) 
SPHERLSanal_CPPFLAGS=$(FFTW_CPPFLAGS) $(HDF_CPPFLAGS) -Isrc/ 

SPHERLSgen_SOURCES	=	\
	src/SPHERLSgen/main.cpp	\
	src/SPHERLSgen/main.h	\
	src/exception2.cpp	\
	src/exception2.h	\
	src/xmlFunctions.cpp	\
	src/xmlFunctions.h	\
	src/eos.h	\
	src/eos.cpp	\
	src/xmlParser.h	\
	src/xmlParser.cpp
SPHERLSgen_CPPFLAGS	=	-Isrc/

.PHONY:docs
docs:	docs/manuals/SPHERLS_Ref_Manual.pdf	docs/manuals/SPHERLSgen_Ref_Manual.pdf	\
	docs/manuals/SPHERLSanal_Ref_Manual.pdf	docs/manuals/Scripts_Ref_Manual.pdf docs/manuals/SPHERLS_User_Guide.pdf

if HDF_ENABLE
hdf.so:
	(cd	src/pythonextensions/hdfwrapper;	python setup.py build;python setup.py install --home=../)
endif

docs/manuals/SPHERLS_Ref_Manual.pdf:	${DOC_FILES}	${SPHERLS_SOURCES}	src/SPHERLS/Doxyfile
	(cd	src/SPHERLS;echo	${PWD};	doxygen)
	(cd	src/SPHERLS/docs/latex;	make	pdf;)
	(cd	src/SPHERLS/docs/latex;mv	refman.pdf	../../../../docs/manuals/SPHERLS_Ref_Manual.pdf)

docs/manuals/SPHERLSgen_Ref_Manual.pdf:	${SPHERLSgen_SOURCES}	src/SPHERLSgen/Doxyfile
	(cd	src/SPHERLSgen;	doxygen)
	(cd	src/SPHERLSgen/docs/latex;	make	pdf;)
	(cd	src/SPHERLSgen/docs/latex;mv	refman.pdf	../../../../docs/manuals/SPHERLSgen_Ref_Manual.pdf)

docs/manuals/SPHERLSanal_Ref_Manual.pdf:	${SPHERLSanal_SOURCES}	src/SPHERLSanal/Doxyfile
	(cd	src/SPHERLSanal;	doxygen)
	(cd	src/SPHERLSanal/docs/latex;	make	pdf;)
	(cd	src/SPHERLSanal/docs/latex;mv	refman.pdf	../../../../docs/manuals/SPHERLSanal_Ref_Manual.pdf)

docs/manuals/Scripts_Ref_Manual.pdf:	scripts/Doxyfile	scripts/*.py
	(cd	scripts;	doxygen)
	(cd	scripts/docs/latex;	make	pdf;)
	(cd	scripts/docs/latex;mv	refman.pdf	../../../docs/manuals/Scripts_Ref_Manual.pdf)
  
docs/manuals/SPHERLS_User_Guide.pdf: src/Doxyfile ${DOC_FILES}
	(cd	src;	doxygen)
	(cd	src/docs/latex;	make	pdf;)
	(cd	src/docs/latex;mv	refman.pdf	../../../docs/manuals/SPHERLS_User_Guide.pdf)