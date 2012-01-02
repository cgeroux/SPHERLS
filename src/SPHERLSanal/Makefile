######################
#### MPI MAKEFILE ####
######################

#program name
TARGET = SPHERLSanal

#files needed to make target other than exercise1.o
LIBDIR=/home/cgeroux/lib
EXTOBJS =${LIBDIR}/exception/exception2.o\
  ${LIBDIR}/fftw/fftw++.o\
  /home/cgeroux/SPHERLS/trunk/eos.o

INTOBJS=\
  main.o

INTSOURCE=main.cpp

INTHEADERS=\
  main.h\
  userguide.h
SCRIPTPATH=./scripts

PYTHONSCRIPTS=\
  ${SCRIPTPATH}/cp_files.py\
  ${SCRIPTPATH}/mv_files.py\
  ${SCRIPTPATH}/rm_files.py
  

#macros to set compiler name/options and file naming conventions
COMPILER  = g++
LINKER    = g++
LINKFLAGS =  -fsignaling-nans -Wall -g -o
COMPBEFFLAGS = -fsignaling-nans -Wall -g -c
COMPAFTFLAGS = -o
OBJECTEXT = .o
SOURCEEXT = .cpp
HEADEREXT = .h
LIBS      = -L${LIBDIR}/fftw-3.3/lib -lfftw3
INCLUDE   = -I${LIBDIR}/exception -I${LIBDIR}/fftw -I${LIBDIR}/fftw-3.3/include -I/home/cgeroux/SPHERLS/trunk

##################
#### make bin ####
##################

#link program to make binary
${TARGET}: ${EXTOBJS} ${INTOBJS} ${INTHEADERS}
	${LINKER} ${LINKFLAGS} ${TARGET} ${EXTOBJS} ${INTOBJS} ${LIBS}

#compiles source files, but doesn't link
%${OBJECTEXT}:%${SOURCEEXT}
	${COMPILER} ${COMPBEFFLAGS} $< ${COMPAFTFLAGS} $@ ${INCLUDE}

####################
##### clean up #####
####################

#cleans up object files, and executable
clean:cleano cleanout cleanx cleansh cleancore

#cleans up object files
cleano:
  #remove object files if they exsist
	@(for i in ${OBJS} ${TARGET}.o; do\
	  if [ -f $$i ]; then rm $$i;echo rm $$i;\
	  fi;\
	  done;)

#cleans up output files of run
cleanout:
  #remove outputfiles
	@(for i in ${TARGET}*.out; do\
	  if [ -f $$i ]; then rm $$i;echo rm $$i;\
	  fi;\
	  done;)
	@(if [ -f ${TARGET}.err ]; then rm ${TARGET}.err;echo rm ${TARGET}.err;fi)

#cleans up executable
cleanx:
  #remove executable
	@(if [ -f ${TARGET} ]; then rm ${TARGET};echo rm ${TARGET};fi)

#cleans up submit script
cleansh:
  #remove submit scripts
	@(if [ -f ${TARGET}_test_que.sh ]; then rm ${TARGET}_test_que.sh;echo rm ${TARGET}_test_que.sh;fi)
	@(if [ -f ${TARGET}.sh ]; then rm ${TARGET}.sh;echo rm ${TARGET}.sh;fi)

cleancore:
#cleans up any core dumps
	@(for i in core*; do\
	  if [ -f $$i ]; then rm $$i;echo rm $$i;\
	  fi;\
	  done;)

#delete documentation
cleandocs:
	rm -r -f ./docs

#################
#### run job ####
#################

#production job
run:submit_script
	@echo "not yet implimented"

#runs on test que
runt:submit_script_testq
	@qsub ${TARGET}_test_que.sh

#run on current machine
rund:
	mpirun -np ${NUMPROCS} ./${TARGET}

#run mpi code with portaland group debugger
rundd:
	mpirun --debug -np ${NUMPROCS} ./${TARGET}
#	pgdbg -text -mpi:mpirun -np ${NUMPROCS} ./${TARGET}

#creates the submit script with correct names
submit_script_testq:
	@echo '#!/bin/sh'           >${TARGET}_test_que.sh
	@echo '#$$ -S /bin/sh'      >>${TARGET}_test_que.sh
	@echo '#$$ -N ${TARGET}'>>${TARGET}_test_que.sh
	@echo '#$$ -cwd'            >>${TARGET}_test_que.sh
	@echo '## output file name'            >>${TARGET}_test_que.sh
	@echo '#$$ -o ${TARGET}.out'>>${TARGET}_test_que.sh
	@echo '## error output file'            >>${TARGET}_test_que.sh
	@echo '#$$ -e ${TARGET}.err'>>${TARGET}_test_que.sh
	@echo '## error is not merged with standard out'            >>${TARGET}_test_que.sh
	@echo '#$$ -j n'            >>${TARGET}_test_que.sh
	@echo '#$$ -pe mpich ${NUMPROCS}'>>${TARGET}_test_que.sh
	@echo '## below is the run time in seconds needed'>>${TARGET}_test_que.sh
	@echo '#$$ -l h_rt=3600'>>${TARGET}_test_que.sh
	@echo '## submit to the test que'>>${TARGET}_test_que.sh
	@echo '#$$ -l test=true'>>${TARGET}_test_que.sh
	@echo '#$$ -l h_vmem=1G'>>${TARGET}_test_que.sh
	@echo '## transfer all environment variables when job is submitted'>>${TARGET}_test_que.sh
	@echo '#$$ -V'>>${TARGET}_test_que.sh
	@echo '#$$ -l h_stack=10240K'>>${TARGET}_test_que.sh
	@echo '## eamil on end, abort, suspend'>>${TARGET}_test_que.sh
	@echo '#$$ -m eas'>>${TARGET}_test_que.sh
	@echo '## who to email'>>${TARGET}_test_que.sh
	@echo '#$$ -M cgeroux@ap.smu.ca'>>${TARGET}_test_que.sh
#	@echo 'hostname'>>${TARGET}.sh
	@echo 'mpirun -np $$NSLOTS -machinefile $$TMPDIR/machines ./${TARGET}'>>${TARGET}_test_que.sh
#use my own machine file to make sure I always get the same node
#	@echo 'mpirun -np $$NSLOTS -machinefile ./machines ./${TARGET}'>>${TARGET}.sh
#	@echo 'date'>>${TARGET}.sh

submit_script:
	@echo "production que script code not written, can not submit"

#make documentation
docs:${INTSOURCE} Doxyfile ${INTHEADERS} ${PYTHONSCRIPTS}
	@echo ${INTSOURCE} Doxyfile ${INTHEADERS} ${PYTHONSCRIPTS}
	@doxygen
	@(cd ./docs/latex; make pdf)