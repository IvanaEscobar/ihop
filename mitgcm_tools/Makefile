## Notes on how to use this Makefile with ihop package
##   ihop is intended to be compiled within mitgcm, for debug purposes this file 
## was created. Now you can compile to core ihop code isolated from mitgcm. To 
## use this Makefile, bellhop.F90 needs to be modified.
##
## 1) Replace these lines:
## MODULE BELLHOP
## ...
## END MODULE BELLHOP
## with 
## PROGRAM BELLHOP
## ...
## END PROGRAM BELLHOP
##
## 2) Remove IHOP_INIT subroutine by deleting/commenting out the lines:
## SUBROUTINE IHOP_INIT
## ...
## END SUBROUTINE IHOP_INIT
##
## 3) Move CONTAINS to the line above SUBROUTINE BellhopCore
##
### Makefile from top level directory ###
export FC=ifort
export F90C=${FC}

export FFLAGS= -g -fPIC -W0 -WB -O3 -funroll-loops

# MITgcm on Sverdrup
FFLAGS+= -convert big_endian -align -traceback -ftz -fp-model precise \
		 -assume byterecl -xCORE-AVX2 -ip
F90FLAGS=${FFLAGS}

# ______________________________________________________________________________
ROOTDIR     = /home/ivana/MITgcm/MITgcm_current
TOOLSDIR    = $(ROOTDIR)/tools

# Dependency generator
MAKEDEPEND = $(TOOLSDIR)/xmakedepend
# C compiler
CC = icc

DEFINES = -DALLOW_USE_MPI -DALWAYS_USE_MPI -DWORDLENGTH=4 -DALLOW_USE_MPI -DHAVE_SYSTEM -DHAVE_FDATE -DHAVE_ETIME_FCT -DHAVE_CLOC -DHAVE_SETRLSTK -DHAVE_SIGREG -DHAVE_STAT -DHAVE_NETCDF -DHAVE_FLUSH 
# Includes for CPP
INCLUDES = -I/opt/ohpc/pub/mpi/openmpi-intel/1.10.7/include -I/opt/ohpc/pub/libs/intel/openmpi/netcdf-fortran/4.4.4/include -I/opt/ohpc/pub/libs/intel/openmpi/netcdf/4.4.1.1/include
# C preprocessing and replacing the _d in constants:
CPPCMD = cat $< |  cpp -traditional -P $(DEFINES) $(INCLUDES) | $(TOOLSDIR)/set64bitConst.sh
# ______________________________________________________________________________

EXECUTABLES = ihop

MISC_OBJECTS 	= 	ihop_fatalerror.o 		attenmod.o 		beampattern.o \
					ihopparams.o 	monotonicmod.o 	pchipmod.o \
					polymod.o 			refcoef.o 		\
					sortmod.o 			srpositions.o \
					splinec.o 			subtabulate.o
BELLHOP_OBJECTS	=	anglemod.o			arrmod.o 		bdrymod.o \
					ihopmod.o 		sspmod.o 		writeray.o \
					influence.o 		step.o 			readenvihop.o \
					bellhop.o	
BELLHOP_LIBS =

all: $(EXECUTABLES)
	@echo "*************"
	@echo "Bellhop built"
	@echo "*************"
	@echo " "

install: all
	for f in $(EXECUTABLES) ; do \
		echo "----- Installing $$f"; cp -p $$f ../bin; \
	done

clean:
	-rm -f bin/*.exe
	find . -name '*.dSYM' -exec rm -r {} +
	find . -name '*.png'  -exec rm -r {} +
	find . -name '*.eps'  -exec rm -r {} +
	find . -name '*.mod'  -exec rm -r {} +
	find . -name '*.grn'  -exec rm -r {} +
	find . -name '*.shd'  -exec rm -r {} +
	find . -name '*.shd.mat'  -exec rm -r {} +
	find . -name '*.prt'  -exec rm -r {} +
	-rm -f *.o *.mod *.exe *_genmod.F90

CLEAN: clean
	find . -name '*.gcno' -exec rm -r {} +
	find . -name '*.gcda' -exec rm -r {} +

# Declare implicit rules (for make utilities with missing or broken build rules)

.SUFFIXES:
.SUFFIXES: .a .o .mod .F90 .f90

.f.o:
	$(FC) -c $(FFLAGS) $< 

.F90.o:
	$(FC) -c $(FFLAGS) $< 

.F90.mod:
	$(FC) -c $(FFLAGS) $< 

# ______________________________________________________________________________
.F.f:
	$(CPPCMD)  > $@
.F.o:
	$(FC) $(FFLAGS) -c $<
.F90.f90:
	$(CPPCMD)  > $@
.FF90.ff90:
	$(CPPCMD)  > $@
# ______________________________________________________________________________

# Declare explicit dependencies (and explicit build rules as needed)

###########################
### from misc directory ###
###########################
ihop_fatalerror.o: ihop_fatalerror.F90 ihopparams.o

attenmod.o: attenmod.F90 ihopparams.o ihop_fatalerror.o
	$(FC) -c $(FFLAGS) $<

beampattern.o: beampattern.F90 ihop_fatalerror.o
	$(FC) -c $(FFLAGS) $<

ihopparams.o: ihopparams.F90
	$(FC) -c $(FFLAGS) $<

monotonicmod.o: monotonicmod.F90
	$(FC) -c $(FFLAGS) $<

pchipmod.o: pchipmod.F90 splinec.o
	$(FC) -c $(FFLAGS) $<

polymod.o: polymod.F90 
	$(FC) -c $(FFLAGS) $<

refcoef.o : refcoef.F90 polymod.o ihop_fatalerror.o
	$(FC) -c $(FFLAGS) $<

sortmod.o: sortmod.F90
	$(FC) -c $(FFLAGS) $<

srpositions.o: srpositions.F90 monotonicmod.o sortmod.o subtabulate.o ihop_fatalerror.o
	$(FC) -c $(FFLAGS) $<

splinec.o: splinec.F90

subtabulate.o: subtabulate.F90
	$(FC) -c $(FFLAGS) $<

##############################
### from bellhop directory ###
##############################
ihop: $(MISC_OBJECTS) $(BELLHOP_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(MISC_OBJECTS) $(BELLHOP_OBJECTS) $(BELLHOP_LIBS)

anglemod.o : anglemod.F90 srpositions.o sortmod.o subtabulate.o ihopparams.o ihop_fatalerror.o
	$(FC) -c $(FFLAGS) $< 

arrmod.o: arrmod.F90 ihopparams.o
	$(FC) -c $(FFLAGS) $< 

bdrymod.o: bdrymod.F90 ihopparams.o ihop_fatalerror.o monotonicmod.o
	$(FC) -c $(FFLAGS) $< 

ihopmod.o: ihopmod.F90 ihopparams.o
	$(FC) -c $(FFLAGS) $< 

influence.o: influence.F90 srpositions.o sspmod.o writeray.o ihopparams.o arrmod.o ihopmod.o
	$(FC) -c $(FFLAGS) $< 

readenvihop.o: readenvihop.F90 sspmod.o attenmod.o srpositions.o anglemod.o ihop_fatalerror.o ihopparams.o ihopmod.o
	$(FC) -c $(FFLAGS) $< 

sspmod.o: sspmod.F90 ihop_fatalerror.o pchipmod.o attenmod.o ihopparams.o splinec.o
	$(FC) -c $(FFLAGS) $< 

step.o: step.F90 sspmod.o ihopparams.o ihopmod.o bdrymod.o
	$(FC) -c $(FFLAGS) $<

writeray.o: writeray.F90 sspmod.o ihopmod.o ihopparams.o
	$(FC) -c $(FFLAGS) $<

bellhop.o:  refcoef.o attenmod.o srpositions.o beampattern.o \
			ihop_fatalerror.o ihopparams.o arrmod.o anglemod.o bdrymod.o \
			readenvihop.o ihopmod.o 	sspmod.o \
			writeray.o 			influence.o   step.o \
			bellhop.F90
