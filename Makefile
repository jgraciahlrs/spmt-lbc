#-------------------------------------------------------------------------
#									
#               M A K E F I L E   F O R   L B C		        
#									
#-------------------------------------------------------------------------

# compile on  Asama/A1/ONTAKE for SX 
# --------------------------------
#PLATFORM?=$(shell uname -m)
#HOSTNAME?=$(shell uname -n|cut --delimiter=. --fields=1)
#

#---------------------
# Define Machine Name here:

HPC_TYPE=nehalem
#HPC_TYPE=cray     # uses Cray compiler wrappers

#---------------------
# Compiler:

COMPILER=gnu_openmp#_profile#_openmp#_debug#
#COMPILER=intel#_profile#_debug#
#COMPILER=intel




#-------------------------------------------------------------------------
# Possible Arguments for Parallel

# (arguments are set using -DARG1 -DARG2

#-----------------------------------------------------
# arg                             # description

# For parallel either use USE_MPI with suboptions
PARALLEL += -DUSE_MPI    # Activate MPI:              
#PARALLEL += -DDEBUG_MPI 
PARALLEL += -DVERBOSE   # Output detailed information about parallel execution 
# MPI Method
#PARALLEL += -DPERSISTENT   # Persistent MPI Communication. Must be used together with BUFFER_EACH 
#PARALLEL += -DPERSISTENT2   # Persistent MPI Communication. Must be used together with BUFFER_EACH 
#PARALLEL += -DISEND_IRECV       # use MPI_IRECV and MPI_ISEND
PARALLEL += -DSENDRECV          # use MPI_SENDRECV

# MPI options
#PARALLEL += -DMPI_SUBARRAY      # use a mpi derived type subarray for the exchange
#PARALLEL += -DBUFFER_EACH       # Use separate buffers for each direction 
PARALLEL += -DCOMM_REDUCED      # Only Exchange necessary densities
#PARALLEL += -DSYNCHRONIZE_COMM  # mpl_barrier before communication       

# StarSs options
#PARALLEL += -DUSE_SMPSS -DUSE_TILES -DTILESIZE=32

# OpenMP options; choose either one
PARALLEL += -DUSE_OMP_FINE	# fine-grained loop level OpenMP version
#PARALLEL += -DUSE_OMP_COARSE -DUSE_TILES -DTILESIZE=32 # coarse-grained tile level OpenMP version, shares code with StarSs version
 
# OmpSs parallelisation strategy, if no strategy is set the default(comm_hiding strategy) is used
#PARALLEL += -DUSE_OMPSS_FORK_JOIN

# Compiles your application with ompss MPI
#PARALLEL += -DUSE_OMPSS_MPI

#-----------------------------------------------------
# Model arguments 
# Lattice Boltzmann Model: 
# 3D:  use three-dimensional model with nineteen velocities
# Standard is D3Q19

# Relaxation modes (optional)
# Standard is BGK

# Use the incompressible model for higher stability
#LBMODEL += -DINCOMPRESSIBLE

# Performance issues (optional)
# Combine stream and collide routine, split relaxation 
LBMODEL += -DCOMBINE_STREAM_COLLIDE
#LBMODEL += -DSPLITLOOPS

# Memory layout:
# Choose ONE of the following
# save array in memory as x-,y-,z- component and l as density (1...9 / 1...19)
#LBMODEL += -DLAYOUT_IJKL
#  save array in memory as l as density (1...9 / 1...19) and x-,y-,z- component 
LBMODEL += -DLAYOUT_LIJK

# Use Fortran 2003 statements
#LBMODEL += -DF2003_ALLOC_DERIVED_TYPE # use allocatables in derived type
#LBMODEL += -DVERBOSE # write out extra information 


CPPDEFINES=$(LBMODEL) $(PARALLEL)

#-----------------------------------------------------

#----------------------------------------------------------------
# Machine specific 
#----------------------------------------------------------------
 

#--------------------------------
# Cray machine
ifeq ($(findstring cray,$(HPC_TYPE)),cray)
   COMP=ftn
   CCOMP=cc
endif


#--------------------------------
# INTEL 
ifeq ($(findstring intel,$(COMPILER)),intel)
   COMP?=mpif90
   CCOMP?=mpicc

   F90_CFLAGS_DEBUG=-fpp -O0 -g -check all -r8 -i4 -traceback -implicitnone -warn all -keep
   F90_LFLAGS_DEBUG=-g -O0 -r8 -i4  -L$(MPI_LIB_DIR) -lmpi_f90  
   #F90_CFLAGS_DEBUG+=-warn alignments,declarations,errors,general,ignore_loc,interfaces,truncated_source,uncalled,nounused,usage

   F90_CFLAGS_OPT=-fpp -r8 -i4 -O3 -i-dynamic 
   F90_LFLAGS_OPT=-r8 -i4 -O3 -i-dynamic   -L$(MPI_LIB_DIR) -lmpi_f90
   #--------------------------------
   # PROFILING 
   ifeq ($(findstring profile,$(COMPILER)),profile)
      F90_CFLAGS_DEBUG+=-p -g
      F90_LFLAGS_DEBUG+=-p
      F90_CFLAGS_OPT+=-p -g
      F90_LFLAGS_OPT+=-p
   endif
endif

#--------------------------------
# GNU

ifeq ($(findstring gnu,$(COMPILER)),gnu)
   COMP?=mpif90 -O2
   CCOMP?=mpicc -O2 

   F90_CFLAGS_DEBUG+=-xf95-cpp-input -O0 -K -g3 -fcheck=all -fcheck=bounds -debug -Wall -fbacktrace -fbounds-check -Wall -pedantic -fcray-pointer -fwhole-file
   F90_LFLAGS_DEBUG+=-O0 -g3 -fbacktrace -K  -fcheck=all -Wall -fcheck=bounds -debug

   #F90_CFLAGS_OPT=-xf95-cpp-input -O3 -fcray-pointer -fwhole-file 
   #-std=f2003
   F90_CFLAGS_OPT+=-O3
   F90_LFLAGS_OPT+=
   #--------------------------------
   # OpenMP
   ifeq ($(findstring openmp,$(COMPILER)),openmp)
      F90_CFLAGS_DEBUG += -fopenmp
      F90_LFLAGS_DEBUG += -fopenmp
      F90_CFLAGS_OPT   += -fopenmp
      F90_LFLAGS_OPT   += -fopenmp
   endif 	
   # PROFILING 
   ifeq ($(findstring profile,$(COMPILER)),profile)
      F90_CFLAGS_DEBUG+=-pg
      F90_LFLAGS_DEBUG+=-pg
      F90_CFLAGS_OPT+=-pg
      F90_LFLAGS_OPT+=-pg
   endif
endif




ifeq ($(findstring debug,$(COMPILER)),debug)
F90CFLAGS= $(F90_CFLAGS_DEBUG) 
F90LFLAGS= $(F90_LFLAGS_DEBUG) 
else
F90CFLAGS= $(F90_CFLAGS_OPT) 
F90LFLAGS= $(F90_LFLAGS_OPT) 
endif

#--------------------------------
# share compiler and flag variables with sub-Makefiles

export COMP      
export F90CFLAGS 

#----------------------------------------------------------------
#


#--------------------------------
# links and includes 

BASELIB_DIR = ./lib
INCLUDE_DIR = ./include

LINK_BASELIB=-L$(BASELIB_DIR) 


LIBS= $(LINK_BASELIB)
ifeq ($(findstring USE_OMPSS_MPI,$(PARALLEL)),USE_OMPSS_MPI)
   INCLUDES+=-I$(OMPSS_MPI_HOME)
   LIBS+=-L$(OMPSS_MPI_HOME) -lompssmpi
endif

# Shell:
# ------
SHELL=/bin/bash

sources := nrt_lib.F90 timing.F90 lbmodel.F90 tools.F90 mpl_set.F90
sources += function_set.F90 lb_geo.F90 lb_init.F90 lb_bc.F90 lbm_functions.F90 
sources += lbm_step.F90 lbm_step_tiled.F90

# Object-Files:
# -------------
OBJS := $(sources:.F90=.o)

main := lbc.F90
mainobject := $(main:.F90=.o) 

CMD := lbc

# Suffixes:
# ---------
#.SUFFIXES:
#.SUFFIXES:.o .c .f .f90 .a

# Defult rules:
# -------------
#.f90.o:
#	$(COMP) $(F90CFLAGS) $(CPPDEFINES) $(INCLUDES) -c $<	
#.c.o:
#	$(CCOMP) $(CFLAGS) $(CPPDEFINES) $(INCLUDES) -c $<	
# sparate preprocessor:
# ---------------------
%.o: %.f90
	$(COMP) $(CPPDEFINES) $(F90CFLAGS) -c $< $(INCLUDES) 
%.o: %.F90
	$(COMP) $(CPPDEFINES) $(F90CFLAGS) -c $< $(INCLUDES) 

# dummy rule to override implicit rule in GNU make; sigh!
%.o : %.mod
	;

# Targets:
# --------
#echo:
#	@echo $(OBJS)
all: $(CMD) 


$(CMD): $(OBJS) $(mainobject)
	$(COMP) $(F90LFLAGS) -o  $@  $^ $(INCLUDES) $(LIBS)
	@echo -e '\033[32m                                                       '
	@echo '      '
	@echo '         oooo   .o8                       '
	@echo '          888   888                       '
	@echo '          888   888oooo.   .ooooo.        '
	@echo '          888   d88   88b d88    Y8       '
	@echo '          888   888   888 888             '
	@echo '          888   888   888 888   .o8       '
	@echo '         o888o  `Y8bod8P  `Y8bod8P        '
	@echo '      '
	@echo -e '\033[0m                                                        '
	@echo '==============================================================='
	@echo '=  L B C   C O M P I L E D                                    ='
	@echo '==============================================================='

.PHONY: $(LIBS)

clean:
	- rm -f *.o *.mod *.i lbc *.mf95ss *-adapters.c *-registration.c 
	@echo '==============================================================='
	@echo '=  c l e a n e d                                              ='
	@echo '==============================================================='

cleanall:
	- rm -f *.o *.mod *.i *.i90 lbc *.mf95ss *-registration.c *-adapters.c *.tasks
	@echo '==============================================================='
	@echo '=  c l e a n e d     u p    a l l                             ='
	@echo '==============================================================='

cleandata:
	- rm -f *.plt *.res *.inf fort.* 
	@echo '==============================================================='
	@echo '=  c l e a n e d     d a t a                                  ='
	@echo '==============================================================='

touch:
	- touch *.F90
#display:
#	@echo
#"==============================================================================="
#	@echo "Plattform:        $(PLATFORM)"
#	@echo "HPC_TYPE:         $(HPC_TYPE)"
#	@echo "Fortran-Compiler: $(COMP)"
#	@echo "Fortran-Optionen: $(F90CFLAGS)"
#	@echo "C-Compiler:       $(CCOMP)"
#	@echo "Lib-Extension:    $(LIBEXT)"
#	@echo "BASELIB_DIR:      $(BASELIB_DIR)"
#	@echo
#"================================================================================"



# Abhaengigkeiten
(lbmodel.o): (nrt_lib.o)
(tools.o): (nrt_lib.o) (lbmodel.o) 
(mpl_set.o): (nrt_lib.o) (lbmodel.o) (tools.o) 
(function_set.o): (nrt_lib.o) (timing.o) (lbmodel.o) (mpl_set.o) 
(init.o): (nrt_lib.o) (timing.o) (lbmodel.o) (mpl_set.o) (lb_geo.o)
(lb_bc.o): (nrt_lib.o) (timing.o) (lbmodel.o) (mpl_set.o) (lb_geo.o)
(lbm_functions.o): (nrt_lib.o) (timing.o) (lbmodel.o) (mpl_set.o) (lb_geo.o) (lb_bc.o)
(lbc.o): (nrt_lib.o) (timing.o) (lbmodel.o) (mpl_set.o) (lbm_functions.o) (function_set.o) (lbm_step.o) (lbm_step_tiled.o)
