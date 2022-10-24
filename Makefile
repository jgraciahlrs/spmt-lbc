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

#HPC_TYPE=sx8
#HPC_TYPE=sx9
#HPC_TYPE=xe6
HPC_TYPE=nehalem

#---------------------
# Compiler:

# Add _debug for debugging and smpss for StarSs/SMPSs
#COMPILER=gnu_debug
COMPILER=gnu#_ompss#_smpss#_omp#_trace#_debug
#COMPILER=sx
#COMPILER=vampir
#COMPILER=intel_debug
#COMPILER=intel#_smpss




#-------------------------------------------------------------------------
# Possible Arguments for Parallel

# (arguments are set using -DARG1 -DARG2

#-----------------------------------------------------
# arg                             # description

# For parallel either use USE_MPI with suboptions
PARALLEL += $(FDEF)-DUSE_MPI    # Activate MPI:              
#PARALLEL += $(FDEF)-DDEBUG_MPI 
PARALLEL += $(FDEF)-DVERBOSE   # Output detailed information about parallel execution 
# MPI Method
#PARALLEL += $(FDEF)-DPERSISTENT   # Persistent MPI Communication. Must be used together with BUFFER_EACH 
#PARALLEL += $(FDEF)-DPERSISTENT2   # Persistent MPI Communication. Must be used together with BUFFER_EACH 
#PARALLEL += $(FDEF)-DISEND_IRECV       # use MPI_IRECV and MPI_ISEND
PARALLEL += $(FDEF)-DSENDRECV          # use MPI_SENDRECV

# MPI options
#PARALLEL += $(FDEF)-DMPI_SUBARRAY      # use a mpi derived type subarray for the exchange
#PARALLEL += $(FDEF)-DBUFFER_EACH       # Use separate buffers for each direction 
PARALLEL += $(FDEF)-DCOMM_REDUCED      # Only Exchange necessary densities
#PARALLEL += $(FDEF)-DSYNCHRONIZE_COMM  # mpl_barrier before communication       

# StarSs options
PARALLEL += $(FDEF)-DUSE_SMPSS $(FDEF)-DUSE_TILES $(FDEF)-DTILESIZE=32

# OpenMP options; choose either one
#PARALLEL += $(FDEF)-DUSE_OMP_FINE	# fine-grained loop level OpenMP version
#PARALLEL += $(FDEF)-DUSE_OMP_COARSE $(FDEF)-DUSE_TILES $(FDEF)-DTILESIZE=32 # coarse-grained tile level OpenMP version, shares code with StarSs version
 
# OmpSs parallelisation strategy, if no strategy is set the default(comm_hiding strategy) is used
PARALLEL += $(FDEF)-DUSE_OMPSS_FORK_JOIN

# Compiles your application with ompss MPI
PARALLEL += $(FDEF)-DUSE_OMPSS_MPI

#-----------------------------------------------------
# Model arguments 
# Lattice Boltzmann Model: 
# 3D:  use three-dimensional model with nineteen velocities
# Standard is D3Q19

# Relaxation modes (optional)
# Standard is BGK

# Use the incompressible model for higher stability
#LBMODEL += $(FDEF)-DINCOMPRESSIBLE

# Performance issues (optional)
# Combine stream and collide routine, split relaxation 
LBMODEL += $(FDEF)-DCOMBINE_STREAM_COLLIDE
#LBMODEL += $(FDEF)-DSPLITLOOPS

# Memory layout:
# Choose ONE of the following
# save array in memory as x-,y-,z- component and l as density (1...9 / 1...19)
#LBMODEL += $(FDEF)-DLAYOUT_IJKL
#  save array in memory as l as density (1...9 / 1...19) and x-,y-,z- component 
LBMODEL += $(FDEF)-DLAYOUT_LIJK

# Use Fortran 2003 statements
#LBMODEL += $(FDEF)-DF2003_ALLOC_DERIVED_TYPE # use allocatables in derived type
#LBMODEL += $(FDEF)-DVERBOSE # write out extra information 


CPPDEFINES=$(LBMODEL) $(PARALLEL)

#-----------------------------------------------------

#----------------------------------------------------------------
# Machine specific 
#----------------------------------------------------------------

#--------------------------------
# Vector Machine

ifeq ($(findstring sx,$(HPC_TYPE)),sx) #vector machine
   COMP=sxmpif90
   CCOMP=sxmpicc

   F90_CFLAGS_DEBUG=-Chopt  -g -EP  -Wf,'-M zdiv,flovf,flunf,fxovf,inv' -Wf,'-init stack=nan,heap=nan'
   F90_LFLAGS_DEBUG=-g -Cdebug

   #F90CFLAGS_OPT=-float0 -f4 -Pmulti -ew -EP -p -f2003 #-Wf"-ptr byte-pvctl noassume -pvctl chgpwr -pvctl loopcnt=200000" #-Wl"-Z 60G"
   F90_CFLAGS_OPT=-Chopt -EP -ftrace #-pvctl fullmsg
   F90_LFLAGS_OPT=-Chopt -ftrace
   ifeq ($(findstring F2003,$(HPC_TYPE)),F2003) #vector machine
      F90_CFLAGS_OPT+= -f2003
      F90_LFLAGS_OPT+= -f2003
      F90_CFLAGS_DEBUG+= -f2003
      F90_LFLAGS_DEBUG+= -f2003
   endif
#   ifeq ($(HPC_TYPE),sx9)
#   #-ftrace -mpitrace
#   else ($(HPC_TYPE),sx8)
#   #
#   endif

else 

#--------------------------------
# no vector machine

#--------------------------------
# INTEL 
   ifeq ($(findstring intel,$(COMPILER)),intel)
      COMP=mpif90
      CCOMP=mpicc

      ifeq ($(findstring smpss,$(COMPILER)),smpss)
        # find parameters with mpif90/cc -showme	
	SFDEF=-Wf,
	SCDEF=-Wc,
	SLDEF=-Wl,
	COMP=$(shell mpif90 -showme | sed 's: \(-[^ILlWDOg]\): $(SFDEF)\1:g' | sed 's/ifort/smpss-cc/')
	CCOMP=$(shell mpicc -showme | sed 's: \(-[^ILlWDOg]\): $(SCDEF)\1:g' | sed 's/icc/smpss-cc/') 
      endif

      F90_CFLAGS_DEBUG=$(SFDEF)-fpp -O0 -g $(SFDEF)-check all $(SFDEF)-r8 $(SFDEF)-i4 $(SFDEF)-traceback $(SFDEF)-implicitnone $(SFDEF)-warn $(SFDEF)all $(SFDEF)-keep
      F90_LFLAGS_DEBUG=-g -O0 $(SLDEF)-r8 $(SLDEF)-i4  -L$(MPI_LIB_DIR) -lmpi_f90  
      #F90_CFLAGS_DEBUG+=-warn alignments,declarations,errors,general,ignore_loc,interfaces,truncated_source,uncalled,nounused,usage

      F90_CFLAGS_OPT=$(SFDEF)-fpp $(SFDEF)-r8 $(SFDEF)-i4 -O3 $(SFDEF)-i-dynamic 
      F90_LFLAGS_OPT=$(SLDEF)-r8 $(SLDEF)-i4 -O3 $(SLDEF)-i-dynamic   -L$(MPI_LIB_DIR) -lmpi_f90
      #--------------------------------
      # PROFILING 
      ifeq ($(findstring profile,$(COMPILER)),profile)
         F90_CFLAGS_DEBUG+=$(SFDEF)-p -g
         F90_LFLAGS_DEBUG+=$(SLDEF)-p
         F90_CFLAGS_OPT+=$(SFDEF)-p -g
         F90_LFLAGS_OPT+=$(SLDEF)-p
      endif
      ifeq ($(findstring vamp,$(COMPILER)),vamp)
      COMP=vtf90 -vt:f90 mpif90
      CCOMP=vtcc -vt:cc mpicc
      endif
   endif

#--------------------------------
# Vampir

   ifeq ($(findstring vampir,$(COMPILER)),vampir)
      COMP=vtf90 -vt:f90 mpif90
      CCOMP=vtcc -vt:cc mpicc

      F90_CFLAGS_DEBUG=-xf95-cpp-input -O0 -g -fbounds-check -Wall -pedantic -fcray-pointer 
      F90_LFLAGS_DEBUG=-O0 -g -fbounds-check -Wall -pedantic 

#      F90_CFLAGS_OPT=-xf95-cpp-input -std=f2003 -O3 -fcray-pointer 
#      F90_LFLAGS_OPT=-std=f2003 -O3

      F90_CFLAGS_OPT=-warn all -fpp -r8 -i4 -O3 -i-dynamic 
      F90_LFLAGS_OPT=-warn all -r8 -i4 -O3 -i-dynamic 
 
      ifeq ($(findstring profile,$(COMPILER)),profile)
         F90_CFLAGS_DEBUG+=-pg
         F90_LFLAGS_DEBUG+=-pg
         F90_CFLAGS_OPT+=-pg
         F90_LFLAGS_OPT+=-pg
      endif
   endif

#--------------------------------
# GNU

    ifeq ($(findstring gnu,$(COMPILER)),gnu)
      ifeq ($(findstring xe6,$(HPC_TYPE)),xe6)
          COMP=ftn 
          CCOMP=cc 
          ifeq ($(findstring smpss,$(COMPILER)),smpss)
              # find parameters with mpif90/cc -showme	
              SFDEF=-Wf,
	      SCDEF=-Wc,
	      SLDEF=-Wl,
              COMP=smpss-cc -I$(MPICH_DIR)/include
	      CCOMP=smpss-cc -I$(MPICH_DIR)/include
          endif
          ifeq ($(findstring ompss,$(COMPILER)),ompss)
              # find parameters with mpif90/cc -showme	
              SFDEF=-Wf,
	      SCDEF=-Wc,
	      SLDEF=-Wl,


	      COMP=mfc --cc=ftn --ompss -I$(MPICH_DIR)/include -L$(MPICH_DIR)/lib -lfmpich -lmpich
	      CCOMP=mcc --cc=CC --ompss -I$(MPICH_DIR)/include -L$(MPICH_DIR)/lib -lfmpich -lmpich
          endif
      else
          COMP=mpif90 -O2
          CCOMP=mpicc -O2 

          ifeq ($(findstring smpss,$(COMPILER)),smpss)
            # find parameters with mpif90/cc -showme	
	    SFDEF=-Wf,
	    SCDEF=-Wc,
	    SLDEF=-Wl,
	    COMP=$(shell mpif90 -showme | sed 's: \(-[^ILlWDOg]\): $(SFDEF)\1:g' | sed 's/gfortran/smpss-cc2/')
	    CCOMP=$(shell mpicc -showme | sed 's: \(-[^ILlWDOg]\): $(SCDEF)\1:g' | sed 's/gcc/smpss-cc2/') 
          endif
          ifeq ($(findstring ompss,$(COMPILER)),ompss)
            # find parameters with mpif90/cc -showme	
	    #SFDEF=-Wf,
	    #SCDEF=-Wc,
	    #SLDEF=-Wl,
	    COMP=$(shell mpif90 -showme | sed 's: \(-[^ILlWDOg]\): $(SFDEF)\1:g' | sed 's/gfortran/mfc --ompss  -g3/') -O2  --instrumentation # -L/opt/mpi/openmpi/1.7.3-gnu-4.8.1/lib -lmpi -g3 -O0  -fPIC 
	    CCOMP=$(shell mpicc -showme | sed 's: \(-[^ILlWDOg]\): $(SCDEF)\1:g' | sed 's/gcc/mcc --ompss -g3 /')  -O2 --instrumentation # -L/opt/mpi/openmpi/1.7.3-gnu-4.8.1/lib -lmpi -g3 -O0  -fPIC  
          endif
      endif

      F90_CFLAGS_DEBUG+=$(SFDEF)-xf95-cpp-input -O0 -K -g3 -fcheck=all -fcheck=bounds -debug -Wall -fbacktrace $(SFDEF)-fbounds-check $(SFDEF)-Wall $(SFDEF)-pedantic $(SFDEF)-fcray-pointer $(SFDEF)-fwhole-file
      F90_LFLAGS_DEBUG+=-O0 -g3 -fbacktrace -K  -fcheck=all -Wall -fcheck=bounds -debug

      #F90_CFLAGS_OPT=$(SFDEF)-xf95-cpp-input -O3 $(SFDEF)-fcray-pointer $(SFDEF)-fwhole-file 
#-std=f2003
      F90_CFLAGS_OPT+=-O3
      F90_LFLAGS_OPT+=
      #--------------------------------
      # OpenMP
      ifeq ($(findstring omp,$(COMPILER)),omp)
	ifeq ($(findstring ompss,$(COMPILER)),ompss)
	# ignore, this is ompss
	else  # this is really omp
	   F90_CFLAGS_DEBUG += -fopenmp
           F90_LFLAGS_DEBUG += -fopenmp
           F90_CFLAGS_OPT   += -fopenmp
           F90_LFLAGS_OPT   += -fopenmp
	endif
      endif 	
      # PROFILING 
      ifeq ($(findstring profile,$(COMPILER)),profile)
         F90_CFLAGS_DEBUG+=$(SFDEF)-pg
         F90_LFLAGS_DEBUG+=$(SLDEF)-pg
         F90_CFLAGS_OPT+=$(SFDEF)-pg
         F90_LFLAGS_OPT+=$(SLDEF)-pg
      endif
     #
     # smpss/paraver tracing
     ifeq ($(findstring trace,$(COMPILER)),trace)
	F90_CFLAGS_DEBUG += -g       # paraver caller stack unwinding
	F90_LFLAGS_DEBUG += -g
	F90_CFLAGS_OPT   += -g
	F90_LFLAGS_OPT   += -g
	ifeq ($(findstring smpss,$(COMPILER)),smpss)  # smpss events
	  F90_CFLAGS_DEBUG += -t -k
	  F90_LFLAGS_DEBUG += -t -k
	  F90_CFLAGS_OPT   += -t -k
	  F90_LFLAGS_OPT   += -t -k
        endif
	ifeq ($(findstring ompss,$(COMPILER)),ompss)  # ompss events
          COMP = mfc --cc=ftn --ompss --instrumentation  -I$(MPICH_DIR)/include -L$(EXTRAE_HOME)/lib -lnanosmpitracef
          COMP = mfc --cc=CC --ompss --instrumentation  -I$(MPICH_DIR)/include -L$(EXTRAE_HOME)/lib -lnanosmpitracef
        endif
     endif
   endif

#--------------------------------
# Scalasca: simply add scalasca -instrument before command

   ifeq ($(findstring scalasca,$(COMPILER)),scalasca)
      COMP:=scalasca -instrument $(COMP) 
      CCOMP:=scalasca -instrument $(CCOMP) 
   endif



#--------------------------------
# LIBRARIES
   ifeq ($(findstring vampir,$(COMPILER)),vampir) 
      PROFEXT+=_vampir
   endif

   endif

   LIBEXT=_$(HPC_TYPE)$(PROFEXT)
   ifeq ($(HPC_TYPE),nehalem)
      LIBEXT=_$(HPC_TYPE)_x86_64$(PROFEXT)
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
CPP=cpp -C -traditional-cpp
%.f90: %.F90
	$(CPP) $(CPPDEFINES)  $< $(INCLUDES) | grep -v "^#" | grep -v "^$$" > $@  # note $$ is an  escaped $

%.o: %.f90
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
	- rm -f *.o *.mod *.i lbc *.mf95ss smpss_*.{f,F}90 *-adapters.c *-registration.c 
	@echo '==============================================================='
	@echo '=  c l e a n e d                                              ='
	@echo '==============================================================='

cleanall:
	- rm -f *.o *.mod *.i *.i90 lbc *.mf95ss smpss_*.{F,f}90 *-registration.c *-adapters.c *.tasks
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
