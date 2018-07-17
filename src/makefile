SMILEICXX		= mpicxx
HDF5_DIR		?=
SuperLU_DIR		?=
PYTHONCONFIG	= python3-config
##################################################
CFLAGS 			+= -std=c++11


DESCRIBE:=$(shell git describe 2>/dev/null || echo '??')
BRANCH:=$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '??')
COMMITDATE:=$(shell git show -s --pretty="%ci" 2>/dev/null || echo '??')
VERSION:="$(DESCRIBE)-$(BRANCH)"

CFLAGS += -D__VERSION=\"$(VERSION)\" -D__COMMITDATE=\""$(COMMITDATE)"\"

####################################################
ifeq ($(SuperLU_serial),yes)
CFLAGS 			   += -DSuperLU_serial
LIB_Blas			= $(SuperLU_dir)/lib/libblas.a
LIB_SuperLU 		= $(SuperLU_dir)/lib/libsuperlu_4.3.a
LIBS 				= $(LIB_SuperLU) $(LIB_Blas)
CFLAGS 			   += -I$(SuperLU_dir)/SRC
$(warning  $(CFLAGS))
endif

ifeq ($(SuperLU_mpi),yes)
CFLAGS 			   += -DSuperLU_mpi
LIB_Blas			= $(SuperLU_DIST_dir)/lib/libblas.a
LIB_SuperLU_DIST 	= $(SuperLU_DIST_dir)/lib/libsuperlu_dist.a
LIBS 				= $(LIB_SuperLU_DIST) $(LIB_Blas)
CFLAGS 		       += -I$(SuperLU_DIST_dir)/SRC
$(warning  $(CFLAGS))
endif

# for PETSc library
ifeq ($(petsc),yes)
include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS 		+= -Dpetsc
CFLAGS 		+= -I${petsc_dir}/include
LDFLAGS 	+= -L${petsc_dir}/lib -lpetsc -lz
endif


# ============== hdf5 libary =========== 
CFLAGS 		+= -I${hdf5_dir}/include
LDFLAGS 	+= -L${hdf5_dir}/lib -lhdf5 -lz

LDFLAGS 	+= $(LIBS)

####################################################
#add subdirs
DIRS 		= $(shell find * -type d)
#add include directives for subdirs
CFLAGS 		+= $(DIRS:%=-I%)



PYSCRIPTS := $(shell find Python -name \*.py)
PYHEADERS := $(PYSCRIPTS:.py=.pyh)

#collect all cpp files
SRCS := $(shell find . -name \*.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(SRCS:.cpp=.d)
PYCS := $(shell find . -name \*.pyc)
PYCS += $(shell find ../tools -name \*.pyc)

CFLAGS		+=$(shell $(PYTHONCONFIG) --includes)
PATH_Python = $(shell $(PYTHONCONFIG) --prefix)
LDFLAGS 	+= -L${PATH_Python}/lib $(shell $(PYTHONCONFIG) --ldflags)


ifeq ($(config),debug)
	CFLAGS += -g -pg -Wall -D__DEBUG -O3# -shared-intel
else
	CFLAGS += -O3#  -xHost -ipo
endif



###############################################################
EXEC		= liz
default: $(EXEC)

help:
	@echo 'Usage                 : make [clean] [openmp={gnu|intel}] [config=debug]'
	@echo 'Environment variables : SMILEICXX (c++ compiler) HDF5_ROOT_DIR (HDF5 dir with /include and /lib)'

release:
	make

debug:
	make config=debug

clean:
	rm -f $(OBJS) $(DEPS) $(PYHEADERS) $(EXEC) $(PYCS)

distclean: clean
	rm -f $(EXEC)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

# this generates a .h file containing a char[] with the python script in binary then
#you can just include this file to get the contents
%.pyh: %.py
	@ echo "Creating binary char for $< : $@"
	@ xxd -i $< > $@

%.d: %.cpp $(PYHEADERS)
	@ echo "Checking dependencies for $<"
# create and modify dependecy file .d to take into account the location subdir
	@ $(SMILEICXX) $(CFLAGS) -MM $< 2>/dev/null | sed -e "s@\(^.*\)\.o:@$(shell  dirname $<)/\1.d $(shell  dirname $<)/\1.o:@" > $@

%.o : %.cpp
	$(SMILEICXX) $(CFLAGS) -c $< -o $@

$(EXEC): $(PYHEADERS) $(OBJS)
	$(SMILEICXX) $(CFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

# DO NOT DELETE
