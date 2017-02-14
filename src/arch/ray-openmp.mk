#-------------------------------------------------------------------------------

CXX=mpiclang++-gpu 
CC=mpiclang-gpu --std=gnu99 
LD=$(CXX)

DFLAGS = -DPELOTON -DWITH_PIO -DWITH_MPI \
	-DADD_ -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

INCLUDE =

CFLAGS_BASE = $(INCLUDE) $(DFLAGS)
CXXFLAGS_BASE = $(INCLUDE) $(DFLAGS)

CFLAGS_OPT =   $(CFLAGS_BASE) -g -O3 -ffp-contract=fast -fopenmp-nonaliased-maps
CFLAGS_DEBUG = $(CFLAGS_BASE) -g -ggdb -O0 -fno-inline
CFLAGS_PROF =  $(CFLAGS_BASE) -g -pg -O3 -DPROFILE

CXXFLAGS_OPT =   $(CXXFLAGS_BASE) -g -O3 -ffp-contract=fast -fopenmp-nonaliased-maps
CXXFLAGS_DEBUG = $(CXXFLAGS_BASE) -g -ggdb -O0 -fno-inline
CXXFLAGS_PROF =  $(CXXFLAGS_BASE) -g -pg -O3 -DPROFILE

LDFLAGS_OPT   = $(LDFLAGS_BASE) $(CFLAGS_OPT) $(CXXFLAGS_OPT)
LDFLAGS_DEBUG = $(LDFLAGS_BASE) $(CFLAGS_DEBUG) $(CXXFLAGS_DEBUG)
LDFLAGS_PROF  = $(LDFLAGS_BASE) $(CFLAGS_PROF) $(CXXFLAGS_PROF)

OMPFLAGS= $(CXX_FLAGS) #-Xcuda-ptxas -maxrregcount=96

#-------------------------------------------------------------------------------

