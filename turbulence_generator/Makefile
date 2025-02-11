#
#  Makefile for turbulence generator
#
#  Written by
#
#          Christoph Federrath
#

# general complile switches
HAVE_MPI = no
HAVE_HDF5 = yes

# C++ compiler (e.g., g++), flags, and HDF5 library path
CCOMP = mpicxx
CFLAGS = -O3
#HDF5_PATH = /opt/local
#HDF5_PATH = /usr/local/lib/HDF5-parallel
#HDF5_PATH = /usr/local/lib/HDF5-serial/lib
#HDF5_PATH = /usr/lib/x86_64-linux-gnu/hdf5/openmpi
HDF5_PATH = /usr/lib/x86_64-linux-gnu/hdf5/serial


ifeq ($(strip $(HAVE_MPI)), yes)
CFLAGS += -DHAVE_MPI
endif
ifeq ($(strip $(HAVE_HDF5)), yes)
CFLAGS += -DHAVE_HDF5 -I$(HDF5_PATH)/include -L$(HDF5_PATH)/lib -lhdf5
#CFLAGS += -DHAVE_HDF5 -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -lcrypto -lcurl -lpthread -lsz -lz -ldl -lm
endif

BINS = TurbGen TurbGenDemo

all: $(BINS)

clean:
	rm -f *.o *~ $(BINS)

$(BINS): %: %.cpp
	$(CCOMP) $(CFLAGS) -o $@ $<

# dependencies
$(BINS).o: TurbGen.h HDFIO.h
