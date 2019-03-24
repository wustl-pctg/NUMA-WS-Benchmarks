-*- mode: makefile-gmake; -*-

#Configure your myconfig.mk
TEST_MK_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
include $(TEST_MK_DIR)/myconfig.mk

OPT_CFLAGS ?= -O3
CILK_CFLAGS = -fcilkplus $(ADDED_FLAGS)
INCLUDE=$(strip $(CILKRTS_HOME))/include -I/usr/local/include
#Werror
ifeq ($(LIKWID), 1)
 BASIC_CFLAGS = $(OPT_CFLAGS) -W -Wall -DLIKWID_PERFMON -I$(INCLUDE)
else
 BASIC_CFLAGS = $(OPT_CFLAGS) -W -Wall -I$(INCLUDE)
endif
BASIC_CXXFLAGS=$(BASIC_CFLAGS)
CFLAGS = $(BASIC_CFLAGS) $(CILK_CFLAGS) -std=c11
CXXFLAGS = $(BASIC_CXXFLAGS) $(CILK_CFLAGS) -std=c++11

ifdef CILKRR_HOME
 export CPATH+=$(CILKRR_HOME)
 CILKRR_SLIB=$(CILKRR_HOME)/libcilkrr.a
 LDFLAGS = $(CILKRR_SLIB)
endif

# CILKRTS related stuff --- for dynalic linking use the following
CILKRTS_DLIB=-Wl,-rpath -Wl,$(strip $(CILKRTS_HOME))/lib
LDFLAGS += -L$(strip $(CILKRTS_HOME))/lib $(CILKRTS_DLIB) -L/usr/local/lib -Wl, -rpath -Wl,/usr/local/lib

#use LIKWID=1 to enable likwid
ifeq ($(LIKWID),1)
 LDLIBS += -lm -lrt -ldl -lpthread -lcilkrts -llikwid -lnuma
else
 LDLIBS += -lm -lrt -ldl -lpthread -lcilkrts -lnuma
endif

#configure the number of CPUs per socket for benchmarks that interleave when
#using the Cilk Plus runtime.
CPUS_PER_SOCKET = 8
CFLAGS += -DCPUS_PER_SOCKET=$(CPUS_PER_SOCKET)
CXXFLAGS += -DCPUS_PER_SOCKET=$(CPUS_PER_SOCKET)

#use NO_PIN=1 to remove NUMA-WS runtime calls
ifndef NO_PIN
NO_PIN = 0
endif

ifeq ($(NO_PIN),1)
        CFLAGS += -DNO_PIN
        CXXFLAGS += -DNO_PIN
endif


ifndef DISABLE_NONLOCAL_STEAL
DISABLE_NONLOCAL_STEAL = 1
endif

ifeq ($(DISABLE_NONLOCAL_STEAL),1)
        CFLAGS += -DDISABLE_NONLOCAL_STEAL
        CXXFLAGS += -DDISABLE_NONLOCAL_STEAL
endif


#use TIMING_COUNT=num_runs to enable timings
ifndef TIMING_COUNT
TIMING_COUNT=0
endif

ifneq ($(TIMING_COUNT),0)
       TIMING_CODE=ktiming.o
       CFLAGS += -DTIMING_COUNT=$(TIMING_COUNT)
       CXXFLAGS += -DTIMING_COUNT=$(TIMING_COUNT)
endif

#use SERIAL=1 to enable serial execution
ifndef SERIAL
SERIAL=0
endif

ifneq ($(SERIAL),0)
   CFLAGS += -DSERIAL
   CXXFLAGS += -DSERIAL
endif

#use POS_2=1 to enable the alternative position for enable nonlocal steal
ifndef POS_2
POS_2=0
endif

ifneq ($(POS_2),0)
   CFLAGS += -DPOS_2
   CXXFLAGS += -DPOS_2
endif

CC = $(COMPILER_HOME)/bin/clang
CXX = $(COMPILER_HOME)/bin/clang++

.PHONY : default clean

default : $(TARGETS)

# Each C source file will have a corresponding file of prerequisites.
# Include the prerequisites for each of our C source files.
-include $(OBJ:.o=.d)

# This rule generates a file of prerequisites (i.e., a makefile)
# called name.d from a C source file name.c.
%.d : CFLAGS += -MM -MP
%.d : %.c
	@set -e; rm -f $@; \
	$(CC) $(CFLAGS) -MF $@.$$$$ $<; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# This rule generates a file of prerequisites (i.e., a makefile)
# called name.d from a CPP source file name.cpp.
%.d : CXXFLAGS += -MM -MP
%.d : %.cpp
	@set -e; rm -f $@; \
	$(CXX) $(CXXFLAGS) -MF $@.$$$$ $<; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%.o : %.c
	$(CC) $(CFLAGS) -c $<

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $<

veryclean::
	cd ../;rm -r bin/
