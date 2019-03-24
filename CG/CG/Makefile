SHELL=/bin/sh
BENCHMARK=cg
BENCHMARKU=CG

include ../config/make.def


#use TIMING_COUNT=num_runs to enable timings
ifndef TIMING_COUNT
TIMING_COUNT=0
endif

ifneq ($(TIMING_COUNT),0)
       CFLAGS += -DTIMING_COUNT=$(TIMING_COUNT)
endif

#use NO_PIN=1 to remove NUMA-WS runtime calls
ifndef NO_PIN
NO_PIN = 0
endif

ifeq ($(NO_PIN),1)
        CFLAGS += -DNO_PIN
endif

#use SERIAL=1 to enable serial execution
ifndef SERIAL
SERIAL=0
endif

ifneq ($(SERIAL),0)
   CFLAGS += -DSERIAL
endif


#use LIKWID=1 to enable likwid
ifeq ($(LIKWID),1)
C_LIB += -llikwid
endif

#configure the number of CPUs per socket for benchmarks that interleave when
#using the Cilk Plus runtime.
CPUS_PER_SOCKET = 8
CFLAGS += -DCPUS_PER_SOCKET=$(CPUS_PER_SOCKET)
CXXFLAGS += -DCPUS_PER_SOCKET=$(CPUS_PER_SOCKET)


OBJS = cg.o ${COMMON}/c_print_results.o  \
       ${COMMON}/c_${RAND}.o ${COMMON}/c_timers.o ${COMMON}/c_wtime.o numa_allocate.o ktiming.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

cg.o:		cg.c  npbparams.h
	${CCOMPILE} $(CFLAGS) cg.c

numa_allocate.o:
	${CCOMPILE} numa_allocate.c

ktiming.o:
	${CCOMPILE} ktiming.c

clean:
	- rm -f *.o *~
	- rm -f npbparams.h core
