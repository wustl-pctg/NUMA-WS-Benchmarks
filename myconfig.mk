COMPLIER ?= Tapir
COMPILER_HOME ?= /usr/local/tapir/build

ifneq ($(CILKPLUS),1)
#Location of NUMA-WS Runtime
CILKRTS_HOME ?= /project/adams/home/j.deters/sandbox/locality_runtime
else
#Location of the Cilk Plus Runtime
CILKRTS_HOME ?= /project/adams/home/j.deters/sandbox/vanilla_runtime
endif
LIKWID_INCLUDE ?= /usr/local/include
CILK ?= “”
