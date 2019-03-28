COMPLIER ?= Tapir
COMPILER_HOME ?= 

ifneq ($(CILKPLUS),1)
#Location of NUMA-WS Runtime
CILKRTS_HOME ?= 
else
#Location of the Cilk Plus Runtime
CILKRTS_HOME ?= 
endif
LIKWID_INCLUDE ?= 
CILK ?= “”
