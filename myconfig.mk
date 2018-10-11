COMPLIER ?= Tapir
COMPILER_HOME ?= #SET ME
ifndef LOCALITY
LOCALITY=0
endif
ifneq ($(CILKPLUS),1)
#Location of NUMA-WS Runtime
CILKRTS_HOME ?= #SET ME
else
#Location of the Cilk Plus Runtime
CILKRTS_HOME ?= #SET ME Optionally
endif
LIKWID_INCLUDE ?= #SET ME Optionally
CILK ?= “”
