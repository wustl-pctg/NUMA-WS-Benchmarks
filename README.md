# NUMA-WS Benchmarks

These are a series of benchmarks to test and demonstrate the capability of the
NUMA-WS runtime system. Read more [here](https://ieeexplore.ieee.org/document/8573486).

## Using the runtime
Currently these are *hand compiled* benchmarks with the runtime calls that manage
the NUMA aware portions of the runtime inserted with great care.

## Getting Started
1. **Setting up the Tapir compiler:** It is highly recommended that these benchmarks
are compiled using the MIT Tapir compiler which can be found
[here](http://cilk.mit.edu/download/). The rest of this document (and the configuration
files in the repo) will assume that you're using this compiler. In **myconfig.mk**
set COMPILER_HOME to the location of your Tapir compiler
(ex: COMPILER_HOME ?= /usr/local/tapir/build).

2. **Setting up the NUMA-WA Runtime:** Please follow the instructions on how to
configure and install the NUMA-WS runtime that are included in its README. Once
you have done that, set the CILKRTS_HOME variable in **myconfig.mk** designated
as *Location of NUMA-WS Runtime* to the location of your Cilk Plus Runtime.

3. **Setting up the Cilk Plus Runtime:** This is an optional step. If you choose
to use these benchmarks with the vanilla Cilk Plus runtime set the CILKRTS_HOME
variable in **myconfig.mk** designated as *Location of the Cilk Plus Runtime* to
the location of your Cilk Plus Runtime.

    If you choose to use this option. When compiling use CILKPLUS=1 compiler command
to link with your Cilk Plus runtime and the NO_PIN=1 command to remove the NUMA-WS
runtime calls from the benchmarks. (ex: make CILKPLUS=1 NO_PIN=1) On our system, we found that cilksort and cg performed best when using interleaved
memory in Cilk Plus. Please set the CPUS_PER_SOCKET variable in common.mk to
make sure that the interleaving works correctly.

4. **Setting up LIKWID:** This is an optional step. If you wish to use the
[LIKWID](https://github.com/RRZE-HPC/likwid) profiling system with these benchmarks
set LIKWID_INCLUDE in **myconfig.mk** to the location of your LIKWID installation
(ex: /usr/local/include). When compiling be sure to set LIKWID=1 to enable the
functionality (ex: make LIKWID=1)

## Compilation Commands
- CILKPLUS
    - If 0 compile with link to the NUMA-WS runtime.
    - If 1 compile with link to the Cilk Plus runtime.
- NO_PIN
    - If 0 compile with NUMA-WS runtime calls.
    - If 1 compile without NUMA-WS runtime calls.
- LIKWID
    - If 0 do nothing.
    - If 1 compile with link to the LIKWID profiler.
- SERIAL
    - If 0 do nothing.
    - If 1 compile with all cilk removed and run serially.
- TIMING_COUNT
    - If 0 do nothing.
    - If > 0 compile with timing enabled and run each benchmark n times, where
    n is the value of TIMING_COUNT.
- make move

    This moves all the compiled binaries into the /bin directory at the top level of the repo. This is used for convenient testing.

- make veryclean

    This will do a normal clean of the repo, but also remove the binaries in /bin.

## Notes on Convex Hull and CG
These two benchmarks exist within their own directories. To use them, you will
need to compile them separately from the benchmarks that are in the /benchmarks
directory. However, the still respond to all the same compilation commands as shown
above.

**Data generation** for Convex Hull can be found in the top level /data_generation
directory.

## Cititation
Please use the following citiation when using this runtime in your work:

J. Deters, J. Wu, Y. Xu and I. Lee, "A NUMA-Aware Provably-Efficient Task-Parallel Platform Based on the Work-First Principle", 
2018 IEEE International Symposium on Workload Characterization (IISWC), 2018. 
Available: https://ieeexplore.ieee.org/document/8573486.

## Acknowledgment
This research was supported in part by National Science Foundation under grant number CCF-1527692 and CCF-1733873.

Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) 
and do not necessarily reflect the views of the National Science Foundation.
