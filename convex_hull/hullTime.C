// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <algorithm>
#include "gettime.h"
#include "utils.h"
#include "geometry.h"
#include "parallel.h"
#include "geometryIO.h"
#include "parseCommandLine.h"
#include "hull.h"
#include <unistd.h>

#ifndef TIMING_COUNT
#define TIMING_COUNT 0
#endif

#if TIMING_COUNT
#include "ktiming.h"
#endif

#ifdef NO_PIN
#define __cilkrts_set_pinning_info(n)
#define __cilkrts_disable_nonlocal_steal()
#define __cilkrts_unset_pinning_info()
#define __cilkrts_enable_nonlocal_steal()
#define __cilkrts_pin_top_level_frame_at_socket(n)
#endif

#include <cilk/cilk_api.h>
#include <cilk/cilk.h>
#include "numa_allocate.h"

using namespace std;
using namespace benchIO;
static point2d *offset_helper(point2d * buf, size_t off){
  char *tmp = (char *)buf;
  tmp = tmp + off;
  return (point2d *)tmp;
}

void timeHull(point2d* P, intT n, int rounds, char* outFile) {
  intT m;
  _seq<intT> I;

  for (intT i=0; i < rounds; i++) {
    if (i != 0) I.del();
    startTime();
    I = hull(P, n);
    nextTimeN();
  }
  cout << endl;
  if (outFile != NULL) writeIntArrayToFile(I.A, I.n, outFile);
}

void timeHullP(point2d* P, intT n, int rounds, char* outFile) {
  intT m;
  _seq<point2d> ptmp;
  int num_pages = (4*n * sizeof(point2d)) / getpagesize() + 1;
#ifdef NO_PIN
  point2d *p = (point2d *)malloc(4 * n * sizeof(point2d));
  point2d *pin = (point2d *)malloc(4 * n * sizeof(point2d));
#else
  int num_sockets = __cilkrts_num_sockets();
  int mem_pattern[] = {0, 0, 0, 0};
  int pin_pattern[] = {0, 0, 0};
  if(num_sockets == 2){
    int mem_patternT[] = {0, 0, 1, 1};
    int pin_patternT[] = {0, 1, 1};
    memcpy(mem_pattern, mem_patternT, 4*sizeof(int));
    memcpy(pin_pattern, pin_patternT, 3*sizeof(int));
  }else if(num_sockets == 3){
    int mem_patternT[] = {0, 1, 2, -1};
    int pin_patternT[] = {1, 2, -1};
    memcpy(mem_pattern, mem_patternT, 4*sizeof(int));
    memcpy(pin_pattern, pin_patternT, 3*sizeof(int));
  }else if(num_sockets == 4){
    int mem_patternT[] = {0, 1, 2, 3};
    int pin_patternT[] = {1, 2, 3};
    memcpy(mem_pattern, mem_patternT, 4*sizeof(int));
    memcpy(pin_pattern, pin_patternT, 3*sizeof(int));
  }

  point2d *p = (point2d *)pattern_bind_memory_numa(num_pages, 4, mem_pattern);
  point2d *pin = (point2d *)pattern_bind_memory_numa(num_pages, 4, mem_pattern);
  __cilkrts_pin_top_level_frame_at_socket(0);
#endif

#if TIMING_COUNT
  clockmark_t begin, end;
  uint64_t elapsed[TIMING_COUNT];
  int i;
  for(i=0; i < TIMING_COUNT; i++) {
  cilk_for(int i = 0; i < n/4; i++){
    pin[i] = P[i];
  }

  cilk_for(int i = 0; i < n/4; i++){
    offset_helper(pin, alloc_size)[i] = P[n/4 + i];
  }

  cilk_for(int i = 0; i < n/4; i++){
    offset_helper(pin, alloc_size*2)[i] = P[n/4*2+i];
  }

  cilk_for(int i = 0; i < n - n/4*3; i++){
    offset_helper(pin, alloc_size*3)[i] = P[n/4*3+i];
  }
  __cilkrts_reset_timing();
  begin = ktiming_getmark();
  ptmp = hullP(pin, n, p);
  end = ktiming_getmark();
  elapsed[i] = ktiming_diff_usec(&begin, &end);
  }
  print_runtime(elapsed, TIMING_COUNT);
#else
  intT alloc_size = (n * sizeof(point2d));
  cilk_for(int i = 0; i < n/4; i++){
    pin[i] = P[i];
  }

  cilk_for(int i = 0; i < n/4; i++){
    offset_helper(pin, alloc_size)[i] = P[n/4 + i];
  }

  cilk_for(int i = 0; i < n/4; i++){
    offset_helper(pin, alloc_size*2)[i] = P[n/4*2+i];
  }

  cilk_for(int i = 0; i < n - n/4*3; i++){
    offset_helper(pin, alloc_size*3)[i] = P[n/4*3+i];
  }
  __cilkrts_reset_timing();
  ptmp = hullP(pin, n, p);
#endif
  __cilkrts_accum_timing();
  if (outFile != NULL) writePointsToFile(ptmp.A, ptmp.n, outFile);
}


int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);

  _seq<point2d> PIn = readPointsFromFile<point2d>(iFile);

  timeHullP(PIn.A, PIn.n, rounds, oFile);
}
