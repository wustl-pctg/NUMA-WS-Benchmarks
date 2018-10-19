/*
 * Heat diffusion (Jacobi-type iteration)
 *
 * Usage: see function usage();
 *
 * Volker Strumpen, Boston                                 August 1996
 *
 * Copyright (c) 1996 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef SERIAL
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#else
#define cilk_spawn 
#define cilk_sync
#define __cilkrts_reset_timing()
#define __cilkrts_accum_timing()
#define __cilkrts_init()
#define __cilkrts_set_pinning_info(n)
#define __cilkrts_disable_nonlocal_steal()
#define __cilkrts_unset_pinning_info()
#define __cilkrts_enable_nonlocal_steal()
#define __cilkrts_pin_top_level_frame_at_socket(n)
#define __cilkrts_get_nworkers() 1
#define __cilkrts_num_sockets() 1
#endif

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>
#include <errno.h>
#include <numa.h>

#include <sys/mman.h>
#include <string.h>
#include <unistd.h>

#include "getoptions.h"
#include "ktiming.h"
#include "numa_allocate.h"

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

int pinning[4] = {0, 1, 2, 3};
//By setting num_sockets to 4 we make the assumption that we would like 
//to do a four way split at the top in all NO_PIN instances.
int num_sockets = 4;
extern int errno;

/* Define ERROR_SUMMARY if you want to check your numerical results */
//#define ERROR_SUMMARY

#define f(x,y)     (sin(x)*sin(y))
#define randa(x,t) (0.0)
#define randb(x,t) (exp(-2*(t))*sin(x))
#define randc(y,t) (0.0)
#define randd(y,t) (exp(-2*(t))*sin(y))
#define solu(x,y,t) (exp(-2*(t))*sin(x)*sin(y))

int nx, ny, nt;
double xu, xo, yu, yo, tu, to;
double dx, dy, dt;

double dtdxsq, dtdysq;
double t;

int leafmaxcol;

/*****************   Initialization of grid partition  NEW VERSION ********************/

void initgrid_v(double *old, int lb, int ub) {

  int a, b, llb, lub;

  llb = (lb == 0) ? 1 : lb;
  lub = (ub == nx) ? nx - 1 : ub;

  for (a=llb, b=0; a < lub; a++){		/* boundary nodes */
    old[a * ny + b] = randa(xu + a * dx, 0);
  }

  for (a=llb, b=ny-1; a < lub; a++){
    old[a * ny + b] = randb(xu + a * dx, 0);
  }

  if (lb == 0) {
    for (a=0, b=0; b < ny; b++){
      old[a * ny + b] = randc(yu + b * dy, 0);
    }
  }
  if (ub == nx) {
    for (a=nx-1, b=0; b < ny; b++){
      old[a * ny + b] = randd(yu + b * dy, 0);
    }
  }
  for (a=llb; a < lub; a++) { /* inner nodes */
    for (b=1; b < ny-1; b++) {
      old[a * ny + b] = f(xu + a * dx, yu + b * dy);
    }
  }
}


/***************** Five-Point-Stencil Computation NEW VERSION ********************/

void compstripe_v(double *neww, double *old, int lb, int ub) {
  int a, b, llb, lub;

  llb = (lb == 0) ? 1 : lb;
  lub = (ub == nx) ? nx - 1 : ub;

  for (a=llb; a < lub; a++) {
    for (b=1; b < ny-1; b++) {
        neww[a * ny + b] =   dtdxsq * (old[(a+1) * ny + b] - 2 * old[a * ny +b] + old[(a-1) * ny + b])
        + dtdysq * (old[a * ny + (b+1)] - 2 * old[a * ny + b] + old[a * ny + (b-1)])
        + old[a * ny + b];
    }
  }

  for (a=llb, b=ny-1; a < lub; a++)
    neww[a * ny + b] = randb(xu + a * dx, t);

  for (a=llb, b=0; a < lub; a++)
    neww[a * ny + b] = randa(xu + a * dx, t);

  if (lb == 0) {
    for (a=0, b=0; b < ny; b++)
      neww[a * ny + b] = randc(yu + b * dy, t);
  }
  if (ub == nx) {
    for (a=nx-1, b=0; b < ny; b++)
      neww[a * ny + b] = randd(yu + b * dy, t);
  }
}

/***************** Decomposition of 2D grids in stripes ********************/

#define ALLC       0
#define INIT       1
#define COMP       2

int divide_v(int lb, int ub, double *neww,
           double *old, int mode, int timestep){
  int l = 0, r = 0;

  if (ub - lb > leafmaxcol) {

    l = cilk_spawn divide_v(lb, (ub + lb) / 2, neww, old, mode, timestep);
    r = divide_v((ub + lb) / 2, ub, neww, old, mode, timestep);
    cilk_sync;

    return(l + r);

  } else {
    switch (mode) {
      case COMP:
        if (timestep % 2) //this sqitches back and forth between the two arrays by exchanging neww and old
          compstripe_v(neww, old, lb, ub);
        else
          compstripe_v(old, neww, lb, ub);
        return 1;

      case INIT:
        initgrid_v(old, lb, ub);
        return 1;
    }

    return 0;
  }
}

int divide_top_level_n_way(int lb, int ub, double *neww,
           double *old, int mode, int timestep){

  int l = 0, r = 0, lm = 0, rm = 0;

    __cilkrts_disable_nonlocal_steal();

    int chunk = (ub - lb) / num_sockets;

    int i; 
    for(i = 1; i < num_sockets; i++) {
       __cilkrts_set_pinning_info(i);
       cilk_spawn divide_v(lb + chunk * (i-1) , lb + chunk * i, neww, old, mode, timestep);
    }
    divide_v(lb + chunk * (i-1) , ub, neww, old, mode, timestep);
    
    __cilkrts_unset_pinning_info();
    __cilkrts_enable_nonlocal_steal();
    cilk_sync;

    return(l + lm + rm + r);
}

int divide_top_level_4_way(int lb, int ub, double *neww,
           double *old, int mode, int timestep){

  int l = 0, r = 0, lm = 0, rm = 0;

    __cilkrts_disable_nonlocal_steal();

    //Split 0
    __cilkrts_set_pinning_info(pinning[1]);
    l = cilk_spawn divide_v(lb, (ub + lb) / 4, neww, old, mode, timestep);

    //Split 1
    __cilkrts_set_pinning_info(pinning[2]);
    lm = cilk_spawn divide_v((ub + lb) / 4, (ub + lb) / 2, neww, old, mode, timestep);

    //Split 2
    __cilkrts_set_pinning_info(pinning[3]);
    rm = cilk_spawn divide_v((ub + lb) / 2, 3*(ub + lb) / 4, neww, old, mode, timestep);

    //Split 3
    __cilkrts_unset_pinning_info();
    r = divide_v(3*(ub + lb) / 4, ub, neww, old, mode, timestep);

    __cilkrts_enable_nonlocal_steal();
    cilk_sync;

    return(l + lm + rm + r);
}

int heat(double *old_v, double *new_v) {
  int  c;

#ifdef ERROR_SUMMARY
  double tmp, *mat;
  double mae = 0.0;
  double mre = 0.0;
  double me = 0.0;
  int a, b;
#endif

  /* Jacobi Iteration (divide x-dimension of 2D grid into stripes) */

  divide_top_level_n_way(0, nx, new_v, old_v, INIT, 0);
  for (c = 1; c <= nt; c++) {
    t = tu + c * dt;
        divide_top_level_n_way(0, nx, new_v, old_v, COMP, c);
  }

#ifdef ERROR_SUMMARY
  /* Error summary computation: Not parallelized! */
  mat = (c % 2) ? old_v : new_v;
  printf("\n Error summary of last time frame comparing with exact solution:");
  for (a=0; a<nx; a++)
    for (b=0; b<ny; b++) {
      tmp = fabs(mat[a * nx + b] - solu(xu + a * dx, yu + b * dy, to));
      //printf("error: %10e\n", tmp);
      if (tmp > mae) 
        mae = tmp;
    }

  printf("\n   Local maximal absolute error  %10e ", mae);

  for (a=0; a<nx; a++)
    for (b=0; b<ny; b++) {
      tmp = fabs(mat[a * nx + b] - solu(xu + a * dx, yu + b * dy, to));
      if (mat[a * nx + b] != 0.0)
        tmp = tmp / mat[a * nx + b];
      if (tmp > mre)
        mre = tmp;
    }

  printf("\n   Local maximal relative error  %10e %s ", mre * 100, "%");

  me = 0.0;
  for (a=0; a<nx; a++)
    for (b=0; b<ny; b++) {
      me += fabs(mat[a * nx + b] - solu(xu + a * dx, yu + b * dy, to));
    }

  me = me / (nx * ny);
  printf("\n   Global Mean absolute error    %10e\n\n", me);
#endif
  return 0;
}

int usage(void) {

  fprintf(stderr, "\nUsage: heat [<cilk-options>] [<options>}\n\n");
  fprintf(stderr, "This program uses a Jacobi-type iteration to "
      "solve a finite-difference\n");
  fprintf(stderr, "approximation of parabolic partial differential "
      "equations that models\n");
  fprintf(stderr, "for example the heat diffusion problem.\n\n");
  fprintf(stderr, "Optional parameter: \n");
  fprintf(stderr, "   -g #     "
      "granularity (columns per partition)  default: 10\n");
  fprintf(stderr, "   -nx #    "
      "total number of columns              default: 4096\n");
  fprintf(stderr, "   -ny #    "
      "total number of rows                 default: 512\n");
  fprintf(stderr, "   -nt #    "
      "total time steps                     default: 100\n");
  fprintf(stderr, "   -f filename    parameter file for nx, ny, ...\n");
  fprintf(stderr, "   -benchmark short/medium/long\n");
  return 1;
}

void read_heatparams(char *filefn) {

  FILE *f;
  int l;

  if ((f = fopen(filefn, "r")) == NULL) {
    printf("\n Can't open %s\n", filefn);
    exit(0);
  }
  l = fscanf(f, "%d %d %d %lf %lf %lf %lf %lf %lf",
      &nx, &ny, &nt, &xu, &xo, &yu, &yo, &tu, &to);
  if (l != 9)
    printf("\n Warning: fscanf errno %d", errno);
  fclose(f);

}

const char *specifiers[] = { "-g", "-nx", "-ny", "-nt", "-xu", "-xo", "-yu", "-yo", "-tu", "-to", "-f", "-benchmark", "-h", 0};
int opt_types[] = {INTARG, INTARG, INTARG, INTARG, DOUBLEARG, DOUBLEARG, DOUBLEARG, DOUBLEARG, DOUBLEARG, DOUBLEARG, STRINGARG, BENCHMARK, BOOLARG, 0 };

int main(int argc, char *argv[]) {

  int ret, benchmark, help;
  char filename[100];

  nx = 512;
  ny = 512;
  nt = 100;
  xu = 0.0;
  xo = 1.570796326794896558;
  yu = 0.0;
  yo = 1.570796326794896558;
  tu = 0.0;
  to = 0.0000001;
  leafmaxcol = 10;
  filename[0]=0;

  // use the math related function before parallel region;
  // there is some benigh race in initalization code for the math functions.
  fprintf(stderr, "Testing exp: %f\n", randb(nx, nt));

  get_options(argc, argv, specifiers, opt_types, &leafmaxcol,
              &nx, &ny, &nt, &xu, &xo, &yu, &yo, &tu, &to,
              filename, &benchmark, &help);

  if (help) return usage();

  if (benchmark) {
    switch (benchmark) {
      case 1:      /* short benchmark options -- a little work*/
        nx = 512;
        ny = 512;
        nt = 1;
        xu = 0.0;
        xo = 1.570796326794896558;
        yu = 0.0;
        yo = 1.570796326794896558;
        tu = 0.0;
        to = 0.0000001;
        leafmaxcol = 10;
        filename[0]=0;
        break;
      case 2:      /* standard benchmark options*/
        nx = 4096;
        ny = 512;
        nt = 40;
        xu = 0.0;
        xo = 1.570796326794896558;
        yu = 0.0;
        yo = 1.570796326794896558;
        tu = 0.0;
        to = 0.0000001;
        leafmaxcol = 10;
        filename[0]=0;
        break;
      case 3:      /* long benchmark options -- a lot of work*/
        nx = 4096;
        ny = 1024;
        nt = 100;
        xu = 0.0;
        xo = 1.570796326794896558;
        yu = 0.0;
        yo = 1.570796326794896558;
        tu = 0.0;
        to = 0.0000001;
        leafmaxcol = 1;
        filename[0]=0;
        break;
    }
  }

  if (filename[0]) read_heatparams(filename);

  dx = (xo - xu) / (nx - 1);
  dy = (yo - yu) / (ny - 1);
  dt = (to - tu) / nt;	/* nt effective time steps! */

  dtdxsq = dt / (dx * dx);
  dtdysq = dt / (dy * dy);

  __cilkrts_init();
  __cilkrts_pin_top_level_frame_at_socket(0);

  /* Memory Allocation */

  double *old_v, *new_v;
#ifndef NO_PIN
  num_sockets = __cilkrts_num_sockets();

  if(num_sockets == 2) { 
      pinning[0] = pinning[1] = 0;
      pinning[2] = pinning[3] = 1;
  } else if (num_sockets == 3) {
      pinning[3] = -1;
  }
#endif
  old_v = (double *) malloc(nx * ny * sizeof(double));
  new_v = (double *) malloc(nx * ny * sizeof(double));

  __cilkrts_reset_timing();
#if TIMING_COUNT
  clockmark_t begin, end;
  uint64_t elapsed[TIMING_COUNT];

  for(int i=0; i < TIMING_COUNT; i++) {
    begin = ktiming_getmark();
    ret = heat(old_v, new_v);
    end = ktiming_getmark();
    elapsed[i] = ktiming_diff_usec(&begin, &end);
  }
  print_runtime(elapsed, TIMING_COUNT);
#else
  ret = heat(old_v, new_v);
#endif

  __cilkrts_accum_timing();
  printf("\nCilk Example: heat\n");
  printf("\n   dx = %f", dx);
  printf("\n   dy = %f", dy);
  printf("\n   dt = %f", dt);

  printf("\n\n Stability Value for explicit method must be > 0:  %f\n\n",
      0.5 - (dt / (dx * dx) + dt / (dy * dy)));
  printf("Options: granularity = %d\n", leafmaxcol);
  printf("         nx          = %d\n", nx);
  printf("         ny          = %d\n", ny);
  printf("         nt          = %d\n", nt);

  return 0;
}
