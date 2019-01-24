/*--------------------------------------------------------------------

  NAS Parallel Benchmarks 3.0 structured OpenMP C versions - CG

  This benchmark is an OpenMP C version of the NPB CG code.

  The OpenMP C 2.3 versions are derived by RWCP from the serial Fortran versions
  in "NPB 2.3-serial" developed by NAS. 3.0 translation is performed by the UVSQ.

  Permission to use, copy, distribute and modify this software for any
  purpose with or without fee is hereby granted.
  This software is provided "as is" without express or implied warranty.

  Information on OpenMP activities at RWCP is available at:

           http://pdplab.trc.rwcp.or.jp/pdperf/Omni/

  Information on NAS Parallel Benchmarks 2.3 is available at:

           http://www.nas.nasa.gov/NAS/NPB/

--------------------------------------------------------------------*/
/*--------------------------------------------------------------------

  Authors: M. Yarrow
           C. Kuszmaul

  OpenMP C version: S. Satoh

  3.0 structure translation: F. Conti

--------------------------------------------------------------------*/

/*
c---------------------------------------------------------------------
c  Note: please observe that in the routine conj_grad three
c  implementations of the sparse matrix-vector multiply have
c  been supplied.  The default matrix-vector multiply is not
c  loop unrolled.  The alternate implementations are unrolled
c  to a depth of 2 and unrolled to a depth of 8.  Please
c  experiment with these to find the fastest for your particular
c  architecture.  If reporting timing results, any of these three may
c  be used without penalty.
c---------------------------------------------------------------------
*/
#ifndef SERIAL
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#else
#define cilk_spawn
#define cilk_sync
#define __cilkrts_accum_timing()
#define __cilkrts_set_pinning_info(n)
#define __cilkrts_disable_nonlocal_steal()
#define __cilkrts_unset_pinning_info()
#define __cilkrts_enable_nonlocal_steal()
#define __cilkrts_pin_top_level_frame_at_socket(n)
#define __cilkrts_init()
#define __cilkrts_get_nworkers() 1
#define __cilkrts_num_sockets() 1
#define __cilkrts_reset_timing()
#define cilk_for for
#endif

#ifndef DISABLE_NONLOCAL_STEAL
#define __cilkrts_disable_nonlocal_steal()
#define __cilkrts_enable_nonlocal_steal()
#endif

#include "npb-C.h"
#include "npbparams.h"
#include "numa_allocate.h"
#include <unistd.h>
#include <sched.h>
#include <numa.h>
#include <numaif.h>

#ifndef TIMING_COUNT
#define TIMING_COUNT 0
#endif

#if TIMING_COUNT
#include "ktiming.h"

clockmark_t begin, end;
uint64_t *elapsed;
#endif

#define	NZ	NA*(NONZER+1)*(NONZER+1)+NA*(NONZER+2)

/* global variables */

/* common /partit_size/ */
static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;

/* common /main_int_mem/ */
static int colidx[NZ+1];	/* colidx[1:NZ] */
static int rowstr[NA+1+1];	/* rowstr[1:NA+1] */
static int iv[2*NA+1+1];	/* iv[1:2*NA+1] */
static int arow[NZ+1];		/* arow[1:NZ] */
static int acol[NZ+1];		/* acol[1:NZ] */

/* common /main_flt_mem/ */
static double v[NA+1+1];	/* v[1:NA+1] */
static double aelt[NZ+1];	/* aelt[1:NZ] */
static double a[NZ+1];		/* a[1:NZ] */
static double x[NA+2+1];	/* x[1:NA+2] */
static double z[NA+2+1];	/* z[1:NA+2] */
static double p[NA+2+1];	/* p[1:NA+2] */
static double q[NA+2+1];	/* q[1:NA+2] */
static double r[NA+2+1];	/* r[1:NA+2] */
//static double w[NA+2+1];	/* w[1:NA+2] */

/* common /urando/ */
static double amult;
static double tran;

/* function declarations */
static void conj_grad (int colidx[], int rowstr[], double x[], double z[],
		       double a[], double p[], double q[], double r[],
		       //double w[],
		       double *rnorm);
static void makea(int n, int nz, double a[], int colidx[], int rowstr[],
		  int nonzer, int firstrow, int lastrow, int firstcol,
		  int lastcol, double rcond, int arow[], int acol[],
		  double aelt[], double v[], int iv[], double shift );
static void sparse(double a[], int colidx[], int rowstr[], int n,
		   int arow[], int acol[], double aelt[],
		   int firstrow, int lastrow,
		   double x[], boolean mark[], int nzloc[], int nnza);
static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[],
		   int mark[]);
static int icnvrt(double x, int ipwr2);
const unsigned long BASE_CASE = 512;

double reduce_add_mul(double init_sum, double *arr1, double *arr2, unsigned long i, unsigned long e){
  if( e - i < BASE_CASE){
    for(int j = i; j < e; j++){
      init_sum += arr1[j] * arr2[j];
    }
    return init_sum;
  }

  unsigned long mid = (i + e) / 2;
  double k = cilk_spawn reduce_add_mul(init_sum, arr1, arr2, i, mid);
  double l = reduce_add_mul(init_sum, arr1, arr2, mid, e);
  cilk_sync;
  return k+l;
}

#ifndef NO_PIN
void map_mul_scalar(double *res, double *A, double alpha, unsigned long i, unsigned long e){
  unsigned long mid = (i + e) / 2;
  if( e - i < BASE_CASE){
    for(int j = i; j < e; j++){
      res[j] = A[j] * alpha;
    }
    return;
  }
  cilk_spawn map_mul_scalar(res, A, alpha, i, mid);
  map_mul_scalar(res, A, alpha, mid, e);
  cilk_sync;
  return;
}
#endif

const unsigned long IDX_MUL_BASE_CASE = BASE_CASE;
double reduce_add_mul_idx(double init_sum, double *arr1, double *arr2, int *idx, unsigned long i, unsigned long e){
  unsigned long mid = (i + e) / 2;
  //printf("here: e- i : %d\n", e-i);
  if( e - i < IDX_MUL_BASE_CASE){
    for(int j = i; j < e; j++){
      init_sum += arr1[j] * arr2[idx[j]];
    }
    return init_sum;
  }
  double k = cilk_spawn reduce_add_mul_idx(init_sum, arr1, arr2, idx, i, mid);
  double l = reduce_add_mul_idx(init_sum, arr1, arr2, idx, mid, e);
  cilk_sync;
  return k+l;
}

static void vecset(int n, double v[], int iv[], int *nzv, int i, double val);

#ifndef NO_PIN
static int sockets;
static int mem_pattern[] = {0, 0, 0, 0, 0, 0, 0, 0};
static int pin_pattern[] = {0, 0, 0, 0, 0, 0, 0, 0};
#endif

/*--------------------------------------------------------------------
      program cg
--------------------------------------------------------------------*/
#ifndef NO_PIN
#define SET_PIN(N) __cilkrts_set_pinning_info(N)
#endif

int fakemain(int argc, char **argv, int run) {
  int	i, j, k, it;
  int nthreads = 1;
    double zeta;
    double rnorm;
    double norm_temp11;
    double norm_temp12;
    double t, mflops;
    char class;
    boolean verified;
    double zeta_verify_value, epsilon;

    firstrow = 1;
    lastrow  = NA;
    firstcol = 1;
    lastcol  = NA;

    if (NA == 1400 && NONZER == 7 && NITER == 15 && SHIFT == 10.0) {
	class = 'S';
	zeta_verify_value = 8.5971775078648;
    } else if (NA == 7000 && NONZER == 8 && NITER == 15 && SHIFT == 12.0) {
	class = 'W';
	zeta_verify_value = 10.362595087124;
    } else if (NA == 14000 && NONZER == 11 && NITER == 15 && SHIFT == 20.0) {
	class = 'A';
	zeta_verify_value = 17.130235054029;
    } else if (NA == 75000 && NONZER == 13 && NITER == 75 && SHIFT == 60.0) {
	class = 'B';
	zeta_verify_value = 22.712745482631;
    } else if (NA == 150000 && NONZER == 15 && NITER == 75 && SHIFT == 110.0) {
	class = 'C';
	zeta_verify_value = 28.973605592845;
    } else {
	class = 'U';
    }

    //printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"
	  // " - CG Benchmark\n");
    //printf(" Size: %10d\n", NA);
    //printf(" Iterations: %5d\n", NITER);

    naa = NA;
    nzz = NZ;

/*--------------------------------------------------------------------
c  Initialize random number generator
c-------------------------------------------------------------------*/
    tran    = 314159265.0;
    amult   = 1220703125.0;
    zeta    = randlc( &tran, amult );

/*--------------------------------------------------------------------
c
c-------------------------------------------------------------------*/
    //printf("data gen\n");
    makea(naa, nzz, a, colidx, rowstr, NONZER,
	  firstrow, lastrow, firstcol, lastcol,
	  RCOND, arow, acol, aelt, v, iv, SHIFT);
    //printf("data gen end\n");

/*---------------------------------------------------------------------
c  Note: as a result of the above call to makea:
c        values of j used in indexing rowstr go from 1 --> lastrow-firstrow+1
c        values of colidx which are col indexes go from firstcol --> lastcol
c        So:
c        Shift the col index vals from actual (firstcol --> lastcol )
c        to local, i.e., (1 --> lastcol-firstcol+1)
c---------------------------------------------------------------------*/
    cilk_for (unsigned long j = 1; j <= lastrow - firstrow + 1; j++) {
	for (unsigned long k = rowstr[j]; k < rowstr[j+1]; k++) {
            colidx[k] = colidx[k] - firstcol + 1;
	}
    }
    //printf("end colidx\n");
/*--------------------------------------------------------------------
c  set starting vector to (1, 1, .... 1)
c-------------------------------------------------------------------*/
    cilk_for (unsigned long i = 1; i <= NA+1; i++) {
	x[i] = 1.0;
    }
    cilk_for (unsigned long j = 1; j <= lastcol-firstcol+1; j++) {
         q[j] = 0.0;
         z[j] = 0.0;
         r[j] = 0.0;
         p[j] = 0.0;
      }
    zeta  = 0.0;

/*-------------------------------------------------------------------
c---->
c  Do one iteration untimed to init all code and data page tables
c---->                    (then reinit, start timing, to niter its)
c-------------------------------------------------------------------*/
#ifndef NO_PIN
    int rNZ = rowstr[lastrow-firstrow+2];
    int sockets = __cilkrts_num_sockets();
    //colidx
    unsigned long long colidx_pages = (rNZ+1)*sizeof(int)/getpagesize()+1;
    int *colidx_numa = (int *) pattern_bind_memory_numa(colidx_pages, sockets, mem_pattern);
    cilk_for(int i =0;i < rNZ+1; i++){
      colidx_numa[i] = colidx[i];
    }
    //rowstr
    unsigned long long rowstr_pages = (NA+1+1)*sizeof(int)/getpagesize()+1;
    //printf("rowstr pages: %d\n", rowstr_pages);
    int *rowstr_numa = (int *) pattern_bind_memory_numa(rowstr_pages, sockets, mem_pattern);
    cilk_for(int i =0;i < NA+1+1; i++){
      rowstr_numa[i] = rowstr[i];
    }

    //x
    unsigned long long x_pages = (NA+2+1)*sizeof(double)/getpagesize()+1;
    double *x_numa = (double *) pattern_bind_memory_numa(x_pages, sockets, mem_pattern);
    cilk_for(int i =0;i < NA+2+1; i++){
      x_numa[i] = x[i];
    }

    //z
    unsigned long long z_pages = (NA+2+1)*sizeof(double)/getpagesize()+1;
    double *z_numa = (double *) pattern_bind_memory_numa(z_pages, sockets, mem_pattern);
    cilk_for(int i =0;i < NA+2+1; i++){
      z_numa[i] = z[i];
    }

    //a
    unsigned long long a_pages = (rNZ+1)*sizeof(double)/getpagesize()+1;
    double *a_numa = (double *) pattern_bind_memory_numa(a_pages, sockets, mem_pattern);
    cilk_for(int i =0;i < rNZ+1; i++){
      a_numa[i] = a[i];
    }
    //p
    unsigned long long p_pages = (NA+2+1)*sizeof(double)/getpagesize()+1;
    double *p_numa = (double *) pattern_bind_memory_numa(p_pages, sockets, mem_pattern);
    cilk_for(int i =0;i < NA+2+1; i++){
      p_numa[i] = p[i];
    }

    //q
    unsigned long long q_pages = (NA+2+1)*sizeof(double)/getpagesize()+1;
    double *q_numa = (double *) pattern_bind_memory_numa(q_pages, sockets, mem_pattern);
    cilk_for(int i =0;i < NA+2+1; i++){
      q_numa[i] = q[i];
    }

    //r
    unsigned long long r_pages = (NA+2+1)*sizeof(double)/getpagesize()+1;
    double *r_numa = (double *) pattern_bind_memory_numa(r_pages, sockets, mem_pattern);
    cilk_for(int i =0;i < NA+2+1; i++){
      r_numa[i] = r[i];
    }
    //unsigned long midt1 = (lastcol-firstcol+2)/2;
    //unsigned long midt11 = (1 + midt1)/2;
    //unsigned long midt12 = (midt1 + lastcol-firstcol+2)/2;
    unsigned long unit =  (lastcol-firstcol+2)/sockets;
#endif
    for (it = 1; it <= 1; it++) {
      //printf("in first iter\n");
/*--------------------------------------------------------------------
c  The call to the conjugate gradient routine:
c-------------------------------------------------------------------*/
#ifndef NO_PIN
      conj_grad (colidx_numa, rowstr_numa, x_numa, z_numa, a_numa, p_numa, q_numa, r_numa,/* w,*/ &rnorm);
#else
	conj_grad (colidx, rowstr, x, z, a, p, q, r,/* w,*/ &rnorm);
#endif

	//printf("end conj grad\n");
/*--------------------------------------------------------------------
c  zeta = shift + 1/(x.z)
c  So, first: (x.z)
c  Also, find norm of z
c  So, first: (z.z)
c-------------------------------------------------------------------*/
	norm_temp11 = 0.0;
	norm_temp12 = 0.0;
#ifndef NO_PIN
	double norm_temp11_s[sockets];
	double norm_temp12_s[sockets];

	__cilkrts_disable_nonlocal_steal();
	int start_spawn0 = 1;
	for(int i = 0; i < sockets; i++){
	  if(i != sockets - 1){
	    norm_temp11_s[i] = cilk_spawn reduce_add_mul(norm_temp11, x_numa, z_numa, start_spawn0, start_spawn0+unit); // @0
	    SET_PIN(pin_pattern[i+1]);
	    norm_temp12_s[i] = cilk_spawn reduce_add_mul(norm_temp12, z_numa, z_numa, start_spawn0, start_spawn0+unit); // @0
	  }else{
	    norm_temp11_s[i] = cilk_spawn reduce_add_mul(norm_temp11, x_numa, z_numa, start_spawn0, lastcol-firstcol+2); // @3
	    __cilkrts_enable_nonlocal_steal();
	    norm_temp12_s[i] = reduce_add_mul(norm_temp12, z_numa, z_numa, start_spawn0, lastcol-firstcol+2); // @3
	  }
	  start_spawn0 += unit;
	}
        SET_PIN(pin_pattern[0]);
	cilk_sync;
	for(int i = 0; i < sockets; i ++){
	  norm_temp11 += norm_temp11_s[i];
	  norm_temp12 += norm_temp12_s[i];
	}
#else
	norm_temp11 = cilk_spawn reduce_add_mul(norm_temp11, x, z, 1, lastcol-firstcol+2);
	norm_temp12 = reduce_add_mul(norm_temp12, z, z, 1, lastcol-firstcol+2);
    cilk_sync;
#endif

	norm_temp12 = 1.0 / sqrt( norm_temp12 );

/*--------------------------------------------------------------------
c  Normalize z to obtain x
c-------------------------------------------------------------------*/
	cilk_for (unsigned long j = 1; j <= lastcol-firstcol+1; j++) {
#ifndef NO_PIN
	x_numa[j] = norm_temp12*z_numa[j];
#else
	x[j] = norm_temp12*z[j];
#endif
	}

    } /* end of do one iteration untimed */

/*--------------------------------------------------------------------
c  set starting vector to (1, 1, .... 1)
c-------------------------------------------------------------------*/
    cilk_for (unsigned long i = 1; i <= NA+1; i++) {
#ifndef NO_PIN
        x_numa[i] = 1.0;
#else
				x[i] = 1.0;
#endif
    }
    zeta  = 0.0;



/*--------------------------------------------------------------------
c---->
c  Main Iteration for inverse power method
c---->
c-------------------------------------------------------------------*/

    //printf("entering main iter\n");
    __cilkrts_reset_timing();
    timer_clear( 1 );
    //timer_start( 1 );
#if TIMING_COUNT
    begin = ktiming_getmark();
#endif
    for (it = 1; it <= NITER; it++) {
/*--------------------------------------------------------------------
c  The call to the conjugate gradient routine:
c-------------------------------------------------------------------*/
#ifndef NO_PIN
	conj_grad(colidx_numa, rowstr_numa, x_numa, z_numa, a_numa, p_numa, q_numa, r_numa/*, w*/, &rnorm);
#else
	conj_grad(colidx, rowstr, x, z, a, p, q, r/*, w*/, &rnorm);
#endif

/*--------------------------------------------------------------------
c  zeta = shift + 1/(x.z)
c  So, first: (x.z)
c  Also, find norm of z
c  So, first: (z.z)
c-------------------------------------------------------------------*/
	norm_temp11 = 0.0;
	norm_temp12 = 0.0;
#ifndef NO_PIN
	double norm_temp11_s[sockets];
	double norm_temp12_s[sockets];

	__cilkrts_disable_nonlocal_steal();
	int start_spawn1 = 1;
	for(int i = 0; i < sockets; i++){
	  if(i != sockets - 1){
	    norm_temp11_s[i] = cilk_spawn reduce_add_mul(norm_temp11, x_numa, z_numa, start_spawn1, start_spawn1+unit); // @0
	    SET_PIN(pin_pattern[i+1]);
	    norm_temp12_s[i] = cilk_spawn reduce_add_mul(norm_temp12, z_numa, z_numa, start_spawn1, start_spawn1+unit); // @0
	  }else{
	    norm_temp11_s[i] = cilk_spawn reduce_add_mul(norm_temp11, x_numa, z_numa, start_spawn1, lastcol-firstcol+2); // @3
	    __cilkrts_enable_nonlocal_steal();
	    norm_temp12_s[i] = reduce_add_mul(norm_temp12, z_numa, z_numa, start_spawn1, lastcol-firstcol+2); // @3
	  }
	  start_spawn1 += unit;
	}
        SET_PIN(pin_pattern[0]);
	cilk_sync;
	for(int i = 0; i < sockets; i ++){
	  norm_temp11 += norm_temp11_s[i];
	  norm_temp12 += norm_temp12_s[i];
	}
#else
	norm_temp11 = cilk_spawn reduce_add_mul(norm_temp11, x, z, 1, lastcol-firstcol+2);
	norm_temp12 = reduce_add_mul(norm_temp12, z, z, 1, lastcol-firstcol+2);
    cilk_sync;
#endif
	norm_temp12 = 1.0 / sqrt( norm_temp12 );

	zeta = SHIFT + 1.0 / norm_temp11;

   /*
	if( it == 1 ) {
	  printf("   iteration           ||r||                 zeta\n");
	}
	printf("    %5d       %20.14e%20.13e\n", it, rnorm, zeta);
    */
/*--------------------------------------------------------------------
c  Normalize z to obtain x
c-------------------------------------------------------------------*/
	/* cilk_for (unsigned long j = 1; j <= lastcol-firstcol+1; j++) { */
        /*     x[j] = norm_temp12*z[j]; */
	/* } */
#ifndef NO_PIN
	__cilkrts_disable_nonlocal_steal();
	int start_spawn2 = 1;
	for(int i = 0; i < sockets; i++){
	  if(i != sockets - 1){
	    SET_PIN(pin_pattern[i+1]);
	    cilk_spawn map_mul_scalar(x_numa, z_numa, norm_temp12, start_spawn2, start_spawn2+unit); // @0
	  }else{
	    __cilkrts_enable_nonlocal_steal();
	    map_mul_scalar(x_numa, z_numa, norm_temp12, start_spawn2, lastcol-firstcol+2); // @3
	  }
	  start_spawn2 += unit;
	}
        SET_PIN(pin_pattern[0]);
	cilk_sync;
#else
	cilk_for (unsigned long j = 1; j <= lastcol-firstcol+1; j++) {
		x[j] = norm_temp12*z[j];
	}
#endif
}
		/* end of main iter inv pow meth */
    nthreads = __cilkrts_get_nworkers();
    //timer_stop( 1 );
#if TIMING_COUNT
    end = ktiming_getmark();
    elapsed[run] = ktiming_diff_usec(&begin, &end);
#endif
    __cilkrts_accum_timing();
/*--------------------------------------------------------------------
c  End of timed section
c-------------------------------------------------------------------*/

    t = timer_read( 1 );

    //printf(" Benchmark completed\n");

    epsilon = 1.0e-10;
    if (class != 'U') {
	if (fabs(zeta - zeta_verify_value) <= epsilon) {
            verified = TRUE;
	    //printf(" VERIFICATION SUCCESSFUL\n");
	    //printf(" Zeta is    %20.12e\n", zeta);
	    //printf(" Error is   %20.12e\n", zeta - zeta_verify_value);
	} else {
            verified = FALSE;
	    //printf(" VERIFICATION FAILED\n");
	    //printf(" Zeta                %20.12e\n", zeta);
	    //printf(" The correct zeta is %20.12e\n", zeta_verify_value);
	}
    } else {
	verified = FALSE;
	//printf(" Problem size unknown\n");
	//printf(" NO VERIFICATION PERFORMED\n");
    }

    if ( t != 0.0 ) {
	mflops = (2.0*NITER*NA)
	    * (3.0+(NONZER*(NONZER+1)) + 25.0*(5.0+(NONZER*(NONZER+1))) + 3.0 )
	    / t / 1000000.0;
    } else {
	mflops = 0.0;
    }

    c_print_results("CG", class, NA, 0, 0, NITER, nthreads, t,
		    mflops, "          floating point",
		    verified, NPBVERSION, COMPILETIME,
		    CS1, CS2, CS3, CS4, CS5, CS6, CS7);

    return 0;
}
const unsigned long INIT_BASE_CASE = BASE_CASE;

void initialize(double *q, double *z, double *r, double *p, double *x, unsigned long i, unsigned long e){
  unsigned long mid = (i + e) / 2;
  if(e - i  <= INIT_BASE_CASE){
    for(int j = i; j < e; j++){
      q[j] = 0.0;
      z[j] = 0.0;
      r[j] = x[j];
      p[j] = r[j];
    }
    return;
  }
  cilk_spawn initialize(q, z, r, p, x, i, mid);
  initialize(q, z, r, p, x, mid, e);
  cilk_sync;
}
      /* cilk_for (unsigned long j = 1; j <= lastrow-firstrow+1; j++) { */
      /* 	sum = 0.0; */
      /* 	for (k = rowstr[j]; k < rowstr[j+1]; k++) { */
      /* 	  sum = sum + a[k]*p[colidx[k]]; */
      /* 	} */
      /* 	//w[j] = sum; */
      /* 	q[j] = sum; */
      /* } */
const unsigned long Q_BASE_CASE = 1;
void compute_q(int *rowstr, int *colidx, double *p, double *q, double *a, unsigned long i, unsigned long e){

  if(e - i >= NA/4-100){
    //start =  ktiming_getmark();
  }
  unsigned long mid = (i+e)/2;
  if(e - i <= Q_BASE_CASE){
    for(int j = i; j < e; j++){
      double sum = 0.0;
      //printf("j: %d, rowstr: %d, rowstr+1: %d, work: %d\n", j, rowstr[j], rowstr[j+1], rowstr[j+1]-rowstr[j]);
      /*
      for (int k = rowstr[j]; k < rowstr[j+1]; k++){
      	//printf("k: %d, colidx: %d, p: %f, NA: %d\n", k, colidx[k], p[colidx[k]], NA);
      	sum += a[k]* p[colidx[k]];
      }*/
/*#ifndef NO_PIN
      sum  = reduce_add_mul_idx(sum, a, p, colidx, rowstr[j], rowstr[j+1]);
#else*/
			for (int k = rowstr[j]; k < rowstr[j+1]; k++){
				sum += a[k] * p[colidx[k]];
			}
//#endif

      q[j] = sum;
    }
    return;
  }
  cilk_spawn compute_q(rowstr, colidx, p, q, a, i, mid);
  compute_q(rowstr, colidx, p, q, a, mid, e);
  cilk_sync;

  if(e - i >= NA/4-100){
    //clockmark_t end =  ktiming_getmark();
    //uint64_t diff = ktiming_diff_usec(&start, &end);
    //printf("NZ: %d, diff : %lu at i: %d, e: %d, start:%d, end: %d, amount of work: %d, a bind at %d, i am at %d\n", NZ, diff, i, e, rowstr[i], rowstr[e], rowstr[e] - rowstr[i], get_mem_binding(&a[(rowstr[e] + rowstr[i])/2]), numa_node_of_cpu(sched_getcpu()));
  }

}
const unsigned long MAP_BASE_CASE = BASE_CASE;
void map_add_mul(double *res, double *A, double *B, double alpha, unsigned long i, unsigned long e){
  unsigned long mid = (i + e) / 2;
  if(e - i < MAP_BASE_CASE){
    for(int j = i; j < e; j++){
      res[j] = A[j] + B[j] * alpha;
    }
    return;
  }

  cilk_spawn map_add_mul(res, A, B, alpha, i, mid);
  map_add_mul(res, A, B, alpha, mid, e);
  cilk_sync;
}

/*#ifdef NO_PIN
const unsigned long IDX_MUL_BASE_CASE = 512;
double reduce_add_mul_idx(double init_sum, double *arr1, double *arr2, int *idx, unsigned long i, unsigned long e){
  unsigned long mid = (i + e) / 2;
  //printf("here: e- i : %d\n", e-i);
  if( e - i < IDX_MUL_BASE_CASE){
    for(int j = i; j < e; j++){
      init_sum += arr1[j] * arr2[idx[j]];
    }
    return init_sum;
  }
  double k = cilk_spawn reduce_add_mul_idx(init_sum, arr1, arr2, idx, i, mid);
  double l = reduce_add_mul_idx(init_sum, arr1, arr2, idx, mid, e);
  cilk_sync;
  return k+l;
}
#endif*/

    /* for (unsigned long j = 1; j <= lastrow-firstrow+1; j++) { */
    /* 	d = 0.0; */
    /* 	for (k = rowstr[j]; k <= rowstr[j+1]-1; k++) { */
    /*         d = d + a[k]*z[colidx[k]]; */
    /* 	} */
    /* 	r[j] = d; */
    /* } */

const unsigned long NORM_BASE_CASE = 1;
void compute_norm(int *rowstr, int *colidx, double *a, double *z, double *r, unsigned long i, unsigned long e){
  unsigned long mid = (i+e)/2;
  if(e - i <= NORM_BASE_CASE){
    for(int j = i; j < e; j++){
      r[j] = reduce_add_mul_idx(0.0, a, z, colidx, rowstr[j], rowstr[j+1]);
    }
    return;
  }
  cilk_spawn compute_norm(rowstr, colidx, a, z, r, i, mid);
  compute_norm(rowstr, colidx, a, z, r, mid, e);
  cilk_sync;
}
#ifndef NO_PIN
long find_num(int *arr, int num, unsigned long i, unsigned long e){
  unsigned long mid = (i + e)/2;
  if(e - i < BASE_CASE){
    long idx = -1;
    for(int j = i; j < e; j++){
      if (arr[j] > num){
	idx= j-1;
	break;
      }
    }
    return idx;
  }
  long res1= cilk_spawn find_num(arr, num, i, mid);
  long res2= find_num(arr, num, mid, e);
  cilk_sync;
  if(res1 != -1){
    return res1;
  }
  if(res2 != -1){
    return res2;
  }
  return -1;
}
#endif

/* #pragma omp for reduction(+:sum) */
/*     for (j = 1; j <= lastcol-firstcol+1; j++) { */
/* 	d = x[j] - r[j]; */
/* 	sum = sum + d*d; */
/*     } */
/*     //} //end omp parallel */
/*     (*rnorm) = sqrt(sum); */
/* } */
const unsigned long COMPUTE_SUM_BASE_CASE = BASE_CASE;
double compute_sum(double *x, double *r, unsigned long i, unsigned long e){
  unsigned long mid = (i + e) / 2;
  if(e - i < COMPUTE_SUM_BASE_CASE){
    double sum = 0.0;
    for(unsigned long j = i; j < e; j++){
      double d = x[j] - r[j];
      sum += d*d;
    }
    return sum;
  }
  double v1 = cilk_spawn compute_sum(x, r, i, mid);
  double v2 = compute_sum(x, r, mid, e);
  cilk_sync;
  return v1 + v2;
}
/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/
static void conj_grad (
    int colidx[],	/* colidx[1:nzz] */
    int rowstr[],	/* rowstr[1:naa+1] */
    double x[],		/* x[*] */
    double z[],		/* z[*] */
    double a[],		/* a[1:nzz] */
    double p[],		/* p[*] */
    double q[],		/* q[*] */
    double r[],		/* r[*] */
    //double w[],		/* w[*] */
    double *rnorm )
/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*---------------------------------------------------------------------
c  Floaging point arrays here are named as in NPB1 spec discussion of
c  CG algorithm
c---------------------------------------------------------------------*/
{
    static int callcount = 0;
    double d, sum, rho, rho0, alpha, beta;
    int i, j, k;
    int cgit, cgitmax = 25;

    rho = 0.0;
    //#pragma omp parallel default(shared) private(j,sum) shared(rho,naa)
/*--------------------------------------------------------------------
c  Initialize the CG algorithm:
c-------------------------------------------------------------------*/
//{
    /* cilk_for (unsigned long j = 1; j <= naa+1; j++) { */
    /* 	q[j] = 0.0; */
    /* 	z[j] = 0.0; */
    /* 	r[j] = x[j]; */
    /* 	p[j] = r[j]; */
    /* 	//w[j] = 0.0; */
    /* } */
#ifndef NO_PIN
    int sockets = __cilkrts_num_sockets();
    unsigned long init_unit = (1 + naa+2)/4;
    __cilkrts_disable_nonlocal_steal();
    int start_spawn1 = 1;
    for(int i = 0; i < sockets; i++){
      if(i != sockets - 1){
	SET_PIN(pin_pattern[i+1]);
	cilk_spawn initialize(q, z, r, p, x, start_spawn1, start_spawn1+init_unit); // @0
      }else{
	__cilkrts_enable_nonlocal_steal();
	initialize(q, z, r, p, x, start_spawn1, naa+2); // @3
      }
      start_spawn1 += init_unit;
    }
    SET_PIN(pin_pattern[0]);
    cilk_sync;
#else
	initialize(q, z, r, p, x, 1, naa+2);
#endif
/*--------------------------------------------------------------------
c  rho = r.r
c  Now, obtain the norm of r: First, sum squares of r elements locally...
c-------------------------------------------------------------------*/
/* #pragma omp for reduction(+:rho) */
/*     for (j = 1; j <= lastcol-firstcol+1; j++) { */
/* 	rho = rho + r[j]*r[j]; */
/*     } */
#ifndef NO_PIN
    unsigned long unit =  (lastcol-firstcol+2)/sockets;
    __cilkrts_disable_nonlocal_steal();
    int start_spawn2 = 1;
    double rhos[sockets];

    for(int i = 0; i < sockets; i++){
      if(i != sockets - 1){
	SET_PIN(pin_pattern[i+1]);
	rhos[i] = cilk_spawn reduce_add_mul(rho, r, r, start_spawn2, start_spawn2+unit);
      }else{
	__cilkrts_enable_nonlocal_steal();
	rhos[i] = reduce_add_mul(rho, r, r, start_spawn2, lastrow-firstrow+2);
      }
      start_spawn2 += unit;
    }
    SET_PIN(pin_pattern[0]);
    cilk_sync;
    for (int i = 0; i < sockets; i++){
      rho += rhos[i];
    }
#else
		rho = reduce_add_mul(rho, r, r, 1, lastcol-firstcol+2);
#endif

    //}/* end omp parallel */
/*--------------------------------------------------------------------
c---->
c  The conj grad iteration loop
c---->
c-------------------------------------------------------------------*/

    for (cgit = 1; cgit <= cgitmax; cgit++) {
      rho0 = rho;
      d = 0.0;
      rho = 0.0;
      int last_ele = rowstr[lastrow-firstrow+2];
      //int rmid1 = cilk_spawn find_num(rowstr, last_ele/4, 0, lastrow-firstrow+2);
      //int rmid = cilk_spawn find_num(rowstr, last_ele/2, 0, lastrow-firstrow+2);
      //int rmid2 = find_num(rowstr, last_ele/4*3, 0, lastrow-firstrow+2);
      //cilk_sync;
      //printf("rmid1: %d, rmid: %d, rmid2: %d, last: %d\n", rmid1, rmid, rmid2, lastrow-firstrow+2);
      //printf("%d, %d, %d, %d\n", rowstr[rmid1] - rowstr[1], rowstr[rmid] - rowstr[rmid1], rowstr[rmid2] - rowstr[rmid], last_ele - rowstr[rmid2]);
      //printf("prev %d, %d, %d, %d\n", rowstr[mid1] - rowstr[1], rowstr[mid] - rowstr[mid1], rowstr[mid2] - rowstr[mid], last_ele - rowstr[mid2]);
#ifndef NO_PIN
      __cilkrts_disable_nonlocal_steal();
      int start_spawn3 = 1;
      for(int i = 0; i < sockets; i++){
	if(i != sockets - 1){
	  SET_PIN(pin_pattern[i+1]);
	  cilk_spawn compute_q(rowstr, colidx, p, q, a, start_spawn3, start_spawn3+unit);
	}else{
	  __cilkrts_enable_nonlocal_steal();
	  compute_q(rowstr, colidx, p, q, a, start_spawn3, lastrow-firstrow+2);
	}
	start_spawn3 += unit;
      }
      SET_PIN(pin_pattern[0]);
      cilk_sync;
#else
      compute_q(rowstr, colidx, p, q, a, 1, lastrow-firstrow+2);
#endif

/*--------------------------------------------------------------------
c  Clear w for reuse...
c-------------------------------------------------------------------*/
/*
#pragma omp for	nowait
	for (j = 1; j <= lastcol-firstcol+1; j++) {
            w[j] = 0.0;
	}
*/
/*--------------------------------------------------------------------
c  Obtain p.q
c-------------------------------------------------------------------*/
/* #pragma omp for reduction(+:d) */
/* 	for (j = 1; j <= lastcol-firstcol+1; j++) { */
/*             d = d + p[j]*q[j]; */
/* 	} */
/* #pragma omp barrier */
#ifndef NO_PIN
      double ds[sockets];
      __cilkrts_disable_nonlocal_steal();
      int start_spawn4 = 1;
      for(int i = 0; i < sockets; i++){
	if(i != sockets - 1){
	  SET_PIN(pin_pattern[i+1]);
	  ds[i] = cilk_spawn reduce_add_mul(d, p, q, start_spawn4, start_spawn4+unit);
	}else{
	  __cilkrts_enable_nonlocal_steal();
	  ds[i] = reduce_add_mul(d, p, q, start_spawn4, lastcol-firstcol+2);
	}
	start_spawn4 += unit;
      }
      SET_PIN(pin_pattern[0]);
      cilk_sync;
      for(int i = 0; i < sockets; i++){
	d+= ds[i];
      }
#else
	d = reduce_add_mul(d, p, q, 1, lastcol-firstcol+2);
#endif
/*--------------------------------------------------------------------
c  Obtain alpha = rho / (p.q)
c-------------------------------------------------------------------*/
//#pragma omp single
	alpha = rho0 / d;

/*--------------------------------------------------------------------
c  Save a temporary of rho
c-------------------------------------------------------------------*/
	/*	rho0 = rho;*/

/*---------------------------------------------------------------------
c  Obtain z = z + alpha*p
c  and    r = r - alpha*q
c---------------------------------------------------------------------*/
//#pragma omp for reduction(+:rho)
	//remember we can actually combine the operations here
	/* cilk_for (unsigned long j = 1; j <= lastcol-firstcol+1; j++) { */
        /*     z[j] = z[j] + alpha*p[j]; */
        /*     r[j] = r[j] - alpha*q[j]; */
	/* } */
#ifndef NO_PIN
	__cilkrts_disable_nonlocal_steal();
	int start_spawn5 = 1;
	for(int i = 0; i < sockets; i++){
	if(i != sockets - 1){
	  cilk_spawn map_add_mul(z, z, p, alpha, start_spawn5, start_spawn5+unit); // @0
	  SET_PIN(pin_pattern[i+1]);
	  cilk_spawn map_add_mul(r, r, q, -alpha, start_spawn5, start_spawn5+unit); // @0
	}else{
	  cilk_spawn map_add_mul(z, z, p, alpha, start_spawn5, lastcol-firstcol+2); // @3
	  __cilkrts_enable_nonlocal_steal();
	  map_add_mul(r, r, q, -alpha, start_spawn5, lastcol-firstcol+2); // @3
	}
	start_spawn5 += unit;
      }
      SET_PIN(pin_pattern[0]);
#else
	cilk_spawn map_add_mul(z, z, p, alpha, 1, lastcol-firstcol+2);
	map_add_mul(r, r, q, -alpha, 1, lastcol-firstcol+2);
#endif
      cilk_sync;

/*---------------------------------------------------------------------
c  rho = r.r
c  Now, obtain the norm of r: First, sum squares of r elements locally...
c---------------------------------------------------------------------*/
/*
#pragma omp for
	for (j = 1; j <= lastcol-firstcol+1; j++) {*/
	//rho = rho + r[j]*r[j];
	//}
//#pragma omp barrier
#ifndef NO_PIN
	double rhos[sockets];
	__cilkrts_disable_nonlocal_steal();
	int start_spawn6 = 1;
	for(int i = 0; i < sockets; i++){
	  if(i != sockets - 1){
	    SET_PIN(pin_pattern[i+1]);
	    rhos[i] = cilk_spawn reduce_add_mul(rho, r, r, start_spawn6, start_spawn6+unit);
	  }else{
	    __cilkrts_enable_nonlocal_steal();
	    rhos[i] = cilk_spawn reduce_add_mul(rho, r, r, start_spawn6, lastcol-firstcol+2);
	  }
	  start_spawn6 += unit;
	}
	SET_PIN(pin_pattern[0]);
	cilk_sync;
	for(int i = 0; i < sockets; i++){
	  rho += rhos[i];
	}
#else
	rho = reduce_add_mul(rho, r, r, 1, lastcol-firstcol+2);
#endif

	/*--------------------------------------------------------------------
c  Obtain beta:
c-------------------------------------------------------------------*/
//#pragma omp single
	beta = rho / rho0;

/*--------------------------------------------------------------------
c  p = r + beta*p
c-------------------------------------------------------------------*/
//#pragma omp for nowait
	/* cilk_for (unsigned long j = 1; j <= lastcol-firstcol+1; j++) { */
        /*     p[j] = r[j] + beta*p[j]; */
	/* } */
#ifndef NO_PIN
	__cilkrts_disable_nonlocal_steal();
	int start_spawn7 = 1;
	for(int i = 0; i < sockets; i++){
	  if(i != sockets - 1){
	    SET_PIN(pin_pattern[i+1]);
	    cilk_spawn map_add_mul(p, r, p, beta, start_spawn7, start_spawn7+unit);
	  }else{
	    __cilkrts_enable_nonlocal_steal();
	    map_add_mul(p, r, p, beta, start_spawn7, lastcol-firstcol+2);
	  }
	  start_spawn7 += unit;
	}
	SET_PIN(pin_pattern[0]);
	cilk_sync;
#else
	map_add_mul(p, r, p, beta, 1, lastcol-firstcol+2);
#endif
    callcount++;
    //} /* end omp parallel */
} /* end of do cgit=1,cgitmax */

/*---------------------------------------------------------------------
c  Compute residual norm explicitly:  ||r|| = ||x - A.z||
c  First, form A.z
c  The partition submatrix-vector multiply
c---------------------------------------------------------------------*/
    sum = 0.0;

    //#pragma omp parallel default(shared) private(j,d) shared(sum)
    //{
    //#pragma omp for //private(d, k)
    /* for (unsigned long j = 1; j <= lastrow-firstrow+1; j++) { */
    /* 	d = 0.0; */
    /* 	for (k = rowstr[j]; k <= rowstr[j+1]-1; k++) { */
    /*         d = d + a[k]*z[colidx[k]]; */
    /* 	} */
    /* 	r[j] = d; */
    /* } */
	#ifndef NO_PIN
    __cilkrts_disable_nonlocal_steal();
    int start_spawn8 = 1;
    for(int i = 0; i < sockets; i++){
      if(i != sockets - 1){
	SET_PIN(pin_pattern[i+1]);
	cilk_spawn compute_norm(rowstr, colidx, a, z, r, start_spawn8, start_spawn8+unit);
      }else{
	__cilkrts_enable_nonlocal_steal();
	compute_norm(rowstr, colidx, a, z, r, start_spawn8, lastrow-firstrow+2);
      }
      start_spawn8 += unit;
    }
    SET_PIN(pin_pattern[0]);
    cilk_sync;
#else
	compute_norm(rowstr, colidx, a, z, r, 1, lastrow-firstrow+2);
#endif

/*--------------------------------------------------------------------
c  At this point, r contains A.z
c-------------------------------------------------------------------*/
//todo: finish the reduction here
/* #pragma omp for reduction(+:sum) */
/*     for (j = 1; j <= lastcol-firstcol+1; j++) { */
/* 	d = x[j] - r[j]; */
/* 	sum = sum + d*d; */
/*     } */
/*     //} //end omp parallel */
/*     (*rnorm) = sqrt(sum); */
/* } */
#ifndef NO_PIN
    double sums[sockets];
    __cilkrts_disable_nonlocal_steal();
    int start_spawn9 = 1;
    for(int i = 0; i < sockets; i++){
      if(i != sockets - 1){
	SET_PIN(pin_pattern[i+1]);
	sums[i] = cilk_spawn compute_sum(x, r, start_spawn9, start_spawn9+unit);
      }else{
	__cilkrts_enable_nonlocal_steal();
	sums[i] = compute_sum(x, r, start_spawn9, lastcol-firstcol+2);
      }
      start_spawn9 += unit;
    }
    SET_PIN(pin_pattern[0]);
    cilk_sync;
    for(int i = 0; i < sockets; i++){
      sum += sums[i];
    }
#else
	sum = compute_sum(x, r, 1, lastcol-firstcol+2);
#endif

    (*rnorm) = sqrt(sum);
}
/*---------------------------------------------------------------------
c       generate the test problem for benchmark 6
c       makea generates a sparse matrix with a
c       prescribed sparsity distribution
c
c       parameter    type        usage
c
c       input
c
c       n            i           number of cols/rows of matrix
c       nz           i           nonzeros as declared array size
c       rcond        r*8         condition number
c       shift        r*8         main diagonal shift
c
c       output
c
c       a            r*8         array for nonzeros
c       colidx       i           col indices
c       rowstr       i           row pointers
c
c       workspace
c
c       iv, arow, acol i
c       v, aelt        r*8
c---------------------------------------------------------------------*/
static void makea(
    int n,
    int nz,
    double a[],		/* a[1:nz] */
    int colidx[],	/* colidx[1:nz] */
    int rowstr[],	/* rowstr[1:n+1] */
    int nonzer,
    int firstrow,
    int lastrow,
    int firstcol,
    int lastcol,
    double rcond,
    int arow[],		/* arow[1:nz] */
    int acol[],		/* acol[1:nz] */
    double aelt[],	/* aelt[1:nz] */
    double v[],		/* v[1:n+1] */
    int iv[],		/* iv[1:2*n+1] */
    double shift )
{
    int i, nnza, iouter, ivelt, ivelt1, irow, nzv;

/*--------------------------------------------------------------------
c      nonzer is approximately  (int(sqrt(nnza /n)));
c-------------------------------------------------------------------*/

    double size, ratio, scale;
    int jcol;

    size = 1.0;
    ratio = pow(rcond, (1.0 / (double)n));
    nnza = 0;

/*---------------------------------------------------------------------
c  Initialize colidx(n+1 .. 2n) to zero.
c  Used by sprnvc to mark nonzero positions
c---------------------------------------------------------------------*/
    //cilkfor
    for (unsigned long i = 1; i <= n; i++) {
	colidx[n+i] = 0;
    }
    for (iouter = 1; iouter <= n; iouter++) {
	nzv = nonzer;
	sprnvc(n, nzv, v, iv, &(colidx[0]), &(colidx[n]));
	vecset(n, v, iv, &nzv, iouter, 0.5);
	for (ivelt = 1; ivelt <= nzv; ivelt++) {
	    jcol = iv[ivelt];
	    if (jcol >= firstcol && jcol <= lastcol) {
		scale = size * v[ivelt];
		for (ivelt1 = 1; ivelt1 <= nzv; ivelt1++) {
	            irow = iv[ivelt1];
                    if (irow >= firstrow && irow <= lastrow) {
			nnza = nnza + 1;
			if (nnza > nz) {
			    //printf("Space for matrix elements exceeded in"
				  // " makea\n");
			    //printf("nnza, nzmax = %d, %d\n", nnza, nz);
			    //printf("iouter = %d\n", iouter);
			    exit(1);
			}
			acol[nnza] = jcol;
			arow[nnza] = irow;
			aelt[nnza] = v[ivelt1] * scale;
		    }
		}
	    }
	}
	size = size * ratio;
    }

/*---------------------------------------------------------------------
c       ... add the identity * rcond to the generated matrix to bound
c           the smallest eigenvalue from below by rcond
c---------------------------------------------------------------------*/
    for (i = firstrow; i <= lastrow; i++) {
	if (i >= firstcol && i <= lastcol) {
	    iouter = n + i;
	    nnza = nnza + 1;
	    if (nnza > nz) {
		//printf("Space for matrix elements exceeded in makea\n");
		//printf("nnza, nzmax = %d, %d\n", nnza, nz);
		//printf("iouter = %d\n", iouter);
		exit(1);
	    }
	    acol[nnza] = i;
	    arow[nnza] = i;
	    aelt[nnza] = rcond - shift;
	}
    }

/*---------------------------------------------------------------------
c       ... make the sparse matrix from list of elements with duplicates
c           (v and iv are used as  workspace)
c---------------------------------------------------------------------*/
    sparse(a, colidx, rowstr, n, arow, acol, aelt,
	   firstrow, lastrow, v, &(iv[0]), &(iv[n]), nnza);
}

/*---------------------------------------------------
c       generate a sparse matrix from a list of
c       [col, row, element] tri
c---------------------------------------------------*/
static void sparse(
    double a[],		/* a[1:*] */
    int colidx[],	/* colidx[1:*] */
    int rowstr[],	/* rowstr[1:*] */
    int n,
    int arow[],		/* arow[1:*] */
    int acol[],		/* acol[1:*] */
    double aelt[],	/* aelt[1:*] */
    int firstrow,
    int lastrow,
    double x[],		/* x[1:n] */
    boolean mark[],	/* mark[1:n] */
    int nzloc[],	/* nzloc[1:n] */
    int nnza)
/*---------------------------------------------------------------------
c       rows range from firstrow to lastrow
c       the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
c---------------------------------------------------------------------*/
{
    int nrows;
    int i, j, jajp1, nza, k, nzrow;
    double xi;

/*--------------------------------------------------------------------
c    how many rows of result
c-------------------------------------------------------------------*/
    nrows = lastrow - firstrow + 1;

/*--------------------------------------------------------------------
c     ...count the number of triples in each row
c-------------------------------------------------------------------*/
    //cilk_for
    for (unsigned long j = 1; j <= n; j++) {
	rowstr[j] = 0;
	mark[j] = FALSE;
    }
    rowstr[n+1] = 0;

    for (nza = 1; nza <= nnza; nza++) {
	j = (arow[nza] - firstrow + 1) + 1;
	rowstr[j] = rowstr[j] + 1;
    }

    rowstr[1] = 1;
    for (j = 2; j <= nrows+1; j++) {
	rowstr[j] = rowstr[j] + rowstr[j-1];
    }

/*---------------------------------------------------------------------
c     ... rowstr(j) now is the location of the first nonzero
c           of row j of a
c---------------------------------------------------------------------*/

/*---------------------------------------------------------------------
c     ... preload data pages
c---------------------------------------------------------------------*/
    //cilkfor
    for(unsigned long j = 0;j <= nrows-1;j++) {
      //cilkfor
      for(unsigned long k = rowstr[j];k <= rowstr[j+1]-1;k++)
	       a[k] = 0.0;
      }
/*--------------------------------------------------------------------
c     ... do a bucket sort of the triples on the row index
c-------------------------------------------------------------------*/
    for (nza = 1; nza <= nnza; nza++) {
	j = arow[nza] - firstrow + 1;
	k = rowstr[j];
	a[k] = aelt[nza];
	colidx[k] = acol[nza];
	rowstr[j] = rowstr[j] + 1;
    }

/*--------------------------------------------------------------------
c       ... rowstr(j) now points to the first element of row j+1
c-------------------------------------------------------------------*/
    for (j = nrows; j >= 1; j--) {
	rowstr[j+1] = rowstr[j];
    }
    rowstr[1] = 1;

/*--------------------------------------------------------------------
c       ... generate the actual output rows by adding elements
c-------------------------------------------------------------------*/
    nza = 0;
    //cilkfor
    for (unsigned long i = 1; i <= n; i++) {
	x[i] = 0.0;
	mark[i] = FALSE;
    }

    jajp1 = rowstr[1];
    for (j = 1; j <= nrows; j++) {
	nzrow = 0;

/*--------------------------------------------------------------------
c          ...loop over the jth row of a
c-------------------------------------------------------------------*/
	for (k = jajp1; k < rowstr[j+1]; k++) {
            i = colidx[k];
            x[i] = x[i] + a[k];
            if ( mark[i] == FALSE && x[i] != 0.0) {
		mark[i] = TRUE;
		nzrow = nzrow + 1;
		nzloc[nzrow] = i;
	    }
	}

/*--------------------------------------------------------------------
c          ... extract the nonzeros of this row
c-------------------------------------------------------------------*/
	for (k = 1; k <= nzrow; k++) {
            i = nzloc[k];
            mark[i] = FALSE;
            xi = x[i];
            x[i] = 0.0;
            if (xi != 0.0) {
		nza = nza + 1;
		a[nza] = xi;
		colidx[nza] = i;
	    }
	}
	jajp1 = rowstr[j+1];
	rowstr[j+1] = nza + rowstr[1];
    }
}

/*---------------------------------------------------------------------
c       generate a sparse n-vector (v, iv)
c       having nzv nonzeros
c
c       mark(i) is set to 1 if position i is nonzero.
c       mark is all zero on entry and is reset to all zero before exit
c       this corrects a performance bug found by John G. Lewis, caused by
c       reinitialization of mark on every one of the n calls to sprnvc
---------------------------------------------------------------------*/
static void sprnvc(
    int n,
    int nz,
    double v[],		/* v[1:*] */
    int iv[],		/* iv[1:*] */
    int nzloc[],	/* nzloc[1:n] */
    int mark[] ) 	/* mark[1:n] */
{
    int nn1;
    int nzrow, nzv, ii, i;
    double vecelt, vecloc;

    nzv = 0;
    nzrow = 0;
    nn1 = 1;
    do {
	nn1 = 2 * nn1;
    } while (nn1 < n);

/*--------------------------------------------------------------------
c    nn1 is the smallest power of two not less than n
c-------------------------------------------------------------------*/

    while (nzv < nz) {
	vecelt = randlc(&tran, amult);

/*--------------------------------------------------------------------
c   generate an integer between 1 and n in a portable manner
c-------------------------------------------------------------------*/
	vecloc = randlc(&tran, amult);
	i = icnvrt(vecloc, nn1) + 1;
	if (i > n) continue;

/*--------------------------------------------------------------------
c  was this integer generated already?
c-------------------------------------------------------------------*/
	if (mark[i] == 0) {
	    mark[i] = 1;
	    nzrow = nzrow + 1;
	    nzloc[nzrow] = i;
	    nzv = nzv + 1;
	    v[nzv] = vecelt;
	    iv[nzv] = i;
	}
    }

    for (ii = 1; ii <= nzrow; ii++) {
	i = nzloc[ii];
	mark[i] = 0;
    }
}

/*---------------------------------------------------------------------
* scale a double precision number x in (0,1) by a power of 2 and chop it
*---------------------------------------------------------------------*/
static int icnvrt(double x, int ipwr2) {
    return ((int)(ipwr2 * x));
}

/*--------------------------------------------------------------------
c       set ith element of sparse vector (v, iv) with
c       nzv nonzeros to val
c-------------------------------------------------------------------*/
static void vecset(
    int n,
    double v[],	/* v[1:*] */
    int iv[],	/* iv[1:*] */
    int *nzv,
    int i,
    double val)
{
    int k;
    boolean set;

    set = FALSE;
    for (k = 1; k <= *nzv; k++) {
	if (iv[k] == i) {
            v[k] = val;
            set  = TRUE;
	}
    }
    if (set == FALSE) {
	*nzv = *nzv + 1;
	v[*nzv] = val;
	iv[*nzv] = i;
    }
}
int main(int argc, char **argv){
  __cilkrts_init();
#if TIMING_COUNT
    elapsed = malloc(TIMING_COUNT * sizeof(uint64_t));
    for(int i=0; i < TIMING_COUNT; i++){
    begin = 0;
    end = 0;
#else
    int i = 0;
#endif
  #ifndef NO_PIN
  __cilkrts_pin_top_level_frame_at_socket(0);
  sockets = __cilkrts_num_sockets();
  if(sockets == 2){
      int mem_patternT[] = {0, 1};
      int pin_patternT[] = {0, 1};
      memcpy(mem_pattern, mem_patternT, 4*sizeof(int));
      memcpy(pin_pattern, pin_patternT, 4*sizeof(int));
    }else if(sockets == 3){
      int mem_patternT[] = {0, 2, 1};
      int pin_patternT[] = {0, 2, 1};
      memcpy(mem_pattern, mem_patternT, 4*sizeof(int));
      memcpy(pin_pattern, pin_patternT, 4*sizeof(int));
    }else if(sockets == 4){
      int mem_patternT[] = {0, 1, 2, 3};
      int pin_patternT[] = {0, 1, 2, 3};
      memcpy(mem_pattern, mem_patternT, 4*sizeof(int));
      memcpy(pin_pattern, pin_patternT, 4*sizeof(int));
    }
  #else
    unsigned long nodemask = 0;
    for(int i = 0; i < __cilkrts_get_nworkers() / CPUS_PER_SOCKET; i++) {
      nodemask |= (1L << i);
    }
    set_mempolicy(MPOL_INTERLEAVE, &nodemask ,sizeof(nodemask)*8);
  #endif
  fakemain(argc, argv, i);
#if TIMING_COUNT
   }
   print_runtime(elapsed, TIMING_COUNT);
#endif
	__cilkrts_accum_timing();
  return 0;
}
