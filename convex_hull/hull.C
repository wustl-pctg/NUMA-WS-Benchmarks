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

#include <algorithm>
#include <cstring>
#include "parallel.h"
#include "geometry.h"
#include "sequence.h"
#include "gettime.h"
#include "ktiming.h"
#include <cilk/cilk_api.h>

#ifdef NO_PIN
#define __cilkrts_set_pinning_info(n)
#define __cilkrts_disable_nonlocal_steal()
#define __cilkrts_unset_pinning_info()
#define __cilkrts_enable_nonlocal_steal()
#define __cilkrts_pin_top_level_frame_at_socket(n)
#define SET_PIN(N)
#else
#define SET_PIN(N) __cilkrts_set_pinning_info(N)
#endif

#include <unistd.h>
#include <numa.h>
#include <stdlib.h>
using namespace std;
using namespace sequence;

template <class ET, class F>
pair<intT,intT> split(ET* A, intT n, F lf, F rf) {
  intT ll = 0, lm = 0;
  intT rm = n-1, rr = n-1;
  while (1) {
    while ((lm <= rm) && !(rf(A[lm]) > 0)) {
      if (lf(A[lm]) > 0) A[ll++] = A[lm];
      lm++;
    }
    while ((rm >= lm) && !(lf(A[rm]) > 0)) {
      if (rf(A[rm]) > 0) A[rr--] = A[rm];
      rm--;
    }
    if (lm >= rm) break; 
    ET tmp = A[lm++];
    A[ll++] = A[rm--];
    A[rr--] = tmp;
  }
  intT n1 = ll;
  intT n2 = n-rr-1;
  return pair<intT,intT>(n1,n2);
}

struct aboveLine {
  intT l, r;
  point2d* P;
  aboveLine(point2d* _P, intT _l, intT _r) : P(_P), l(_l), r(_r) {}
  bool operator() (intT i) {return triArea(P[l], P[r], P[i]) > 0.0;}
};


struct aboveLineP {
  point2d l, r;
  point2d* P;
  aboveLineP(point2d* _P, point2d &_l, point2d &_r) : P(_P), l(_l), r(_r) {}
  bool operator() (point2d &i) {return triArea(l, r, i) > 0.0;}
};


intT serialQuickHull(intT* I, point2d* P, intT n, intT l, intT r) {
  if (n < 2) return n;
  intT maxP = I[0];
  double maxArea = triArea(P[l],P[r],P[maxP]);
  for (intT i=1; i < n; i++) {
    intT j = I[i];
    double a = triArea(P[l],P[r],P[j]);
    if (a > maxArea) {
      maxArea = a;
      maxP = j;
    }
  }

  pair<intT,intT> nn = split(I, n, aboveLine(P,l,maxP), aboveLine(P,maxP,r));
  intT n1 = nn.first;
  intT n2 = nn.second;

  intT m1, m2;
  m1 = serialQuickHull(I,      P, n1, l,   maxP);
  m2 = serialQuickHull(I+n-n2, P, n2, maxP,r);
  for (intT i=0; i < m2; i++) I[i+m1+1] = I[i+n-n2];
  I[m1] = maxP;
  return m1+1+m2;
}

intT serialQuickHullP(point2d* P, intT n, point2d l, point2d r) {
  if (n < 2) return n;
  //point2d &maxP = P[0];//check invariance here
  point2d maxP = P[0];
  double maxArea = triArea(l,r,maxP);
  for (intT i=1; i < n; i++) {
    point2d j = P[i];
    double a = triArea(l,r,j);
    if (a > maxArea) {
      maxArea = a;
      maxP = j;
    }
  }

  pair<intT,intT> nn = split(P, n, aboveLineP(P,l,maxP), aboveLineP(P,maxP,r));
  intT n1 = nn.first;
  intT n2 = nn.second;

  intT m1, m2;
  m1 = serialQuickHullP(P, n1, l,   maxP);
  m2 = serialQuickHullP(P+n-n2, n2, maxP,r);
  for (intT i=0; i < m2; i++) P[i+m1+1] = P[i+n-n2];
  P[m1] = maxP;
  return m1+1+m2;
}


struct triangArea {
  intT l, r;
  point2d* P;
  intT* I;
  triangArea(intT* _I, point2d* _P, intT _l, intT _r) : I(_I), P(_P), l(_l), r(_r) {}
  double operator() (intT i) {return triArea(P[l], P[r], P[I[i]]);}
};

struct triangAreaP {
  point2d l, r;
  point2d* P;
  triangAreaP(point2d* _P, point2d &_l, point2d &_r) : P(_P), l(_l), r(_r) {}
  double operator() (intT i) {return triArea(l, r, P[i]);}
};

// intT quickHull(intT* I, intT* Itmp, point2d* P, intT n, intT l, intT r, intT depth, bool reverse) {
//   if(!reverse){
//     if (n < 2 || depth == 0) {
//       intT res = serialQuickHull(I, P, n, l, r);
//       memcpy(Itmp, I, n*sizeof(intT));
//       return res;
//     } else {
//       intT m1, m2;

//       intT idx = maxIndex<double>((intT)0,n,greater<double>(),triangArea(I,P,l,r));
//       intT maxP = I[idx];

//       intT n1 = filter(I, Itmp,    n, aboveLine(P, l, maxP));
//       intT n2 = filter(I, Itmp+n1+1, n, aboveLine(P, maxP, r));

//       m1 = cilk_spawn quickHull(I, Itmp ,P, n1, l, maxP, depth-1, !reverse);
//       m2 = quickHull(I+n1+1, Itmp+n1+1, P, n2, maxP, r, depth-1, !reverse);
//       cilk_sync;
//       Itmp[m1] = maxP;

//       return m1+1+m2;
//     }
//   }else{
//     if (n < 2 || depth == 0) {
//       intT res = serialQuickHull(Itmp, P, n, l, r);
//       memcpy(I, Itmp, n*sizeof(intT));
//       return res;
//     }else{
//       intT m1, m2;

//       intT idx = maxIndex<double>((intT)0,n,greater<double>(),triangArea(Itmp,P,l,r));
//       intT maxP = Itmp[idx];

//       intT n1 = filter(Itmp, I,    n, aboveLine(P, l, maxP));
//       intT n2 = filter(Itmp, I+n1+1, n, aboveLine(P, maxP, r));

//       m1 = cilk_spawn quickHull(I, Itmp ,P, n1, l, maxP, depth-1, !reverse);
//       m2 = quickHull(I+n1+1, Itmp+n1+1, P, n2, maxP, r, depth-1, !reverse);
//       cilk_sync;
//       I[m1] = maxP;
//       return m1+1+m2;
//     }
//   } 
//     //parallel_for (intT i=0; i < m1; i++) I[i] = Itmp[i];
//     //I[m1] = maxP;
//     //parallel_for (intT i=0; i < m2; i++) I[i+m1+1] = Itmp[i+n1];    
// }

// intT quickHull(intT* I, intT* Itmp, point2d* P, intT n, intT l, intT r, intT depth) {
//   if (n < 2 || depth == 0) 
//     return serialQuickHull(I, P, n, l, r);
//   else {
    
//     intT idx = maxIndex<double>((intT)0,n,greater<double>(),triangArea(I,P,l,r));
//     intT maxP = I[idx];

//     intT n1 = filter(I, Itmp,    n, aboveLine(P, l, maxP));
//     intT n2 = filter(I, Itmp+n1, n, aboveLine(P, maxP, r));

//     intT m1, m2;
//     m1 = cilk_spawn quickHull(Itmp, I ,P, n1, l, maxP, depth-1);
//     m2 = quickHull(Itmp+n1, I+n1, P, n2, maxP, r, depth-1);
//     cilk_sync;

//     //parallel_for (intT i=0; i < m1; i++) I[i] = Itmp[i];
//     //I[m1] = maxP;
//     //parallel_for (intT i=0; i < m2; i++) I[i+m1+1] = Itmp[i+n1];
    
//     return m1+1+m2;
//   }
// }

intT quickHull(intT* I, intT* Itmp, point2d* P, intT n, intT l, intT r, intT depth) {
  if (n < 2 || depth == 0) 
    return serialQuickHull(I, P, n, l, r);
  else {
    intT idx = maxIndex<double>((intT)0,n,greater<double>(),triangArea(I,P,l,r));
    intT maxP = I[idx];

    intT n1 = filter(I, Itmp,    n, aboveLine(P, l, maxP));
    intT n2 = filter(I, Itmp+n1, n, aboveLine(P, maxP, r));

    intT m1, m2;
    m1 = cilk_spawn quickHull(Itmp, I ,P, n1, l, maxP, depth-1);
    m2 = quickHull(Itmp+n1, I+n1, P, n2, maxP, r, depth-1);
    cilk_sync;

    parallel_for (intT i=0; i < m1; i++) I[i] = Itmp[i];
    I[m1] = maxP;
    parallel_for (intT i=0; i < m2; i++) I[i+m1+1] = Itmp[i+n1];
    
    return m1+1+m2;
  }
}
const int coersen_num = 10000;
intT quickHullP(point2d* P, point2d* Ptmp, intT n, point2d l, point2d r, intT depth) {
  if (n < coersen_num){// ||depth == 0){
    //cout << depth << endl;
      return serialQuickHullP(P, n, l, r);
  }else {
    if (depth == 5){
      //cout << "in. spawning thread is: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) << " at " << P<<endl;
    }

    intT idx = maxIndex<double>((intT)0,n,greater<double>(),triangAreaP(P,l,r));
    point2d maxP = P[idx];
    //cout << depth << endl;
    
    //intT maxP = idx;

    intT n1 = filter(P, Ptmp,    n, aboveLineP(P, l, maxP));
    intT n2 = filter(P, Ptmp+n1, n, aboveLineP(P, maxP, r));

    intT m1, m2;
    
    if (depth == 5){
      //cout << "before recursion. spawning thread is: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) << " at " << P<<endl;
    }

    m1 = cilk_spawn quickHullP(Ptmp, P, n1, l, maxP, depth-1);
    m2 = quickHullP(Ptmp+n1, P+n1, n2, maxP, r, depth-1);
    cilk_sync;
    parallel_for (intT i=0; i < m1; i++) P[i] = Ptmp[i];
    P[m1] = maxP;
    parallel_for (intT i=0; i < m2; i++) P[i+m1+1] = Ptmp[i+n1];

    if (depth == 5){
      //cout << "before return. spawning thread is: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) << " at " << P<<endl;
    }

    return m1+1+m2;
  }
}


struct makePair {
  pair<intT,intT> operator () (intT i) { return pair<intT,intT>(i,i);}
};

struct minMaxIndex {
  point2d* P;
  minMaxIndex (point2d* _P) : P(_P) {}
  pair<intT,intT> operator () (pair<intT,intT> l, pair<intT,intT> r) {
    intT minIndex = 
      (P[l.first].x < P[r.first].x) ? l.first :
      (P[l.first].x > P[r.first].x) ? r.first :
      (P[l.first].y < P[r.first].y) ? l.first : r.first;
    intT maxIndex = (P[l.second].x > P[r.second].x) ? l.second : r.second;
    return pair<intT,intT>(minIndex, maxIndex);
  }
};
bool compare_x(const point2d &p1, const point2d &p2){
  if(p1.x < p2.x){
    return true;
  }else if (p1.x == p2.x){
    return p1.y < p2.y;
  }else{
    return false;
  }
  //return p1.x < p2.x;
}
bool compare_y(const point2d &p1, const point2d &p2){
  if(p1.y < p2.y){
    return true;
  }else if(p1.x == p2.x){
    return p1.x < p2.x;
  }else{
    return false;
  }
}
static const intT minmax_base = 2000;
pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > merge(pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > minmax_e1, pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > minmax_e2){

  bool cmp_x_min = compare_x(*minmax_e1.first.first, *minmax_e2.first.first);
  bool cmp_x_max = !compare_x(*minmax_e1.first.second, *minmax_e2.first.second);
  bool cmp_y_min = compare_y(*minmax_e1.second.first, *minmax_e2.second.first);
  bool cmp_y_max = !compare_y(*minmax_e1.second.second, *minmax_e2.second.second);
  
  point2d *min_x = cmp_x_min ? minmax_e1.first.first : minmax_e2.first.first;
  
  point2d *max_x = cmp_x_max ? minmax_e1.first.second : minmax_e2.first.second;
  
  point2d *min_y = cmp_y_min ? minmax_e1.second.first : minmax_e2.second.first;
    
  point2d *max_y = cmp_y_max ? minmax_e1.second.second : minmax_e2.second.second;

  return make_pair(make_pair(min_x, max_x), make_pair(min_y, max_y));
}
pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > find_minmax_xy(point2d *p,  intT n){
  if(n < minmax_base){
    pair<point2d *, point2d *> minmax_ex = minmax_element(p, p+n, compare_x);
    pair<point2d *, point2d *> minmax_ey = minmax_element(p, p+n, compare_y);
    return make_pair(minmax_ex, minmax_ey);
  }else{
    pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > minmax_e1 = cilk_spawn find_minmax_xy(p, n/2);
    pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > minmax_e2 = find_minmax_xy(p+n/2, n - n/2);
    cilk_sync;
    return merge(minmax_e1, minmax_e2);
  }
}

// pair<pair<intT, intT>, pair<intT, intT> > find_minmax_xy(point2d *p, intT n){
//   return find_min_max_xy_r(p,);
// }

_seq<intT> hull(point2d* P, intT n) {
  pair<intT,intT> minMax = reduce<pair<intT,intT> >((intT)0,n,minMaxIndex(P), makePair());
  intT l = minMax.first;
  intT r = minMax.second;
  bool* fTop = newA(bool,n);
  bool* fBot = newA(bool,n);
  intT* I = newA(intT, n);
  intT* Itmp = newA(intT, n);
  parallel_for(intT i=0; i < n; i++) {
    Itmp[i] = i;
    double a = triArea(P[l],P[r],P[i]);
    fTop[i] = a > 0;
    fBot[i] = a < 0;
  }

  intT n1 = pack(Itmp, I, fTop, n);
  intT n2 = pack(Itmp, I+n1, fBot, n);
  free(fTop); free(fBot);

  intT m1; intT m2;
  m1 = cilk_spawn quickHull(I, Itmp, P, n1, l, r, 5);
  m2 = quickHull(I+n1, Itmp+n1, P, n2, r, l, 5);
  cilk_sync;

  parallel_for (intT i=0; i < m1; i++) Itmp[i+1] = I[i];
  parallel_for (intT i=0; i < m2; i++) Itmp[i+m1+2] = I[i+n1];
  free(I);
  Itmp[0] = l;
  Itmp[m1+1] = r;
  return _seq<intT>(Itmp, m1+2+m2);
}
point2d *offset_helper(point2d * buf, size_t off){
  char *tmp = (char *)buf;
  tmp = tmp + off;
  return (point2d *)tmp;
}

template <class ET, class intT, class PRED> 
intT wrapped_filter1(ET* In, ET* Out, intT n, PRED p) {
  intT num_pages = (n * sizeof(point2d));
  int n1 = filter(offset_helper(In, 0), Out, n/4, p);
  //cout << "wrapped11: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu())<<endl;
  int n2 = filter(offset_helper(In, num_pages), Out+n1, n/4, p);
  //cout << "wrapped12: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu())<<endl;
  int n3 = filter(offset_helper(In, num_pages*2), Out+n1+n2, n/4, p);
  //cout << "wrapped13: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu())<<endl;
  int n4 = filter(offset_helper(In, num_pages*3), Out+n1+n2+n3, n-n/4, p);
  //cout << "wrapped14: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu())<<endl;
  return n1+n2+n3+n4;
}

template <class ET, class intT, class PRED> 
intT wrapped_filter2(ET* In, ET* Out, intT n, PRED p) {
  intT num_pages = (n * sizeof(point2d));
  int n1 = filter(offset_helper(In, num_pages), Out, n/4, p);
  //cout << "wrapped21: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n2 = filter(offset_helper(In, 0), Out+n1, n/4, p);
  //cout << "wrapped22: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n3 = filter(offset_helper(In, num_pages*3), Out+n1+n2, n/4, p);
  //cout << "wrapped23: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n4 = filter(offset_helper(In, num_pages*2), Out+n1+n2+n3, n-n/4, p);
  //cout << "wrapped24: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  return n1+n2+n3+n4;
}

template <class ET, class intT, class PRED> 
intT wrapped_filter3(ET* In, ET* Out, intT n, PRED p) {
  intT num_pages = (n * sizeof(point2d));
  int n1 = filter(offset_helper(In, num_pages*2), Out, n/4, p);
  //cout << "wrapped31: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n2 = filter(offset_helper(In, num_pages*3), Out+n1, n/4, p);
  //cout << "wrapped32: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n3 = filter(offset_helper(In, num_pages*1), Out+n1+n2, n/4, p);
  //cout << "wrapped33: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n4 = filter(offset_helper(In, 0), Out+n1+n2+n3, n-n/4, p);
  //cout << "wrapped34: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  return n1+n2+n3+n4;
}
template <class ET, class intT, class PRED> 
intT wrapped_filter4(ET* In, ET* Out, intT n, PRED p) {
  intT num_pages = (n * sizeof(point2d));
  int n1 = filter(offset_helper(In, num_pages*3), Out, n/4, p);
  //cout << "wrapped41: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n2 = filter(offset_helper(In, num_pages*2), Out+n1, n/4, p);
  //cout << "wrapped42: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n3 = filter(offset_helper(In, 0), Out+n1+n2, n/4, p);
  //cout << "wrapped43: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n4 = filter(offset_helper(In, num_pages*1), Out+n1+n2+n3, n-n/4, p);
  //cout << "wrapped44: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  return n1+n2+n3+n4;
}

template <class ET, class intT, class PRED> 
intT wrapped_filter12(ET* In, ET* Out, intT n, PRED p) {
  intT num_pages = (n * sizeof(point2d));
  int n1 = cilk_spawn filter(offset_helper(In, 0), Out, n/4, p);
  //cout << "wrapped11: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu())<<endl;
  int n2 = cilk_spawn filter(offset_helper(In, num_pages), Out+n/4, n/4, p);
  //cout << "wrapped12: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu())<<endl;
  int n3 = cilk_spawn filter(offset_helper(In, num_pages*2), Out+n/4*2, n/4, p);
  //cout << "wrapped13: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu())<<endl;
  int n4 = cilk_spawn filter(offset_helper(In, num_pages*3), Out+n/4*3, n-n/4, p);
  //cout << "wrapped14: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu())<<endl;
  cilk_sync;
  return n1+n2+n3+n4;
}

template <class ET, class intT, class PRED> 
intT wrapped_filter22(ET* In, ET* Out, intT n, PRED p) {
  intT num_pages = (n * sizeof(point2d));
  int n1 = cilk_spawn filter(offset_helper(In, num_pages), Out, n/4, p);
  //cout << "wrapped21: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n2 = cilk_spawn filter(offset_helper(In, 0), Out+n/4, n/4, p);
  //cout << "wrapped22: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n3 = cilk_spawn filter(offset_helper(In, num_pages*3), Out+n/4*2, n/4, p);
  //cout << "wrapped23: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n4 = cilk_spawn filter(offset_helper(In, num_pages*2), Out+n/4*3, n-n/4, p);
  //cout << "wrapped24: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  cilk_sync;
  return n1+n2+n3+n4;
}

template <class ET, class intT, class PRED> 
intT wrapped_filter32(ET* In, ET* Out, intT n, PRED p) {
  intT num_pages = (n * sizeof(point2d));
  int n1 = cilk_spawn filter(offset_helper(In, num_pages*2), Out, n/4, p);
  //cout << "wrapped31: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n2 = cilk_spawn filter(offset_helper(In, num_pages*3), Out+n/4, n/4, p);
  //cout << "wrapped32: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n3 = cilk_spawn filter(offset_helper(In, num_pages*1), Out+n/4*2, n/4, p);
  //cout << "wrapped33: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n4 = cilk_spawn filter(offset_helper(In, 0), Out+n/4*3, n-n/4, p);
  //cout << "wrapped34: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  cilk_sync;
  return n1+n2+n3+n4;
}

template <class ET, class intT, class PRED> 
intT wrapped_filter42(ET* In, ET* Out, intT n, PRED p) {
  intT num_pages = (n * sizeof(point2d));
  int n1 = cilk_spawn filter(offset_helper(In, num_pages*3), Out, n/4, p);
  //cout << "wrapped41: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n2 = cilk_spawn filter(offset_helper(In, num_pages*2), Out+n/4, n/4, p);
  //cout << "wrapped42: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n3 = cilk_spawn filter(offset_helper(In, 0), Out+n/4*2, n/4, p);
  //cout << "wrapped43: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  int n4 = cilk_spawn filter(offset_helper(In, num_pages*1), Out+n/4*3, n-n/4, p);
  //cout << "wrapped44: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  cilk_sync;
  return n1+n2+n3+n4;
}
void block_for_func(intT *Sums, bool *Fl, intT s, intT e){
  
}

template <class ET, class intT, class F> 
_seq<ET> packSerialm(ET* Out, bool* Fl, intT s, intT e, F f) {
  if (Out == NULL) {
    intT m = sumFlagsSerial(Fl+s, e-s);
    Out = newA(ET,m);
  }
  intT k = 0;
  for (intT i=s; i < e; i++) {
    if (Fl[i]) {
      Out[k++] = f(i);
    }
  }
  return _seq<ET>(Out,k);
}

template <class ET, class intT, class F> 
void out_for_func(ET *Out, bool *Fl, intT s, intT e, intT *Sums, F f){
  intT n = e-s;
  //intT l = el - sl;
  intT l = nblocks(n, _F_BSIZE);
  parallel_for (intT i = 0; i < l; i++) {
    intT cs = (i) * (_F_BSIZE);
    intT ce;
    ce = min(cs + (_F_BSIZE), e - s);

    packSerialm(Out+Sums[i], Fl+s, cs, ce, f);
  }						
}
template <class intT> 
void sum_for_func(intT *Sums, bool *Fl, intT s, intT e){
  blocked_for (i, s, e, _F_BSIZE, Sums[i] = sumFlagsSerial(Fl+s, e-s););  
}
template <class ET, class intT> 
_seq<ET> packm(ET *In, ET* Out, bool* Fl, intT n) {
  //intT l = nblocks(n, _F_BSIZE);
  //if (l <= 1) return packSerial(Out, Fl, s, e, f);
  intT l14 = nblocks(n/4, _F_BSIZE);
  intT ll = nblocks(n-n/4*3, _F_BSIZE);

  intT *Sums = newA(intT,l14*3+ll);
  intT s = 0;
  intT e = n;
  cilk_spawn sum_for_func(Sums, Fl, (intT)0, n/4);
  cilk_spawn sum_for_func(Sums+l14, Fl, n/4, n/4*2);
  cilk_spawn sum_for_func(Sums+l14*2, Fl, n/4*2, n/4*3);
             sum_for_func(Sums+l14*3, Fl, n/4*3, n);
  cilk_sync;

  intT m = plusScan(Sums, Sums, l14*3+ll);
  intT num_pages = (n * sizeof(ET));

  cilk_spawn out_for_func(Out, Fl, (intT)0,     n/4,   Sums, getA<ET,intT>(offset_helper(In, num_pages*0)));
  cilk_spawn out_for_func(Out, Fl, (intT)n/4,   n/4*2, Sums+l14, getA<ET,intT>(offset_helper(In, num_pages*1)));
  cilk_spawn out_for_func(Out, Fl, (intT)n/4*2, n/4*3, Sums+l14*2, getA<ET,intT>(offset_helper(In, num_pages*2)));
  out_for_func(Out,            Fl, (intT)n/4*3, n,     Sums+l14*3, getA<ET,intT>(offset_helper(In, num_pages*3)));
  cilk_sync;
  free(Sums);
  return _seq<ET>(Out,m);
}

template <class ET, class intT> 
intT packmw(ET* In, ET* Out, bool* Fl, intT n) {
  // intT num_pages = (n * sizeof(ET));
  // int m1= cilk_spawn pack(Out, Fl, (intT) 0, n/4, getA<ET,intT>(offset_helper(In, 0))).n;
  // int m2= cilk_spawn pack(Out, Fl, n/4, n/4*2, getA<ET,intT>(offset_helper(In, num_pages))).n;
  // int m3= cilk_spawn pack(Out, Fl, n/4*2, n/4*3, getA<ET,intT>(offset_helper(In, num_pages*2))).n;
  // int m4= pack(Out, Fl, n/4*3, n-n/4*3, getA<ET,intT>(offset_helper(In, num_pages*3))).n;
  // cilk_sync;
  // return m1+m2+m3+m4;
  return packm(In, Out, Fl, n).n;
}

template<class ET, class intT, class PRED>
void parallel_for_func(ET *In, bool*Fl, intT n, PRED p){
    parallel_for (intT i=0; i < n; i++) Fl[i] = (bool) p(In[i]);
}
template <class ET, class intT, class PRED> 
intT wrapped_filter_new(ET* In, ET* Out, intT n, PRED p) {
  intT num_pages = (n * sizeof(ET));
  bool *Fl = newA(bool,n);
  cilk_spawn parallel_for_func(offset_helper(In, 0), Fl, n/4, p);
  cilk_spawn parallel_for_func(offset_helper(In, num_pages), Fl+n/4, n/4, p);
  cilk_spawn parallel_for_func(offset_helper(In, num_pages*2), Fl+n/4*2, n/4, p);
  parallel_for_func(offset_helper(In, num_pages*3), Fl+n/4*3, n-n/4*3, p);
  cilk_sync;
  
  intT  m = packmw(In, Out, Fl, n);
  free(Fl);
  return m;
}

_seq<point2d> hullP(point2d* P, intT n, point2d *Ptmp) {
  // cout << "n is " << n << endl;
  //cout << "point size is " << sizeof(point2d) << endl;
  intT num_pages = (n * sizeof(point2d));
  #ifndef NO_PIN
  intT num_sockets = __cilkrts_num_sockets();
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
  #endif
  clockmark_t start = ktiming_getmark();
  startTime();
  //printf("pin patttern: %d, %d, %d\n", pin_pattern[0], pin_pattern[1], pin_pattern[2]);
  __cilkrts_disable_nonlocal_steal();
  SET_PIN(pin_pattern[0]);
  pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > minmaxxy1 = cilk_spawn find_minmax_xy(P, n/4);

  SET_PIN(pin_pattern[1]);
  pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > minmaxxy2 = cilk_spawn find_minmax_xy(offset_helper(P, num_pages), n/4);

  SET_PIN(pin_pattern[2]);
  pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > minmaxxy3 = cilk_spawn find_minmax_xy(offset_helper(P, num_pages*2), n/4);

  __cilkrts_unset_pinning_info();
  __cilkrts_enable_nonlocal_steal();
  pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > minmaxxy4 = find_minmax_xy(offset_helper(P, num_pages*3), n - n/4*3);

  __cilkrts_set_pinning_info(0);
  cilk_sync;
  pair<pair<point2d *, point2d *>, pair<point2d *, point2d *> > minmaxxy = merge(minmaxxy1, minmaxxy2);
  minmaxxy = merge(minmaxxy, minmaxxy3);
  minmaxxy = merge(minmaxxy, minmaxxy4);
  
  //nextTime("findmin");
  
  point2d l = *minmaxxy.first.first;
  point2d r = *minmaxxy.first.second;

  point2d b = *minmaxxy.second.first;
  point2d t = *minmaxxy.second.second;


  //cout << "l " << l << " r " << r << " b " << b << " t " << t << endl;
  // __cilkrts_disable_nonlocal_steal();
  // __cilkrts_set_pinning_info(1);
  // intT n11 = cilk_spawn filter(offset_helper(P, 0), offset_helper(Ptmp, 0),  n/4, aboveLineP(P, l, t));
  // __cilkrts_set_pinning_info(2);
  // intT n21 = cilk_spawn filter(offset_helper(P, num_pages), offset_helper(Ptmp, num_pages), n/4, aboveLineP(P, t, r));
  // __cilkrts_set_pinning_info(3);
  // intT n31 = cilk_spawn filter(offset_helper(P, num_pages*2), offset_helper(Ptmp, num_pages*2), n/4, aboveLineP(P, r, b));
  // __cilkrts_unset_pinning_info();
  // __cilkrts_enable_nonlocal_steal();
  // intT n41 = filter(offset_helper(P, num_pages*3), offset_helper(Ptmp, num_pages*3), n-n/4*3, aboveLineP(P, b, l));
  // cilk_sync;

  // __cilkrts_disable_nonlocal_steal();
  // __cilkrts_set_pinning_info(1);
  // intT n12 = cilk_spawn filter(offset_helper(P, num_pages), offset_helper(Ptmp, 0)+n11,  n/4, aboveLineP(P, l, t));
  // __cilkrts_set_pinning_info(2);
  // intT n22 = cilk_spawn filter(offset_helper(P, 0), offset_helper(Ptmp, num_pages)+n21, n/4, aboveLineP(P, t, r));
  // __cilkrts_set_pinning_info(3);
  // intT n32 = cilk_spawn filter(offset_helper(P, num_pages*3), offset_helper(Ptmp, num_pages*2)+n31, n/4, aboveLineP(P, r, b));
  // __cilkrts_unset_pinning_info();
  // __cilkrts_enable_nonlocal_steal();
  // intT n42 = filter(offset_helper(P, num_pages*2), offset_helper(Ptmp, num_pages*3)+n41, n-n/4*3, aboveLineP(P, b, l));
  // cilk_sync;

  // __cilkrts_disable_nonlocal_steal();
  // __cilkrts_set_pinning_info(1);
  // intT n13 = cilk_spawn filter(offset_helper(P, num_pages*2), offset_helper(Ptmp, 0)+n11+n12,  n/4, aboveLineP(P, l, t));
  // __cilkrts_set_pinning_info(2);
  // intT n23 = cilk_spawn filter(offset_helper(P, num_pages*3), offset_helper(Ptmp, num_pages)+n21+n22, n/4, aboveLineP(P, t, r));
  // __cilkrts_set_pinning_info(3);
  // intT n33 = cilk_spawn filter(offset_helper(P, num_pages), offset_helper(Ptmp, num_pages*2)+n31+n32, n/4, aboveLineP(P, r, b));
  // __cilkrts_unset_pinning_info();
  // __cilkrts_enable_nonlocal_steal();
  // intT n43 = filter(offset_helper(P, 0), offset_helper(Ptmp, num_pages*3)+n41+n42, n-n/4*3, aboveLineP(P, b, l));
  // cilk_sync;

  // __cilkrts_disable_nonlocal_steal();
  // __cilkrts_set_pinning_info(1);
  // intT n14 = cilk_spawn filter(offset_helper(P, num_pages*3), offset_helper(Ptmp, 0)+n11+n12+n13,  n/4, aboveLineP(P, l, t));
  // __cilkrts_set_pinning_info(2);
  // intT n24 = cilk_spawn filter(offset_helper(P, num_pages*2), offset_helper(Ptmp, num_pages)+n21+n22+n23, n/4, aboveLineP(P, t, r));
  // __cilkrts_set_pinning_info(3);
  // intT n34 = cilk_spawn filter(offset_helper(P, 0), offset_helper(Ptmp, num_pages*2)+n31+n32+n32, n/4, aboveLineP(P, r, b));
  // __cilkrts_unset_pinning_info();
  // __cilkrts_enable_nonlocal_steal();
  // intT n44 = filter(offset_helper(P, num_pages), offset_helper(Ptmp, num_pages*3)+n41+n42+n43, n-n/4*3, aboveLineP(P, b, l));
  // cilk_sync;

  // int n1 = n11+ n12 + n13 + n14;
  // int n2 = n21+ n22 + n23 + n24;
  // int n3 = n31+ n32 + n33 + n34;
  // int n4 = n41+ n42 + n43 + n44;

  /// function based!!!!!!

  __cilkrts_disable_nonlocal_steal();
  SET_PIN(pin_pattern[0]);
  intT n1 = cilk_spawn wrapped_filter_new(P, offset_helper(Ptmp, 0),  n, aboveLineP(P, l, t));
  SET_PIN(pin_pattern[1]);
  intT n2 = cilk_spawn wrapped_filter_new(P, offset_helper(Ptmp, num_pages), n, aboveLineP(P, t, r));
  SET_PIN(pin_pattern[2]);
  intT n3 = cilk_spawn wrapped_filter_new(P, offset_helper(Ptmp, num_pages*2), n, aboveLineP(P, r, b));
  __cilkrts_unset_pinning_info();
  __cilkrts_enable_nonlocal_steal(); 
  intT n4 = wrapped_filter_new(P, offset_helper(Ptmp, num_pages*3), n, aboveLineP(P, b, l));
  __cilkrts_set_pinning_info(0);
  cilk_sync;
  //function based !!!!!

  //nextTime("filter");

  //cout << "n1 " << n1 << " n2 " << n2 << " n3 "<< n3 << " n4 " << n4 << endl;

  // for(int i = 0; i < n1; i++){
  //   cout << Ptmp[i] << endl;
  // }
  // intT m1; intT m2; 
  intT m1; intT m2; intT m3; intT m4;

  //cout << "top level: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) <<endl;
  //usable here
  __cilkrts_disable_nonlocal_steal();
  SET_PIN(pin_pattern[0]);

  //print_mem_binding(P, 1);
  //print_mem_binding(Ptmp, 1);
  //cout << "1. spawning thread is: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) << " at " << Ptmp<<endl;
  m1 = cilk_spawn quickHullP(Ptmp, P, n1, l, t, 5);
  SET_PIN(pin_pattern[1]);
  //print_mem_binding(offset_helper(P, num_pages), 1);
  //print_mem_binding(offset_helper(Ptmp, num_pages), 1);
  //cout << "2. spawning thread is: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) << endl;
  m2 = cilk_spawn quickHullP(offset_helper(Ptmp, num_pages), offset_helper(P, num_pages), n2, t, r, 5);
  SET_PIN(pin_pattern[2]);
  
  //print_mem_binding(offset_helper(P, num_pages*2), 1);
  //print_mem_binding(offset_helper(Ptmp, num_pages*2), 1);
  //cout << "3. spawning thread is: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) << endl;
  m3 = cilk_spawn quickHullP(offset_helper(Ptmp, num_pages*2), offset_helper(P, num_pages*2), n3, r, b, 5);
  __cilkrts_unset_pinning_info();
  __cilkrts_enable_nonlocal_steal();//enable non local is here!!!!!!!!  
  //print_mem_binding(offset_helper(P, num_pages*3), 1);
  //print_mem_binding(offset_helper(Ptmp, num_pages*3), 1);
  //cout << "4. spawning thread is: " << sched_getcpu() << "which is on: " << numa_node_of_cpu(sched_getcpu()) << endl;
  m4 = quickHullP(offset_helper(Ptmp, num_pages*3), offset_helper(P, num_pages*3), n4, b, l, 5);
  __cilkrts_set_pinning_info(0);
  cilk_sync;
  //nextTime("spawn");

  //cout << "m1 " << m1 << " m2 " << m2 << " m3 " << m3 << " m4 " << m4 << endl; 
  int offset = 0;
  if (l.x != t.x || l.y != t.y){
    offset++;
  }
  parallel_for (intT i=0; i < m1; i++) P[i+offset] = Ptmp[i];
  if (t.x != r.x || t.y != r.y){
    offset++;
  }
  parallel_for (intT i=0; i < m2; i++) P[i+m1+offset] = offset_helper(Ptmp, num_pages)[i];
  if (r.x != b.x || r.y != b.y){
    offset++;
  }
  parallel_for (intT i=0; i < m3; i++) P[i+m1+m2+offset] = offset_helper(Ptmp, num_pages*2)[i];
  if (b.x != l.x || b.y != l.y){
    offset++;
  }
  parallel_for (intT i=0; i < m4; i++) P[i+m1+m2+m3+offset] = offset_helper(Ptmp, num_pages*3)[i];

  //nextTime("copy back");
  //clockmark_t end = ktiming_getmark();
  //uint64_t time = ktiming_diff_usec(&start, &end);
  //print_runtime(&time, 1);
  //cout << "finish" << endl;

  //free(Ptmp); //WARNING: Ptmp is still allocated!!!!!
  int offset2 = 0;

  P[0] = l;
  offset2 += m1;
  if(l.x != t.x || l.y != t.y){
    offset2++;
    P[offset2] = t;
  }
  
  offset2 += m2;
  if(t.x != r.x || t.y != r.y){
    offset2++;
    P[offset2] = r;
  }

  offset2 += m3;
  if(t.x != r.x || t.y != r.y){
    offset2++;
    P[offset2] = b;
  }
  return _seq<point2d>(P, offset2+1);
}
