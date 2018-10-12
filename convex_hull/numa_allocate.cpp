/* Numa Allocate */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <numa.h>
#include <numaif.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <sys/mman.h>

#include "numa_allocate.h"

int get_node(void *p)
{
  int status[1];
  void *pa;
  unsigned long a;

  // round p down to the nearest page boundary

  a  = (unsigned long) p;
  a  = a - (a % ((unsigned long) getpagesize()));
  pa = (void *) a;

  if (move_pages(0,1,&pa,NULL,status,0) != 0) {
    abort();
  }

  return(status[0]);

}

/**
 * allocate num_pgs * pages (4K) of memory
 * returns pointer to the beginning of allocated memory
 * or NULL if allocation fails.
 **/
char * allocate_pages(unsigned long long num_pgs) {
  void *addr = NULL;
  addr = mmap(NULL, NUM_PGS_TO_BYTES(num_pgs), PROT_READ | PROT_WRITE,
	      MAP_PRIVATE | MAP_ANONYMOUS, -1/*fd*/, 0/*offset*/);
  assert(addr != NULL && addr != (void *)-1);
  return (char *)addr;
}

void free_pages(void *addr, unsigned long long num_pgs) {
  munmap(addr, NUM_PGS_TO_BYTES(num_pgs));
}

/**
 * Given addr and num_pgs (which specifies the memory range),
 * print where each page within the memory range is located.
 * NOTE: can also use move_pages to get the same info if NULL is passed for
 * node to move to.
 **/ 
void print_mem_binding(void *addr_in, int num_pgs){
  char *addr = (char *) addr_in;
  int node[num_pgs];
  char *curr_addr = addr;
  for(int i=0; i < num_pgs; i++) {
    node[i] = -1;
    int ret = get_mempolicy(&node[i], NULL, 0, curr_addr,
			    MPOL_F_NODE | MPOL_F_ADDR);
    assert(ret == 0);
    curr_addr = curr_addr + NUM_PGS_TO_BYTES(1);
    printf("Page %d allocated on socket %d.\n", i, node[i]);
  }
}

int print_node(void *addr) {

  int node;
  int ret = get_mempolicy(&node, NULL, 0, addr, MPOL_F_NODE | MPOL_F_ADDR);
  assert(ret == 0);

  return node;
}

void print_mem_policy_on_thread() {

  int possible_nodes = numa_num_possible_nodes();
  unsigned long nodemask[possible_nodes/64];
  errno = 0;
  int ret = get_mempolicy(NULL, nodemask, possible_nodes,
			  NULL, MPOL_F_MEMS_ALLOWED);
  assert(ret == 0);
  // we really only care about the first element, sinc we have e
  // only 4 nodes, but get_mempolicy forces us to use max possible nodes
  fprintf(stderr, "nmask = %lu.\n", nodemask[0]);
  if(errno) perror("Fail to get_mempolicy");
}

/**
 * addr: beginning of memory to be pinned
 * partition: an array of size num_sockets; partition[i] specifies the
 *            number of pages to pin on socket i.
 * Given beginning addr and partition, bind the memory according
 * to the partition, in the order of low to high addresses.
 **/
void bind_memory(char *addr, unsigned long long *partition, int num_sockets) {

  char *curr_addr = addr;
  // to get around an off-by-one bug
  // see https://www.spinics.net/lists/linux-mm/msg07089.html
  int num_nodes = numa_max_node() + 2;
  assert(num_sockets <= num_nodes);

  for(int i=0; i < num_sockets; i++) {
    unsigned long bytes = NUM_PGS_TO_BYTES(partition[i]);
    unsigned long nmask = 0x0;
    nmask = nmask | (0x1 << i);
    // fprintf(stderr, "binding addr %p of %lu bytes to %lx.\n",
    //        curr_addr, bytes, nmask);
    errno = 0;
    int ret = mbind(curr_addr, bytes, MPOL_BIND, &nmask,
		    num_nodes, 0/*flags*/);
    if(ret != 0) {
      if(errno) {
	perror("Fail to mbind");
	exit(1);
      }
    }
    // write to the first type to get the actual page allocation
    for(char *tmp = curr_addr;
	tmp < (curr_addr+bytes); tmp += NUM_PGS_TO_BYTES(1)) {
      tmp[0] = '0';
    }
    curr_addr += bytes;
  }
}

char * bind_memory_numa(unsigned long long num_pgs, int num_sockets) {

  assert(num_pgs > 0);
  unsigned long long partition[num_sockets];
  unsigned long long smallest_part = num_pgs / num_sockets;
  unsigned long long left_over = num_pgs - (smallest_part * num_sockets);
  unsigned long long sum = 0;
  // allocate a big chunck of memory
  char *addr = allocate_pages(num_pgs);
  // check the addr returned is page aligned
  assert(addr != NULL && ((unsigned long)addr & PAGE_OFFSET_MASK) == 0);

  // print_mem_policy_on_thread();

  // calculate how we want to partition the memory across sockets
  for(int i=num_sockets-1; i >= 0; i--) {
    if(left_over > 0) {
      partition[i] = smallest_part + 1;
      left_over = left_over - 1;
    } else {
      partition[i] = smallest_part;
    }
    sum += partition[i];
  }
  assert(sum == num_pgs);

  bind_memory(addr, partition, num_sockets); // perform memory binding
  //print_mem_binding(addr, num_pgs); // check memory binding

  return addr;
}

/**
 * addr: beginning of memory to be pinned
 * partition: an array of size num_sockets; partition[i] specifies the
 *            number of pages to pin on socket i.
 * Given beginning addr and partition, bind the memory according
 * to the partition, in the order of low to high addresses.
 **/
void pattern_bind_memory(char *addr, int num_blocks, int *pattern_array, int *partition) {

  char *curr_addr = addr;
  // to get around an off-by-one bug
  // see https://www.spinics.net/lists/linux-mm/msg07089.html
  int num_nodes = numa_max_node() + 2;

  for(int i=0; i < num_blocks; i++) {
    unsigned long bytes = NUM_PGS_TO_BYTES(partition[i]);
    unsigned long nmask = 0x0;
    nmask = nmask | (0x1 << pattern_array[i]);
    // fprintf(stderr, "binding addr %p of %lu bytes to %lx.\n",
    //        curr_addr, bytes, nmask);
    if(pattern_array[i] != -1) {
        errno = 0;
        int ret = mbind(curr_addr, bytes, MPOL_BIND, &nmask,
                num_nodes, 0/*flags*/);
        if(ret != 0) {
          if(errno) perror("Fail to mbind");
        }
    }
    // write to the first type to get the actual page allocation
    for(char *tmp = curr_addr;
	tmp < (curr_addr+bytes); tmp += NUM_PGS_TO_BYTES(1)) {
      tmp[0] = '0';
    }
    curr_addr += bytes;
  }
}

char * pattern_bind_memory_numa(int num_pgs, int num_blocks, int *pattern_array) {

  assert(num_pgs > 0);
  int partition[num_blocks];
  int smallest_part = num_pgs / num_blocks;
  int left_over = num_pgs - (smallest_part * num_blocks);
  int sum = 0;

  // allocate a big chunck of memory
  char *addr = allocate_pages(num_pgs);
  // check the addr returned is page aligned
  assert(addr != NULL && ((unsigned long)addr & PAGE_OFFSET_MASK) == 0);

  // print_mem_policy_on_thread();

  // calculate how we want to partition the memory across sockets
  for(int i=num_blocks-1; i >= 0; i--) {
    if(left_over > 0) {
      partition[i] = smallest_part + 1;
      left_over = left_over - 1;
    } else {
      partition[i] = smallest_part;
    }
    sum += partition[i];
  }
  assert(sum == num_pgs);

  pattern_bind_memory(addr, num_blocks, pattern_array, partition); // perform memory binding

  return addr;
}
