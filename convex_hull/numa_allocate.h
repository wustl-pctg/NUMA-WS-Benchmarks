/* Numa Allocate  */
#ifndef NUMA_ALLOCATE
#define NUMA_ALLOCATE

#define BITS_PER_PAGE 12
#define NUM_PGS_TO_BYTES(num_pgs) (num_pgs << BITS_PER_PAGE)
//#define NUM_PGS_TO_BYTES(num_pgs) (num_pgs * 4096)
//#define PAGE_OFFSET_MASK (~((0x1UL << BITS_PER_PAGE) - 1))
#define PAGE_OFFSET_MASK (~(0xffffffffffffffffL << BITS_PER_PAGE))

int get_node(void *p);

void free_pages(void *p, unsigned long long num_pgs);

char * allocate_pages(unsigned long long num_pgs);

void print_mem_binding(void *addr, int num_pgs);

int print_node(void *addr);

void print_mem_policy_on_thread();

void bind_memory(char *addr, unsigned long long *partition, int num_sockets);

char * bind_memory_numa(unsigned long long num_pgs, int num_sockets);

void pattern_bind_memory(char *addr, int num_blocks, int *pattern_array, int *partition);

char * pattern_bind_memory_numa(int num_pgs, int num_blocks, int *pattern_array);
#endif
