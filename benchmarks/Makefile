TARGETS := cilksort heat matmul matmul_z strassen strassen_z

default: $(TARGETS)

include ../common.mk

cilksort: getoptions.o cilksort.o ktiming.o numa_allocate.o
heat: getoptions.o heat.o ktiming.o numa_allocate.o
matmul: getoptions.o matmul_rm.o ktiming.o numa_allocate.o
matmul_z: getoptions.o matmul_z.o ktiming.o numa_allocate.o
strassen_z: getoptions.o strassen_z.o ktiming.o numa_allocate.o
strassen: getoptions.o strassen.o ktiming.o numa_allocate.o

%.s: %.cpp
	$(CXX) -S -fcilkplus -I$(INCLUDE) $^ -o $@

%: %.o
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@

move::
	mkdir -p ../bin
	mv $(TARGETS) ../bin
clean::
	rm -f $(TARGETS) *.o *.s
	cd ../bin; rm $(TARGETS)
