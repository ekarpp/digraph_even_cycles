SHELL := /bin/bash -O extglob

CXX := g++
CXXFLAGS := -g -std=c++1z -O3 -Wall -Wextra -march=native -fopenmp
LDFLAGS := -fopenmp

VPATH = src:tests/unit:tests/perf

BIN := digraph digraph-tests extension-perf gf-perf mem-bench

BASE_OBJ := gf.o extension.o fmatrix.o ematrix.o polynomial.o util.o solver.o graph.o
TEST_OBJ := gf_test.o extension_test.o fmatrix_test.o util_test.o solver_test.o ematrix_test.o geng_test.o
PERF_OBJ := extension.o polynomial.o gf.o util.o


###########
# RECIPES #
###########

.PHONY: clean all help
all: $(BIN) nauty/geng nauty/directg nauty/listg

clean:
	rm -f *.o *.s *.asm1 *.asm $(BIN)
	cd nauty && git clean -xf && git checkout .

help:
	@echo 'list of available commands:'
	@echo ''
	@echo '  solver:'
	@echo '    make digraph'
	@echo ''
	@echo '  unit testing:'
	@echo '    make digraph-test'
	@echo ''
	@echo '  runs predefined tests:'
	@echo '    make test'
	@echo ''
	@echo '  tests GF(2^$${b}) solver on all digraphs of $${V} vertices (up to isomorphism):'
	@echo '    EXP=$${b} VERT=$${V} make geng-test'
	@echo ''
	@echo '  GF performance benchamrking:'
	@echo '    make gf-perf'
	@echo ''
	@echo '  extension performance benchmarking:'
	@echo '    make extension-perf'
	@echo ''
	@echo '  memory bandwidth benchmarking:'
	@echo '    make mem-bench'
	@echo ''
	@echo '  clean-up artifacts and binaries:'
	@echo '    make clean'
	@echo ''


##########
# SOLVER #
##########

digraph: main.o $(BASE_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

#########
# TESTS #
#########

digraph-tests: tests.o $(TEST_OBJ) $(BASE_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

###########
# GF PERF #
###########

gf-perf: gf_perf.o $(PERF_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

############
# EXT PERF #
############

extension-perf: extension_perf.o $(PERF_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

#############
# MEM BENCH #
#############

mem-bench: mem_bench.o
	$(CXX) $^ $(LDFLAGS) -o $@

#########
# NAUTY #
#########

nauty/geng:
	cd nauty && [ -f config.log ] || ./configure && make geng

nauty/listg:
	cd nauty && [ -f config.log ] || ./configure && make listg

nauty/directg:
	cd nauty && [ -f config.log ] || ./configure && make directg


#################
# TEST COMMANDS #
#################

test32: digraph-tests
	./digraph-tests -n32 -efux -d20 -t10000
	./digraph-tests -n32 -s -d10 -t100

test16: digraph-tests
	./digraph-tests -n16 -egfux -d20 -t10000
	./digraph-tests -n16 -s -d10 -t100

test24: digraph-tests
	./digraph-tests -n24 -egfux -d20 -n20 -t10000
	./digraph-tests -n24 -s -d10 -n20 -t100

test: test16 test24 test32

geng-test: digraph-tests nauty/geng nauty/directg nauty/listg
	mkdir -p geng-fail/$(VERT)
	nauty/geng -q $(VERT) | nauty/directg -q | nauty/listg -aq | ./digraph-tests -n $(EXP) -c


#############
# ASM STUFF #
#############

%.s: %.cc
	$(CXX) -S $(CXXFLAGS) -fverbose-asm $^

%.asm1: %.s
	c++filt < $^ > $@

%.asm: %.o
	objdump -d -S $^ > $@


#################
# OBJECT RECIPE #
#################

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $^
