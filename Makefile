SHELL := /bin/bash -O extglob

CXX := g++
CXXFLAGS := -g -std=c++1z -O3 -Wall -Wextra -march=native -fopenmp
LDFLAGS := -fopenmp

VPATH = src:tests/unit:tests/perf

BIN := digraph digraph-tests extension-perf gf-perf mem-bench

OBJ := graph util gf extension fmatrix ematrix polynomial solver
OBJ := $(addsuffix $(bits).o, $(OBJ))

PERF_OBJ := extension polynomial gf util
PERF_OBJ := $(addsuffix $(bits).o, $(PERF_OBJ))

TEST_OBJ := gf_test extension_test fmatrix_test util_test solver_test ematrix_test geng_test
TEST_OBJ := $(addsuffix $(bits).o, $(TEST_OBJ))

ALL_OBJ := $(OBJ) $(PERF_OBJ) $(TEST_OBJ) solver$(bits)PAR.o extension_perf$(bits)PAR.o extension_perf$(bits).o main$(bits).o gf_perf$(bits).o tests$(bits).o

all: $(BIN) nauty/geng nauty/directg nauty/listg

CLEAN_REGEX := {0,16,32}?(-PAR)

.PHONY: clean clean-bin clean-obj
clean-obj:
	rm -f *.o
clean-bin:
	rm -f $(addsuffix $(CLEAN_REGEX), $(BIN)) mem-bench
clean: clean-obj clean-bin
	rm -f *.o *.s *.asm1 *.asm2
	cd nauty && git clean -xf && git checkout .

objectsX: $(ALL_OBJ)

objects:
	$(MAKE) objectsX bits=32
	$(MAKE) objectsX bits=16
	$(MAKE) objectsX bits=0

###################
# SOLVER BINARIES #
###################

digraph32:
	$(MAKE) digraphX bits=32

digraph16:
	$(MAKE) digraphX bits=16
	$(MAKE) digraph16-NOVEC bits=16

digraph16-NOVEC: main16.o $(OBJ) fmatrix16NOVEC.o
	$(CXX) main16.o $(subst fmatrix16.o, fmatrix16NOVEC.o, $(OBJ)) -o $@ $(LDFLAGS)

digraph0:
	$(MAKE) digraphX bits=0

digraphX: main$(bits).o $(OBJ) solver$(bits)PAR.o
	$(CXX) main$(bits).o $(OBJ) -o $@ $(LDFLAGS)
	mv $@ digraph$(bits)

	$(CXX) main$(bits).o $(subst solver$(bits).o, solver$(bits)PAR.o, $(OBJ)) -o $@ $(LDFLAGS)
	mv $@ digraph$(bits)-PAR

digraph: digraph32 digraph16 digraph0


#################
# TEST BINARIES #
#################

digraph-tests32:
	$(MAKE) digraph-testsX bits=32

digraph-tests16:
	$(MAKE) digraph-testsX bits=16

digraph-tests0:
	$(MAKE) digraph-testsX bits=0

digraph-testsX: tests$(bits).o $(OBJ) $(TEST_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)
	mv $@ digraph-tests$(bits)

digraph-tests: digraph-tests32 digraph-tests16 digraph-tests0


####################
# GF PERF BINARIES #
####################

gf-perf32:
	$(MAKE) gf-perfX bits=32

gf-perf16:
	$(MAKE) gf-perfX bits=16

gf-perf0:
	$(MAKE) gf-perfX bits=0

gf-perfX: gf_perf$(bits).o $(PERF_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)
	mv $@ gf-perf$(bits)

gf-perf: gf-perf32 gf-perf16 gf-perf0


#####################
# EXT PERF BINARIES #
#####################

extension-perf32:
	$(MAKE) extension-perfX bits=32

extension-perf16:
	$(MAKE) extension-perfX bits=16

extension-perf0:
	$(MAKE) extension-perfX bits=0

extension-perfX: extension_perf$(bits).o $(PERF_OBJ) extension_perf$(bits)PAR.o
	$(CXX) extension_perf$(bits).o $(PERF_OBJ) -o $@ $(LDFLAGS)
	mv $@ extension-perf$(bits)

	$(CXX) extension_perf$(bits)PAR.o $(PERF_OBJ) -o $@ $(LDFLAGS)
	mv $@ extension-perf$(bits)-PAR

extension-perf: extension-perf32 extension-perf16 extension-perf0

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

test32: digraph-tests32
	./digraph-tests32 -efux -d20 -t10000
	./digraph-tests32 -s -d10 -t100

test16: digraph-tests16
	./digraph-tests16 -egfux -d20 -t10000
	./digraph-tests16 -s -d10 -t100

test0: digraph-tests0
	./digraph-tests0 -egfux -d20 -n20 -t10000
	./digraph-tests0 -s -d10 -n20 -t100

test: test32 test16 test0

geng-test: digraph-tests nauty/geng nauty/directg nauty/listg
	mkdir -p geng-fail/$(vert)
	nauty/geng -q $(vert) | nauty/directg -q | nauty/listg -aq | ./digraph-tests$(BITS) -c


#############
# ASM STUFF #
#############

%.s: %.cc
	$(CXX) -D GF2_bits=16 -S $(CXXFLAGS) -fverbose-asm $^

%.asm1: %.s
	c++filt < $^ > $@

%.asm2: %16.o
	objdump -d -S $^ > $@


##################
# OBJECT RECIPES #
##################

%32.o: %.cc
	$(CXX) $(CXXFLAGS) -D GF2_bits=32 -c -o $@ $^

%16.o: %.cc
	$(CXX) $(CXXFLAGS) -D GF2_bits=16 -c -o $@ $^

%0.o: %.cc
	$(CXX) $(CXXFLAGS) -D GF2_bits=0 -c -o $@ $^

%32PAR.o: %.cc
	$(CXX) $(CXXFLAGS) -D GF2_bits=32 -D PARALLEL=1 -c -o $@ $^

%16PAR.o: %.cc
	$(CXX) $(CXXFLAGS) -D GF2_bits=16 -D PARALLEL=1 -c -o $@ $^

%16NOVEC.o: %.cc
	$(CXX) $(CXXFLAGS) -D GF2_bits=16 -D NOVEC=1 -c -o $@ $^

%0PAR.o: %.cc
	$(CXX) $(CXXFLAGS) -D GF2_bits=0 -D PARALLEL=1 -c -o $@ $^
