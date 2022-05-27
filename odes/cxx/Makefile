#!/usr/bin/make
#
# source: Makefile
# author: misael-diaz
# date:   2021/07/20
#
# Synopsis:
# Defines the Makefile for building the program with GNU make.
#
# Copyright (c) 2021 Misael Diaz-Maldonado
#
# This file is released under the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

include make-inc

all: $(TESTS)

$(TEST_NUM_INTEGRATORS_EXE): $(OBJECTS) $(TEST_NUM_INTEGRATORS_OBJ)
	$(CXX) $(CXXOPT) $(OBJECTS) $(TEST_NUM_INTEGRATORS_OBJ) -o\
	       $(TEST_NUM_INTEGRATORS_EXE) -Wl,-rpath -Wl,$(LIBSTDCXX_DIR)

$(TEST_RAMP_EXE): $(OBJECTS) $(TEST_RAMP_OBJ)
	$(CXX) $(CXXOPT) $(OBJECTS) $(TEST_RAMP_OBJ) -o\
	       $(TEST_RAMP_EXE) -Wl,-rpath -Wl,$(LIBSTDCXX_DIR)

$(TEST_IMPULSE_EXE): $(OBJECTS) $(TEST_IMPULSE_OBJ)
	$(CXX) $(CXXOPT) $(OBJECTS) $(TEST_IMPULSE_OBJ) -o\
	       $(TEST_IMPULSE_EXE) -Wl,-rpath -Wl,$(LIBSTDCXX_DIR)

$(TEST_SINUSOID_EXE): $(OBJECTS) $(TEST_SINUSOID_OBJ)
	$(CXX) $(CXXOPT) $(OBJECTS) $(TEST_SINUSOID_OBJ) -o\
	       $(TEST_SINUSOID_EXE) -Wl,-rpath -Wl,$(LIBSTDCXX_DIR)

$(TEST_KINETICS_EXE): $(OBJECTS) $(TEST_KINETICS_OBJ)
	$(CXX) $(CXXOPT) $(OBJECTS) $(TEST_KINETICS_OBJ) -o\
	       $(TEST_KINETICS_EXE) -Wl,-rpath -Wl,$(LIBSTDCXX_DIR)

$(NL_SOLVER_FZERO_OBJ): $(NL_SOLVER_FZERO_MOD)
	$(CXX) $(CXXOPT) -c $(NL_SOLVER_FZERO_MOD) -o $(NL_SOLVER_FZERO_OBJ)

$(NUM_INTEGRATORS_OBJ): $(NL_SOLVER_FZERO_OBJ) $(NUM_INTEGRATORS_MOD)
	$(CXX) $(CXXOPT) -c $(NUM_INTEGRATORS_MOD) -o $(NUM_INTEGRATORS_OBJ)

$(TEST_NUM_INTEGRATORS_OBJ): $(NUM_INTEGRATORS_MOD)\
       	                     $(TEST_NUM_INTEGRATORS_CXX)
	$(CXX) $(CXXOPT) -c $(TEST_NUM_INTEGRATORS_CXX) -o\
	       	            $(TEST_NUM_INTEGRATORS_OBJ)

$(TEST_RAMP_OBJ): $(NUM_INTEGRATORS_MOD) $(TEST_RAMP_CXX)
	$(CXX) $(CXXOPT) -c $(TEST_RAMP_CXX) -o $(TEST_RAMP_OBJ)

$(TEST_IMPULSE_OBJ): $(NUM_INTEGRATORS_MOD) $(TEST_IMPULSE_CXX)
	$(CXX) $(CXXOPT) -c $(TEST_IMPULSE_CXX) -o $(TEST_IMPULSE_OBJ)

$(TEST_SINUSOID_OBJ): $(NUM_INTEGRATORS_MOD) $(TEST_SINUSOID_CXX)
	$(CXX) $(CXXOPT) -c $(TEST_SINUSOID_CXX) -o $(TEST_SINUSOID_OBJ)

$(TEST_KINETICS_OBJ): $(NUM_INTEGRATORS_MOD) $(TEST_KINETICS_CXX)
	$(CXX) $(CXXOPT) -c $(TEST_KINETICS_CXX) -o $(TEST_KINETICS_OBJ)

clean:
	/bin/rm -rf *.o $(TESTS) $(GCM_CACHE)
