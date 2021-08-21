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

$(TEST_NUM_INTEGRATORS_EXE): $(OBJECTS)
	$(CXX) $(CXXOPT) $(OBJECTS) -o $(TEST_NUM_INTEGRATORS_EXE)\
	       	-Wl,-rpath -Wl,$(LIBSTDCXX_DIR)

$(NL_SOLVER_FZERO_OBJ): $(NL_SOLVER_FZERO_MOD)
	$(CXX) $(CXXOPT) -c $(NL_SOLVER_FZERO_MOD) -o $(NL_SOLVER_FZERO_OBJ)

$(NUM_INTEGRATORS_OBJ): $(NL_SOLVER_FZERO_OBJ) $(NUM_INTEGRATORS_MOD)
	$(CXX) $(CXXOPT) -c $(NUM_INTEGRATORS_MOD) -o $(NUM_INTEGRATORS_OBJ)

$(TEST_NUM_INTEGRATORS_OBJ): $(NUM_INTEGRATORS_MOD)\
       	                     $(TEST_NUM_INTEGRATORS_CXX)
	$(CXX) $(CXXOPT) -c $(TEST_NUM_INTEGRATORS_CXX) -o\
	       	            $(TEST_NUM_INTEGRATORS_OBJ)

clean:
	/bin/rm -rf $(OBJECTS) $(TESTS) $(GCM_CACHE)
