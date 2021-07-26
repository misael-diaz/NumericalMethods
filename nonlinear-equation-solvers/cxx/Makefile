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

$(TEST_BRACKETING_EXE): $(OBJECTS)
	$(CXX) $(CXXOPT) $(OBJECTS) -o $(TEST_BRACKETING_EXE)\
	       	-Wl,-rpath -Wl,$(LIBSTDCXX_DIR)

$(NL_SOLVERS_OBJ): $(NL_SOLVERS_MOD)
	$(CXX) $(CXXOPT) -c $(NL_SOLVERS_MOD) -o $(NL_SOLVERS_OBJ)

$(TEST_BRACKETING_OBJ): $(NL_SOLVERS_MOD) $(TEST_BRACKETING_CXX)
	$(CXX) $(CXXOPT) -c $(TEST_BRACKETING_CXX) -o $(TEST_BRACKETING_OBJ)

clean:
	/bin/rm -rf $(OBJECTS) $(TESTS) $(GCM_CACHE)
