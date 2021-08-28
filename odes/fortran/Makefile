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

$(TEST_NUMINT_EXE): $(OBJECTS) $(TEST_NUMINT_OBJ)
	$(FC) $(FCOPT) $(OBJECTS) $(TEST_NUMINT_OBJ) -o $(TEST_NUMINT_EXE)\
                       $(LIBS)

$(TEST_IMPULSE_EXE): $(OBJECTS) $(TEST_IMPULSE_OBJ)
	$(FC) $(FCOPT) $(OBJECTS) $(TEST_IMPULSE_OBJ) -o $(TEST_IMPULSE_EXE)\
                       $(LIBS)

$(NL_SOLVER_FZERO_OBJ): $(NL_SOLVER_FZERO_SRC)
	$(FC) $(FCOPT) -c $(NL_SOLVER_FZERO_SRC) -o $(NL_SOLVER_FZERO_OBJ)

$(NUM_INTEGRATORS_OBJ): $(NL_SOLVER_FZERO_OBJ) $(NUM_INTEGRATORS_SRC)
	$(FC) $(FCOPT) -c $(NUM_INTEGRATORS_SRC) -o $(NUM_INTEGRATORS_OBJ)

$(TEST_NUMINT_OBJ): $(NUM_INTEGRATORS_SRC) $(TEST_NUMINT_SRC)
	$(FC) $(FCOPT) -c $(TEST_NUMINT_SRC) -o $(TEST_NUMINT_OBJ)

$(TEST_IMPULSE_OBJ): $(NUM_INTEGRATORS_SRC) $(TEST_IMPULSE_SRC)
	$(FC) $(FCOPT) -c $(TEST_IMPULSE_SRC) -o $(TEST_IMPULSE_OBJ)

clean:
	/bin/rm -rf *.o *.mod $(TESTS)
