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
	$(CC) $(CCOPT) $(OBJECTS) $(TEST_NUMINT_OBJ) -o $(TEST_NUMINT_EXE)\
              $(LIBS)

$(TEST_RAMP_EXE): $(OBJECTS) $(TEST_RAMP_OBJ)
	$(CC) $(CCOPT) $(OBJECTS) $(TEST_RAMP_OBJ) -o $(TEST_RAMP_EXE)\
	      $(LIBS)

$(NL_SOLVER_FZERO_OBJ): $(NL_SOLVER_FZERO_SRC)
	$(CC) $(CCOPT) -c $(NL_SOLVER_FZERO_SRC) -o $(NL_SOLVER_FZERO_OBJ)

$(NUM_INTEGRATION_OBJ): $(NUM_INTEGRATION_SRC)
	$(CC) $(CCOPT) -c $(NUM_INTEGRATION_SRC) -o $(NUM_INTEGRATION_OBJ)

$(TEST_NUMINT_OBJ): $(NUM_INTEGRATION_SRC) $(TEST_NUMINT_SRC)
	$(CC) $(CCOPT) -c $(TEST_NUMINT_SRC) -o $(TEST_NUMINT_OBJ)

clean:
	/bin/rm -rf *.o $(TESTS)
