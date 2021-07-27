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
	$(CC) $(CCOPT) $(OBJECTS) -o $(TEST_BRACKETING_EXE) $(LIBS)

$(NL_SOLVERS_OBJ): $(NL_SOLVERS_SRC)
	$(CC) $(CCOPT) -c $(NL_SOLVERS_SRC) -o $(NL_SOLVERS_OBJ)

$(TEST_BRACKETING_OBJ): $(NL_SOLVERS_SRC) $(TEST_BRACKETING_SRC)
	$(CC) $(CCOPT) -c $(TEST_BRACKETING_SRC) -o $(TEST_BRACKETING_OBJ)

clean:
	/bin/rm -rf *.o $(TESTS)
