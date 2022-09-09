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

$(TEST_VECTOR_EXE): $(OBJECTS)
	$(CC) $(CCOPT) $(OBJECTS) -o $(TEST_VECTOR_EXE) $(LIBS)

$(VECTOR_OBJ): $(VECTOR_SRC)
	$(CC) $(CCOPT) -c $(VECTOR_SRC) -o $(VECTOR_OBJ)

$(TEST_VECTOR_OBJ): $(VECTOR_SRC) $(TEST_VECTOR_SRC)
	$(CC) $(CCOPT) -c $(TEST_VECTOR_SRC) -o $(TEST_VECTOR_OBJ)

clean:
	/bin/rm -rf *.o $(TESTS)
