# $Id$

# This makefile automatically finds, calculates dependencies, compiles
# and makes a static library from all C++ files found in the directories
# under ${L_SRC_DIR} (see below). Files and even whole directory trees
# under ${L_SRC_DIR} can be skipped by postfixing the name with a "-".

# Matt Davies 04/08/03

DEFAULT_TARGET := libfinleyC.a libfinley.so

PACKAGES := escript mmio esysUtils python23 boost

L_SRC_DIR := ./src/finleyC ./src/finleyC/Solvers ./src/finleyC/SCSL ./src/CPPAdapter

L_DEFS := ITERATIVE_SOLVER=NO_LIB DIRECT_SOLVER=NO_LIB

include $(ESYS_ROOT)/make/Makefile.default

