# $Id$

# This makefile automatically finds, calculates dependencies, compiles
# and makes a static library from all C++ files found in the directories
# under ${L_SRC_DIR} (see below). Files and even whole directory trees
# under ${L_SRC_DIR} can be skipped by postfixing the name with a "-".

# Matt Davies 04/08/03

DODEBUG = NO

DEFAULT_TARGET := libfinleyC.a libfinleycpp.so

INSTALL_PYTH := ./py_src/finley.py

INSTALL_LIB := ./lib/libfinleycpp.so

L_DEFS := ITERATIVE_SOLVER=NO_LIB DIRECT_SOLVER=SGI_SCSL
# L_DEFS := ITERATIVE_SOLVER=NO_LIB DIRECT_SOLVER=NO_LIB

PACKAGES := escript mmio esysUtils python23 boost scsl141pre

L_SRC_DIR := ./src/finleyC ./src/finleyC/Solvers ./src/finleyC/SCSL ./src/CPPAdapter

#L_PYTH_DIR := /raid2/tools/esys/esys

#L_INSTLIB_DIR := /raid2/tools/esys/lib

include $(ESYS_ROOT)/make/Makefile.default

