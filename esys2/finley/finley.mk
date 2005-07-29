# $Id$

# This makefile automatically finds, calculates dependencies, compiles
# and makes a static library from all C++ files found in the directories
# under ${L_SRC_DIR} (see below). Files and even whole directory trees
# under ${L_SRC_DIR} can be skipped by postfixing the name with a "-".

# Matt Davies 04/08/03

DEFAULT_TARGET := libfinleyC.a libfinleycpp.so pyc

# Python install directory
ifeq ($(strip $(L_PYTH_DIR)),)
   L_PYTH_DIR := $(ESYS_ROOT)/esys
endif

L_PYTH_DIR := $(L_PYTH_DIR)/finley

INSTALL_PYTH := ./lib/py_src/__init__.pyc ./lib/py_src/finley.pyc

INSTALL_LIB := ./lib/libfinleycpp.so

# These settings are now specified in Makefile.host
#L_DEFS := ITERATIVE_SOLVER=NO_LIB DIRECT_SOLVER=NO_LIB

PACKAGES := escript mmio esysUtils python23 boost scsl141pre

L_SRC_DIR := ./src/finleyC ./src/finleyC/Solvers ./src/finleyC/SCSL ./src/CPPAdapter ./py_src

include $(ESYS_ROOT)/make/Makefile.default
