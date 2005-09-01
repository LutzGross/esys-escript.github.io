# $Id$

# This makefile automatically finds, calculates dependencies, compiles
# and makes a static library from all C++ files found in the directories
# under ${L_SRC_DIR} (see below). Files and even whole directory trees
# under ${L_SRC_DIR} can be skipped by postfixing the name with a "-".

# Matt Davies 04/08/03

DEFAULT_TARGET := libfinleyC.a libfinleycpp.so pyc

ifeq ($(strip $(L_PYTH_DIR)),)
   L_PYTH_DIR := $(ESYS_ROOT)/esys
endif

L_PYTH_DIR := $(L_PYTH_DIR)/finley

INSTALL_PYTH := ./lib/py_src/__init__.pyc ./lib/py_src/finley.pyc 

INSTALL_LIB := ./lib/libfinleycpp.so

PACKAGES := escript mmio esysUtils python boost scsl141pre

L_SRC_DIR := ./src/finleyC ./src/finleyC/Solvers ./src/finleyC/SCSL ./src/CPPAdapter ./py_src

# determine solver packages to use
ifeq ($(strip $(USE_GCC)),)
   HOSTNAME := ${shell hostname | sed -e 's/\..*//'}
else
   HOSTNAME := gcc
endif

ifeq ($(strip $(HOSTNAME)),ess)
   L_DEFS := ITERATIVE_SOLVER=NO_LIB DIRECT_SOLVER=SGI_SCSL
endif

ifeq ($(strip $(HOSTNAME)),gcc)
   L_DEFS := ITERATIVE_SOLVER=NO_LIB DIRECT_SOLVER=NO_LIB
endif

include $(ESYS_ROOT)/make/Makefile.default
