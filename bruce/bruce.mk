# $Id$

# This makefile automatically finds, calculates dependencies, compiles
# and makes a static library from all C++ files found in the directories
# under ${L_SRC_DIR} (see below). Files and even whole directory trees
# under ${L_SRC_DIR} can be skipped by postfixing the name with a "-".

# Matt Davies 04/08/03

DEFAULT_TARGET := libbrucecpp.so pyc

ifeq ($(strip $(L_PYTH_DIR)),)
   L_PYTH_DIR := $(ESYS_ROOT)/esys
endif

L_PYTH_DIR := $(L_PYTH_DIR)/bruce

INSTALL_PYTH := ./lib/py_src/__init__.pyc ./lib/py_src/bruce.pyc

INSTALL_LIB := ./lib/libbrucecpp.so

PACKAGES := escript esysUtils python boostStatic

L_SRC_DIR:= ./src ./py_src

include $(ESYS_ROOT)/make/Makefile.default
