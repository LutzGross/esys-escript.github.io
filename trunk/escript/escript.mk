# $Id$

# This makefile automatically finds, calculates dependencies, compiles
# and makes a static library from all C++ files found in the directories
# under ${L_SRC_DIR} (see below). Files and even whole directory trees
# under ${L_SRC_DIR} can be skipped by postfixing the name with a "-".

# Matt Davies 04/08/03

DEFAULT_TARGET := libescriptcpp.so pyc

# Python install directory
ifeq ($(strip $(L_PYTH_DIR)),)
   L_PYTH_DIR := $(ESYS_ROOT)/esys
endif

L_PYTH_DIR := $(L_PYTH_DIR)/escript

INSTALL_PYTH := ./lib/py_src/__init__.pyc ./lib/py_src/linearPDEs.pyc ./lib/py_src/pdetools.pyc ./lib/py_src/util.pyc ./lib/py_src/escript.pyc ./lib/py_src/modelframe.pyc ./lib/py_src/timeseries.pyc ./lib/py_src/esysXML.pyc ./lib/py_src/test_linearPDEs.pyc ./lib/py_src/test_util.pyc ./lib/py_src/symbols.pyc

INSTALL_LIB := ./lib/libescriptcpp.so

PACKAGES := esysUtils python boostStatic

L_SRC_DIR:= ./src ./py_src

include $(ESYS_ROOT)/make/Makefile.default
