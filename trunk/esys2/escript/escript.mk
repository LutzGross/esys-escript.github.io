# $Id$

# This makefile automatically finds, calculates dependencies, compiles
# and makes a static library from all C++ files found in the directories
# under ${L_SRC_DIR} (see below). Files and even whole directory trees
# under ${L_SRC_DIR} can be skipped by postfixing the name with a "-".

# Matt Davies 04/08/03

DEFAULT_TARGET := libescriptcpp.so pyc

INSTALL_PYTH := ./lib/py_src/__init__.pyc ./lib/py_src/linearPDEs.pyc ./lib/py_src/pdetools.pyc ./lib/py_src/util.pyc ./lib/py_src/escript.pyc ./lib/py_src/modelframe.pyc ./lib/py_src/timeseries.pyc # ./lib/py_src/modellib/*.pyc

INSTALL_LIB := ./lib/libescriptcpp.so

PACKAGES := esysUtils python23 boostStatic

L_SRC_DIR:= ./src ./py_src ./py_src/modellib

include $(ESYS_ROOT)/make/Makefile.default
