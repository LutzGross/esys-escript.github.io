# $Id$

DODEBUG := YES

DEFAULT_TARGET := SystemMatrixAdapterTest.exe

L_SRC_DIR := .

PACKAGES := esysUtils finley escript CppUnitTest pythonStatic boost scsl141pre mmio

include $(ESYS_ROOT)/make/Makefile.default
