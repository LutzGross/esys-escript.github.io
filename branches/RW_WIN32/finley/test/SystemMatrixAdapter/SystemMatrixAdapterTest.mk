# $Id$

DODEBUG := YES

DEFAULT_TARGET := SystemMatrixAdapterTest.exe

L_SRC_DIR := .

PACKAGES := esysUtils finley paso scsl escript CppUnitTest pythonStatic boost mmio

include $(ESYS_ROOT)/make/Makefile.default
