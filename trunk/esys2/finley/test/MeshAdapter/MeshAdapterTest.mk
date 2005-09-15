# $Id$

DODEBUG := YES

DEFAULT_TARGET := MeshAdapterTest.exe

L_SRC_DIR := .

PACKAGES := esysUtils finley paso scsl escript CppUnitTest pythonStatic boost mmio mkl umfpack

include $(ESYS_ROOT)/make/Makefile.default
