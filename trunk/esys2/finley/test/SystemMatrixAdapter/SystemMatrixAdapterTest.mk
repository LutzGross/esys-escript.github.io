# $Id$

DODEBUG = YES

DEFAULT_TARGET := SystemMatrixAdapterTest.exe

L_SRC_DIR := .

L_DEFS := DOASSERT

PACKAGES := esysUtils finley escript CppUnitTest python23Static boost scsl141pre mmio

include $(ESYS_ROOT)/make/Makefile.default

