# $Id$

DODEBUG := YES

DEFAULT_TARGET := MeshAdapterTest.exe

L_SRC_DIR := .

PACKAGES := esysUtils finley escript CppUnitTest python23Static boost scsl141pre mmio

include $(ESYS_ROOT)/make/Makefile.default
