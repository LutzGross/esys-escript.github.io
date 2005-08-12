# $Id$

DODEBUG := YES

DEFAULT_TARGET := DataTest.exe

L_SRC_DIR := .

L_DEFS := DOASSERT

PACKAGES := escript CppUnitTest esysUtils boostStatic pythonStatic

include $(ESYS_ROOT)/make/Makefile.default
