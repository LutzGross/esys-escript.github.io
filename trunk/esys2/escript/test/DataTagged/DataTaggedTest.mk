# $Id$

DODEBUG := YES

DEFAULT_TARGET := DataTaggedTest.exe

L_SRC_DIR := .

L_DEFS := DOASSERT

PACKAGES := escript CppUnitTest esysUtils pythonStatic boost

include $(ESYS_ROOT)/make/Makefile.default
