# $Id$

DODEBUG := YES

DEFAULT_TARGET := DataTest.exe

L_SRC_DIR := .

L_DEFS := DOASSERT

PACKAGES := escript CppUnitTest esysUtils boostStatic python23Static

include $(ESYS_ROOT)/make/Makefile.default
