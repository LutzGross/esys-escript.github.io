# $Id$

DODEBUG := YES

DEFAULT_TARGET := DataTest.exe

L_SRC_DIR := .

PACKAGES := escript CppUnitTest esysUtils boostStatic pythonStatic

include $(ESYS_ROOT)/make/Makefile.default
