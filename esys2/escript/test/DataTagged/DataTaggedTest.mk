# $Id$

DODEBUG := YES

DEFAULT_TARGET := DataTaggedTest.exe

L_SRC_DIR := .

L_DEFS := DOASSERT

PACKAGES := escript CppUnitTest esysUtils python23Static boost

include $(ESYS_ROOT)/make/Makefile.default
