
DEFAULT_TARGET := DataBlocks2DTest.exe

L_SRC_DIR := .

L_DEFS := DOASSERT

L_EXT_INC_DIRS := ../../inc

L_EXT_OBJS := ../../obj/src/Data/DataBlocks2D.lo ../../obj/src/Data/DataVector.lo

PACKAGES := CppUnitTest esysUtils python23Static boost

include $(ESYS_ROOT)/make/Makefile.default

