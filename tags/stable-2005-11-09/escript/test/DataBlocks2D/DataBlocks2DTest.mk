
DODEBUG := YES

DEFAULT_TARGET := DataBlocks2DTest.exe

L_SRC_DIR := .

L_EXT_INC_DIRS := ../../inc

L_EXT_OBJS := ../../obj/src/Data/DataBlocks2D.lo ../../obj/src/Data/DataVector.lo ../../obj/src/Data/Taipan.lo

PACKAGES := CppUnitTest esysUtils pythonStatic boost

include $(ESYS_ROOT)/make/Makefile.default