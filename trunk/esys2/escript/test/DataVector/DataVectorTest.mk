
DODEBUG := YES

DEFAULT_TARGET := DataVectorTest.exe

L_SRC_DIR := .

PACKAGES := CppUnitTest esysUtils pythonStatic boost

L_EXT_INC_DIRS := ../../inc

L_EXT_OBJS := ../../obj/src/Data/DataVector.lo ../../obj/src/Data/Taipan.lo

include $(ESYS_ROOT)/make/Makefile.default
