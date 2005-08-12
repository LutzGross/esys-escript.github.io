
DODEBUG := YES

DEFAULT_TARGET := TaipanTest.exe

L_SRC_DIR := .

L_DEFS := DOASSERT

PACKAGES := CppUnitTest esysUtils pythonStatic boost

L_EXT_INC_DIRS := ../../inc

L_EXT_OBJS := ../../obj/src/Data/Taipan.lo

include $(ESYS_ROOT)/make/Makefile.default

