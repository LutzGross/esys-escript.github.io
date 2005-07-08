
DODEBUG := YES

DEFAULT_TARGET := DataProfTest.exe

L_SRC_DIR := .

L_DEFS := DOASSERT

PACKAGES := CppUnitTest esysUtils python23Static boost

L_EXT_INC_DIRS := ../../inc

L_EXT_OBJS := ../../obj/src/Data/DataProf.lo

include $(ESYS_ROOT)/make/Makefile.default

