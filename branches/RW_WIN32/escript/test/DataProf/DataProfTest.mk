
DODEBUG := YES

DEFAULT_TARGET := DataProfTest.exe

L_SRC_DIR := .

PACKAGES := CppUnitTest esysUtils pythonStatic boost

L_EXT_INC_DIRS := ../../inc

L_EXT_OBJS := ../../obj/src/Data/DataProf.lo

include $(ESYS_ROOT)/make/Makefile.default

