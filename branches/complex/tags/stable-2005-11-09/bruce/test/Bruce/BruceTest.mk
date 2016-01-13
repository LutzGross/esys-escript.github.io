
DODEBUG := YES

DEFAULT_TARGET := BruceTest.exe

L_SRC_DIR := .

PACKAGES := bruce escript CppUnitTest esysUtils boostStatic pythonStatic

include $(ESYS_ROOT)/make/Makefile.default
