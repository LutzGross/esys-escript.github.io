
DODEBUG := YES

DEFAULT_TARGET := BruceFactoryTest.exe

L_SRC_DIR := .

PACKAGES := bruce escript CppUnitTest esysUtils boostStatic pythonStatic

include $(ESYS_ROOT)/make/Makefile.default
