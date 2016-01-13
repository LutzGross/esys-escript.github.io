
DODEBUG := YES 

DEFAULT_TARGET := DataAlgorithmAdapterTest.exe

L_SRC_DIR := .

PACKAGES := escript CppUnitTest esysUtils pythonStatic boost

include $(ESYS_ROOT)/make/Makefile.default
