
DODEBUG := YES 

DEFAULT_TARGET := DataAlgorithmAdapterTest.exe

L_SRC_DIR := .

PACKAGES := escript CppUnitTest esysUtils python23Static boost

L_DEFS := DOASSERT

include $(ESYS_ROOT)/make/Makefile.default
