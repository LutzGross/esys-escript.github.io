
DEFAULT_TARGET := SystemMatrixAdapterTest.exe

L_SRC_DIR := .

L_DEFS := DOASSERT

PACKAGES := finley escript CppUnitTest esysUtils python23Static boost

include $(ESYS_ROOT)/make/Makefile.default

