# $Id$

# Matt Davies 04/08/03

DEFAULT_TARGET := libpaso.a

INSTALL_LIB := ./lib/libpaso.a

PACKAGES := mmio scsl mkl umfpack

L_SRC_DIR := ./src ./src/Solvers ./src/SCSL 

include $(ESYS_ROOT)/make/Makefile.default
