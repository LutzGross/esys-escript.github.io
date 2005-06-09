# $Id$

# escript

PKG_CC_INC_DIRS  += $(ESYS_ROOT)/escript/inc
PKG_CPP_INC_DIRS += $(ESYS_ROOT)/escript/inc

PKG_LD_LIB_DIRS += $(ESYS_ROOT)/escript/lib

PKG_LD_LIBS += escriptcpp

# $Log$
# Revision 1.2  2005/06/09 05:38:03  jgs
# Merge of development branch back to main trunk on 2005-06-09
#
# Revision 1.1.2.2  2005/05/17 01:31:32  gross
# fedora_gcc compile environment added
#
#
