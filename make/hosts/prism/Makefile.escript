# $Id$

# escript

PKG_CC_INC_DIRS  += $(ESYS_ROOT)/escript/inc
PKG_CPP_INC_DIRS += $(ESYS_ROOT)/escript/inc

PKG_LD_LIB_DIRS += $(ESYS_ROOT)/escript/lib

PKG_LD_LIBS += escriptcpp

# $Log$
# Revision 1.2  2005/09/15 03:44:35  jgs
# Merge of development branch dev-02 back to main trunk on 2005-09-15
#
# Revision 1.1.2.1  2005/09/15 01:07:11  jgs
# added host specific makefiles and environment setup scripts for the prism
#
# Revision 1.2  2005/06/09 05:38:04  jgs
# Merge of development branch back to main trunk on 2005-06-09
#
# Revision 1.1.2.1  2005/05/16 03:23:13  imran
# Added generic host 'gcc'
#
# Revision 1.1.2.1  2004/12/17 02:50:44  imran
# added files for host 'darra'
# compile w/ gcc by default
#
# Revision 1.1.1.1.2.1  2004/12/09 06:22:21  jgs
# *** empty log message ***
#
# Revision 1.1.1.1  2004/10/26 06:53:58  jgs
# initial import of project esys2
#
# Revision 1.1.1.1.2.2  2004/09/28 07:03:14  jgs
# *** empty log message ***
#
# Revision 1.1.1.1.2.1  2004/09/27 06:37:54  jgs
# restructured make files
#
# Revision 1.1.1.1  2004/06/24 04:00:39  johng
# Initial version of eys using boost-python.
#
# Revision 1.1  2003/09/11 02:03:52  davies
# Added makefile configurations for several platforms.
#
