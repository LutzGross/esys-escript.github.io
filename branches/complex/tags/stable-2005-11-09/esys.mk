# $Id$

# no default implicit rules
.SUFFIX:

.PHONY: all clean install

MODULES := escript esysUtils finley tools

RECURSIVE_TARGETS := all-recursive clean-recursive install-recursive

all: all-recursive

clean: clean-recursive

install: install-recursive

${RECURSIVE_TARGETS}:
	@TARGET=`echo $@ | sed s/-recursive//`; \
	for MODULE in ${MODULES}; do \
		pwd; \
		${MAKE} ${MAKEFLAGS} -C $${MODULE} $${TARGET}; \
	done

# $Log$
# Revision 1.1  2004/10/26 06:53:54  jgs
# Initial revision
#
# Revision 1.1  2004/09/23 00:13:29  jgs
# rewrote mk, new Makefile->esys.mk regime
#
# Revision 1.2  2004/07/19 01:12:01  johng
# Added new functions to Data (such as Lsup).
#
# Revision 1.1.1.1  2004/06/24 04:00:38  johng
# Initial version of eys using boost-python.
#
