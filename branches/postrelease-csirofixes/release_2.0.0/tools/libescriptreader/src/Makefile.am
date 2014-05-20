SUBDIRS = escriptreader

bin_PROGRAMS =

if BUILD_ESCRIPT2SILO
bin_PROGRAMS += escript2silo
AM_CPPFLAGS = -I@srcdir@ $(ESCRIPT_CFLAGS) $(NETCDF_CFLAGS) $(SILO_CFLAGS)
escript2silo_LDADD = escriptreader/libescriptreader.la $(SILO_LIBS)
escript2silo_SOURCES = escript2silo.cpp
endif

if BUILD_ESDCREATE
bin_PROGRAMS += esdcreate
esdcreate_SOURCES = esdcreate.cpp
endif

