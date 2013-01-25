lib_LTLIBRARIES = libescriptreader.la

AM_CPPFLAGS = -I@top_srcdir@/src $(ESCRIPT_CFLAGS) $(NETCDF_CFLAGS) $(SILO_CFLAGS) 

libescriptreader_la_LIBADD = $(NETCDF_LIBS)

libescriptreader_la_SOURCES = DataVar.cpp \
                              ElementData.cpp \
                              Mesh.cpp \
                              MeshWithElements.cpp \
                              MPDataSet.cpp

libescriptreader_la_includedir = $(includedir)/escriptreader

libescriptreader_la_include_HEADERS = DataVar.h \
                                      ElementData.h \
                                      Mesh.h \
                                      MeshWithElements.h \
                                      MPDataSet.h \
                                      common.h

