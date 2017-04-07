


################################################################
# the local SGD library assembly (based on SGD-1.3)
################################################################

lib_LIBRARIES = libsgd.a


################################################################
# the VF library headers to distribute
################################################################

lib_sgd_HEADERS = \
					./sgd/timer.h \
					./sgd/wrapper.h \
					./sgd/vectors.h \
					./sgd/matrices.h \
					./sgd/svmmodel.h \
					./sgd/svmsgd.h

# header install directory via dummy lib
lib_sgddir = $(includedir)/sgd


################################################################
# the SGD library objects
################################################################

libsgd_a_SOURCES = \
					$(lib_sgd_HEADERS) \
					./sgd/timer.cpp \
					./sgd/vectors.cpp \
					./sgd/matrices.cpp \
					./sgd/svmsgd.cpp

################################################################

AM_LDFLAGS = 
LIBS =

bin_PROGRAMS = svmsgd

svmsgd_SOURCES = ./bin/run_svmsgd.cc
svmsgd_LDADD = ./libsgd.a

SUBDIRS = 
