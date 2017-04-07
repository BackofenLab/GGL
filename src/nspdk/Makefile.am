
# additional include paths necessary to compile the C++ programs
AM_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/@LOCAL_SGD_FOLDER@ @AM_CXXFLAGS@

################################################################
# the NSPDK library assemblies 
################################################################

lib_LIBRARIES = libnspdk.a

################################################################
# the NSPDK library headers to distribute
################################################################

lib_nspdk_HEADERS = \
					nspdk/Utility.h \
					nspdk/SVector.h \
					nspdk/BaseGraphClass.h \
					nspdk/GraphClass.h \
					nspdk/FlagsService.h \
					nspdk/Histogram.h \
					nspdk/FeatureGenerator.h \
					nspdk/NSPDK_FeatureGenerator.h
					 


# header install directory via dummy lib
lib_nspdkdir = $(includedir)/nspdk


################################################################
# the NSPDK library objects
################################################################

libnspdk_a_SOURCES = \
					$(lib_nspdk_HEADERS) \
					nspdk/BaseGraphClass.cc \
					nspdk/GraphClass.cc \
					nspdk/FlagsService.cc \
					nspdk/Histogram.cc \
					nspdk/NSPDK_FeatureGenerator.cc
					
					
#					nspdk/NSPDK_quickcluster.cc \
					nspdk/NSPDK.cc \
					nspdk/NSPDK_removed_stems.cc \
					nspdk/NSPDK_weighted_edges.cc 


################################################################








################################################################

SUBDIRS = 
