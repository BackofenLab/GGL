
# additional include paths necessary to compile the C++ programs
AM_CXXFLAGS = -I$(top_srcdir)/src \
				-I$(top_builddir)/src/ \
				-I$(top_srcdir)/src/@LOCAL_VF2_FOLDER@ \
				-I$(top_srcdir)/src/@LOCAL_SGD_FOLDER@ \
				-I$(top_srcdir)/src/@LOCAL_NSPDK_FOLDER@ \
				@AM_CXXFLAGS@

################################################################
# the SGM library assemblies 
################################################################

lib_LIBRARIES = libsgm.a libsgm_vf2.a libsgm_nspdk.a


################################################################
# the SGM CORE library headers to distribute
################################################################

lib_sgm_HEADERS = \
					HashMap.hh \
					Graph_Interface.hh \
					Graph_Interface.icc \
					Pattern.hh \
					Pattern_Automorphism.hh \
					PA_OrderCheck.hh \
					PA_OrderCheck.icc \
					Graph_boost.hh \
					Graph_boost.icc \
					Graph_boostV_p.hh \
					Graph_boostV_p.icc \
					SubGraph.hh \
					SubGraph.icc \
					Match.hh \
					MC_Node.hh \
					MC_Edge.hh \
					Match_Reporter.hh \
					SubGraphMatching.hh \
					GraphMatching.hh \
					MR_Counting.hh \
					MR_Storing.hh \
					MR_SymmBreak.hh \
					MR_stream.hh \
					RingReporter.hh \
					RingPerception.hh \
					RP_Hanser96.hh \
					GraphScaffold.hh \
					Graph_NodeLabelPrefix.hh \
					MR_restoreNodeLabel.hh \
					MR_restoreNodeLabel.icc

# header install directory via dummy lib
lib_sgmdir = $(includedir)/sgm


################################################################
# the SGM CORE library objects
################################################################

libsgm_a_SOURCES = \
					$(lib_sgm_HEADERS) \
					Graph_Interface.cc \
					Pattern.cc \
					SubGraph.cc \
					RingReporter.cc \
					RP_Hanser96.cc \
					GraphScaffold.cc
					

################################################################





################################################################
# the SGM VF2 library port 
################################################################

################################################################
# the SGM VF2 library port headers to distribute
################################################################

lib_sgm_vf2_HEADERS = \
					VF2_MatchingHandler.hh \
					VF2_MatchingHandler.icc \
					SGM_vf2.hh \
					GM_vf2.hh 

# header install directory via dummy lib
lib_sgm_vf2dir = $(includedir)/sgm


################################################################
# the SGM VF2 library port objects
################################################################

libsgm_vf2_a_SOURCES = \
					$(lib_sgm_vf2_HEADERS) \
					VF2_MatchingHandler.cc \
					SGM_vf2.cc \
					GM_vf2.cc
					

################################################################






################################################################
# the SGM NSPDK library port 
################################################################

################################################################
# the SGM NSPDK library port headers to distribute
################################################################

lib_sgm_nspdk_HEADERS = \
					NSPDK_port.hh
					 


# header install directory via dummy lib
lib_sgm_nspdkdir = $(includedir)/sgm


################################################################
# the SGM NSPDK library port objects
################################################################

libsgm_nspdk_a_SOURCES = \
					$(lib_sgm_nspdk_HEADERS) \
					NSPDK_port.cc
					

################################################################








################################################################

SUBDIRS = 
