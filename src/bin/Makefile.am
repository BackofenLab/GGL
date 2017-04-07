
# additional libs and include paths necessary to compile the C++ programs
AM_CXXFLAGS = 	-I$(top_srcdir)/src \
				-I$(top_builddir)/src/ \
				-I$(top_srcdir)/src/@LOCAL_VF2_FOLDER@ \
				-I$(top_srcdir)/src/@LOCAL_SGD_FOLDER@ \
				-I$(top_srcdir)/src/@LOCAL_NSPDK_FOLDER@ \
				@AM_CXXFLAGS@

AM_LDFLAGS  = @AM_LDFLAGS@ 
LIBS        = @LIBS@ 

bin_PROGRAMS = 	gmlRings \
				aromModelNSPDK aromModel2impl \
				molTool reactionTool toyChem \
				molMatch molRings \
				gameOfLife sudoku 
# catalan hanser # disabled due to a BGL bug for v < 1.50.0 
# ggrIt

#############   TOYCHEM PROGRAM   ##################################

toyChem_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					biu/Timer.hh \
					version.hh \
					RRC_QuantumMechanics.hh \
					RRC_QuantumMechanics.cc \
					toyChemUtil.hh \
					toyChemUtil.cc \
					toyChem.cc
				
toyChem_LDADD = \
				$(top_builddir)/src/ggl/libgglchem.a \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_nspdk.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_NSPDK_FOLDER@/libnspdk.a \
				$(top_builddir)/src/@LOCAL_SGD_FOLDER@/libsgd.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

####################################################################

###################   MOLECULE CHECKER   ###########################

molTool_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					version.hh \
					toyChemUtil.hh \
					toyChemUtil.cc \
					molTool.cc
				
molTool_LDADD = \
				$(top_builddir)/src/ggl/libgglchem.a \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_nspdk.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_NSPDK_FOLDER@/libnspdk.a \
				$(top_builddir)/src/@LOCAL_SGD_FOLDER@/libsgd.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

####################################################################

###################   MOLECULE MATCHER   ###########################

molMatch_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					version.hh \
					toyChemUtil.hh \
					toyChemUtil.cc \
					molMatch.cc
				
molMatch_LDADD = \
				$(top_builddir)/src/ggl/libgglchem.a \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_nspdk.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_NSPDK_FOLDER@/libnspdk.a \
				$(top_builddir)/src/@LOCAL_SGD_FOLDER@/libsgd.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

####################################################################

###################   GML RING FINDER   ###########################

gmlRings_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					version.hh \
					gmlRings.cc
				
gmlRings_LDADD = \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

####################################################################

###################   MOLECULE RING FINDER   ###########################

molRings_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					version.hh \
					toyChemUtil.hh \
					toyChemUtil.cc \
					molRings.cc
				
molRings_LDADD = \
				$(top_builddir)/src/ggl/libgglchem.a \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_nspdk.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_NSPDK_FOLDER@/libnspdk.a \
				$(top_builddir)/src/@LOCAL_SGD_FOLDER@/libsgd.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

####################################################################

###################   REACTION CHECKER   ###########################

reactionTool_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					version.hh \
					toyChemUtil.hh \
					toyChemUtil.cc \
					reactionTool.cc
				
reactionTool_LDADD = \
				$(top_builddir)/src/ggl/libgglchem.a \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_nspdk.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_NSPDK_FOLDER@/libnspdk.a \
				$(top_builddir)/src/@LOCAL_SGD_FOLDER@/libsgd.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

####################################################################

#############   AROMATICITY MODEL PROGRAM   ########################

aromModelNSPDK_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					biu/Timer.hh \
					version.hh \
					aromModelNSPDK.cc
				
aromModelNSPDK_LDADD = \
				$(top_builddir)/src/ggl/libgglchem.a \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_nspdk.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_NSPDK_FOLDER@/libnspdk.a \
				$(top_builddir)/src/@LOCAL_SGD_FOLDER@/libsgd.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

####################################################################

###################   MODEL CONVERTION   ###########################

aromModel2impl_SOURCES =	\
					aromModel2impl.cc
				
aromModel2impl_LDADD = \
				$(top_builddir)/src/@LOCAL_SGD_FOLDER@/libsgd.a

####################################################################

#######################   GAME OF LIFE   ###########################

gameOfLife_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					version.hh \
					gameOfLife.cc
				
gameOfLife_LDADD = \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

####################################################################

#########################   SUDOKU   ###############################

sudoku_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					version.hh \
					sudoku.cc
				
sudoku_LDADD =  \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

####################################################################

#########################   CATALAN   ##############################

#catalan_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					version.hh \
					catalan.cc
				
#catalan_LDADD =  \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 


####################################################################

#########################   HANSER   ###############################

#hanser_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					version.hh \
					hanser.cc
				
#hanser_LDADD =  \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

####################################################################

#############   GRAPH GRAMMAR ITERATOR PROGRAM   ###################
#
# ggrIt_SOURCES =	\
#					biu/OptionParser.hh \
#					biu/OptionParser.cc \
#					version.hh \
#					ggrIt.cc
#
# ggrIt_LDADD =  \
#				$(top_builddir)/src/ggl/libggl.a \
#				$(top_builddir)/src/sgm/libsgm_vf2.a \
#				$(top_builddir)/src/sgm/libsgm.a \
#				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 
#
####################################################################

