
# additional include paths necessary to compile the C++ programs
AM_CXXFLAGS = -I$(top_srcdir)/src \
				-I$(top_builddir)/src/ \
				-I$(top_srcdir)/src/@LOCAL_VF2_FOLDER@ \
				-I$(top_srcdir)/src/@LOCAL_SGD_FOLDER@ \
				-I$(top_srcdir)/src/@LOCAL_NSPDK_FOLDER@ \
				@AM_CXXFLAGS@

################################################################
# the GGL (graph grammar) library assemblies 
################################################################

lib_LIBRARIES = libggl.a libgglchem.a


################################################################
# the GGL CORE library headers to distribute
################################################################

lib_ggl_HEADERS = \
					Rule.hh \
					Rule.icc \
					RuleGraph.hh \
					RuleGraph.icc \
					Graph.hh \
					Graph_Storage.hh \
					GS_stream.hh \
					GS_stream.icc \
					GS_STL.hh \
					GS_STL.icc \
					MR_ApplyRule.hh \
					MR_ApplyRule.icc \
					MR_ApplyRule_NodeLabelPrefix.hh \
					MR_ApplyRule_NodeLabelPrefix.icc \
					DFS_ApplyRule.hh \
					DFS_ApplyRule.icc \
					Rule_GML_error.hh \
					Rule_GML_writer.hh \
					Rule_GML_writer.icc \
					Rule_GML_grammar.hh \
					Rule_GML_grammar.icc \
					Rule_GMLparser.hh \
					Graph_GML_grammar.hh \
					Graph_GML_grammar.icc \
					Graph_GMLparser.hh \
					Graph_GML_writer.hh \
					Graph_GML_writer.icc \
					Graph_GXL_writer.hh \
					Graph_GXL_writer.icc \
					Graph_gSpan_writer.hh

# header install directory via dummy lib
lib_ggldir = $(includedir)/ggl


################################################################
# the GGL CORE library objects
################################################################

libggl_a_SOURCES = \
					$(lib_ggl_HEADERS) \
					Graph.cc \
					Graph_GMLparser.cc \
					Graph_GXL_writer.cc \
					Graph_GML_writer.cc \
					Graph_gSpan_writer.cc \
					Rule.cc \
					RuleGraph.cc \
					Rule_GML_writer.cc \
					Rule_GMLparser.cc \
					MR_ApplyRule.cc \
					DFS_ApplyRule.cc \
					GS_stream.cc 
					


################################################################
# the GGL CHEMISTRY library headers to distribute
################################################################

# set the test routines for OpenBabel depending classes
gglchem_OB_headerfiles = \
						chem/MoleculeOB.hh \
						chem/MoleculeOB.icc \
						chem/SMILESwriterOB.hh \
						chem/GS_SMILES_OB.hh \
						chem/GS_SMILES_OB.icc

lib_gglchem_HEADERS = \
					chem/Molecule.hh \
					chem/Molecule_Graph_noClass.hh \
					chem/MoleculeUtil.hh \
					chem/MoleculeUtil.icc \
					chem/ClassIDGraph.hh \
					chem/AromaticityPerception.hh \
					chem/AromaticityPerception.icc \
					chem/AP_disabled.hh \
					chem/AP_NSPDK_Model.hh \
					chem/AP_NSPDK_Model_OpenBabel_2013.icc \
					chem/AP_NSPDK_Model_Marvin_general_2013.icc \
					chem/AP_NSPDK.hh \
					chem/AP_NSPDK.icc \
					chem/GS_chem.hh \
					chem/GS_chem.icc \
					chem/GS_MolCheck.hh \
					chem/GS_MolCheck.icc \
					chem/ChemRule.hh \
					chem/ChemRule.icc \
					chem/ChemRuleGraph.hh \
					chem/SMILES_grammar.hh \
					chem/SMILES_grammar.icc \
					chem/SMILESparser.hh \
					chem/SMILESwriter.hh \
					chem/SMILESwriter.icc \
					chem/GS_SMILES.hh \
					chem/GS_SMILES.icc \
					chem/Reaction.hh \
					chem/ReactionTransitionState.hh \
					chem/ReactionTransitionState.icc \
					chem/ReactionRateCalculation.hh \
					chem/MR_Reactions.hh \
					chem/MR_Reactions.icc \
					chem/MoleculeComponent.hh \
					chem/MC_MC_Node.hh \
					chem/MoleculeComponent_GML_grammar.hh \
					chem/MoleculeComponent_GML_grammar.icc \
					chem/MoleculeComponent_GMLparser.hh \
					chem/EnergyCalculationConstants.hh \
					chem/EnergyCalculation.hh \
					chem/MoleculeDecomposition.hh \
					chem/MoleculeDecomposition.icc \
					chem/MoleculeDecomposition_groups.icc \
					chem/MoleculeDecomposition_interactions.icc \
					chem/MoleculeDecomposition_singleAtoms.icc \
					chem/MoleculeDecomposition_smallMolecules.icc \
					chem/EC_MoleculeDecomposition.hh \
					chem/EC_MoleculeDecomposition.icc \
					chem/RRC_ArrheniusLaw.hh \
					$(gglchem_OB_headerfiles)

# header install directory via dummy lib
lib_gglchemdir = $(includedir)/ggl/chem


################################################################
# the GGL CHEMISTRY library objects
################################################################

libgglchem_a_SOURCES = \
					$(lib_gglchem_HEADERS) \
					chem/MoleculeUtil.cc \
					chem/SMILESparser.cc \
					chem/SMILESwriter.cc \
					chem/ChemRule.cc \
					chem/Reaction.cc \
					chem/ReactionTransitionState.cc \
					chem/MR_Reactions.cc \
					chem/GS_MolCheck.cc \
					chem/GS_chem.cc \
					chem/ReactionRateCalculation.cc \
					chem/AromaticityPerception.cc \
					chem/AP_NSPDK.cc \
					chem/AP_NSPDK_Model_OpenBabel_2013.cc \
					chem/AP_NSPDK_Model_Marvin_general_2013.cc \
					chem/AP_NSPDK_Model.cc \
					chem/EnergyCalculationConstants.cc \
					chem/MoleculeComponent.cc \
					chem/MoleculeComponent_GMLparser.cc \
					chem/MoleculeDecomposition.cc \
					chem/RRC_ArrheniusLaw.cc
					

################################################################

################################################################

SUBDIRS = 
