
#################################################################
# TEST CASES TO RUN : 
#################################################################


TESTINCLUDES =	\
				dataTargetGraph_1.icc \
				dataTargetGraph_2.icc \
				dataTargetGraph_3.icc \
				dataTargetGraph_4.icc \
				dataRule_1.icc \
				dataRule_2.icc \
				dataRule_3.icc \
				dataRule_4.icc \
				dataChemRule_1.icc \
				dataSGMgraphs_1.icc \
				dataSGMgraphs_2.icc \
				dataSGMgraphs_3.icc \
				utilPrintGraph_Interface.icc \
				utilPrintRule.icc

# set the test routines for OpenBabel depending classes
OB_FILES = \
				testggl_chem_SMILESwriterOB.cc \
				testggl_chem_MoleculeOB.cc 
		
if with_openbabel
TESTSOURCES_OB = $(OB_FILES)
else
TESTSOURCES_OB = 
endif

# set the list of all test routines to compile
TESTSOURCES =   \
				testsgm_Graph_boost.cc \
				testsgm_Graph_NodeLabelPrefix.cc \
				testsgm_Graph_boostV_p.cc \
				testsgm_SubGraph.cc \
				testsgm_MR_stream.cc \
				testsgm_MR_Counting.cc \
				testsgm_MR_Storing.cc \
				testsgm_SGM_vf2.cc \
				testsgm_GM_vf2.cc \
				testsgm_MC_Edge.cc \
				testsgm_MC_Node.cc \
				testggl_CE_Horton.cc \
				testggl_GS_stream.cc \
				testggl_Rule.cc \
				testggl_Rule_GML_grammar.cc \
				testggl_Rule_GMLparser.cc \
				testggl_Rule_GML_writer.cc \
				testggl_MR_ApplyRule.cc \
				testggl_MR_ApplyRule_NodeLabelPrefix.cc \
				testggl_GS_STL_pushAll.cc \
				testggl_GS_STL_pushAllP.cc \
				testggl_GS_STL_pushUnique.cc \
				testggl_GS_STL_pushUniqueP.cc \
				testggl_Graph_GML_writer.cc \
				testggl_Graph_GML_grammar.cc \
				testggl_Graph_GXL_writer.cc \
				testsgm_RP_Hanser96.cc \
				testsgm_GraphScaffold.cc \
				testggl_chem_MoleculeUtil.cc \
				testggl_chem_SMILES_grammar.cc \
				testggl_chem_SMILESparser.cc \
				testggl_chem_SMILESwriter.cc \
				testggl_chem_GS_SMILES.cc \
				testggl_chem_MR_Reactions.cc \
				testggl_chem_AP_NSPDK.cc \
				testggl_chem_MoleculeComponent_GML_grammar.cc \
				testggl_chem_MoleculeGroups.cc \
				testggl_chem_Molecule_Graph_noClass.cc \
				testggl_chem_MoleculeDecomposition.cc \
				testggl_chem_RRC_ArrheniusLaw.cc \
				$(TESTSOURCES_OB)

TESTSCRIPTS =	\
				testbin_toyChem.sh \
				testbin_molTool.sh \
				testbin_aromModelNSPDK.sh \
				testbin_gameOfLife.sh \
				testbin_sudoku.sh 
#				testbin_hanser.sh testbin_catalan.sh # currently disabled due to the BGL bug for version < 1.50.0


TESTDATA =		\
				dielsalder/da_rule.gml \
				dielsalder/da1.smi \
				dielsalder/da1.gml \
				dielsalder/da2.smi \
				dielsalder/da2.gml \
				dielsalder/da3.smi \
				toyChem/groups.gml \
				toyChem/lactat.smi \
				toyChem/NADH.smi \
				toyChem/lactat-dehydrogenase.gml \
				toyChem/purine.smi \
				game-of-life/gof.rules \
				game-of-life/gof.glider \
				sudoku/very-easy.input \
				sudoku/easy.input \
				sudoku/medium.input \
				sudoku/hard.input \
				catalan/2-step.input \
				catalan/3-step.input \
				hanser/single-loop.input \
				hanser/single-ring.input \
				hanser/two-fused-rings.input \
				hanser/three-fused-rings.input \
				hanser/four-fused-rings.input \
				hanser/cube.input


# set test specific compiler and linker options
TESTCXXFLAGS =	-I$(top_srcdir)/src \
				-I$(top_srcdir)/src/@LOCAL_VF2_FOLDER@ \
				-I$(top_srcdir)/src/@LOCAL_SGD_FOLDER@ \
				-I$(top_srcdir)/src/@LOCAL_NSPDK_FOLDER@ \
				@AM_CXXFLAGS@

TESTLDFLAGS  = 

TESTLIBS     = 	$(top_builddir)/src/ggl/libgglchem.a \
				$(top_builddir)/src/ggl/libggl.a \
				$(top_builddir)/src/sgm/libsgm_vf2.a \
				$(top_builddir)/src/sgm/libsgm_nspdk.a \
				$(top_builddir)/src/sgm/libsgm.a \
				$(top_builddir)/src/@LOCAL_NSPDK_FOLDER@/libnspdk.a \
				$(top_builddir)/src/@LOCAL_SGD_FOLDER@/libsgd.a \
				$(top_builddir)/src/@LOCAL_VF2_FOLDER@/libvf2.a 

EXTRA_DIST =	$(TESTSOURCES) $(TESTSOURCES:.cc=.verified) \
				$(OB_FILES) $(OB_FILES:.cc=.verified) \
				$(TESTSCRIPTS) $(TESTSCRIPTS:.sh=.verified) \
				$(TESTINCLUDES) $(TESTDATA)


#################################################################
# AUXILIARY VARIABLES FOR TEST BUILD AND EXECUTION
#################################################################

# final compiler and linker flag adds
AM_CXXFLAGS = $(TESTCXXFLAGS) @AM_CXXFLAGS@
AM_LDFLAGS = $(TESTLDFLAGS) @AM_LDFLAGS@
LIBS = $(TESTLIBS) @LIBS@

DIFFS=$(TESTSOURCES:.cc=.diff) $(TESTSCRIPTS:.sh=.diff)

 # build without warnings
TESTCXXCOMPILE = $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) -w
TESTCXXLINK = $(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -w \
	-o $@

RUNDIFF = $(DIFF) -w

#################################################################
# MAKE TARGETS FOR TEST BUILD AND EXECUTION
#################################################################

all:

tests: $(DIFFS)


.cc.diff:
	@name=$*; name=$${name/test/}; \
	$(PRINTF) " Building test %-38s" $$name;
	    @$(TESTCXXCOMPILE) -c -o $*.o $<
	    @$(TESTCXXLINK) $*.o $(LIBS) -o $*.test
	@$(PRINTF) " --> execution :";
	    @./$*.test > $*.lastout
	@name=$*; name=$${name/test/}; \
	if [ -e $*.verified ] ; then \
	    if $(RUNDIFF) $*.verified $*.lastout >$@ ; then \
	        $(PRINTF) " OK\n"; \
	    else \
	        $(PRINTF) " FAILED\n";\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(PRINTF) " Diff for $$name (verified '<', current '>')\n" ;\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(CAT) $@ ;\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(PRINTF) " [ If new output is correct: mv $*.lastout $*.verified ]\n" ;\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(PRINTF) "\n" ;\
	    fi ;\
	    $(RM) $@ ;\
	else \
	    $(PRINTF) " DONE\n"; \
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(PRINTF) " Output of $$name :\n" ;\
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(CAT) $*.lastout ;\
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(PRINTF) " WARNING: verified output of $$name is missing. Please generate\n" ;\
	    $(PRINTF) " a file $*.verified that contains the verfied output of $$name.\n"; \
	    $(PRINTF) " [ If output is correct: mv $*.lastout $*.verified ]\n" ;\
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(PRINTF) "\n" ;\
	fi

.sh.diff:
	@name=$*; name=$${name/test/}; \
	$(PRINTF) " Running  test %-38s" $$name;
	@$(PRINTF) " --> execution :";
	    @$(SHELL) ./$*.sh > $*.lastout
	@name=$*; name=$${name/test/}; \
	if [ -e $*.verified ] ; then \
	    if $(RUNDIFF) $*.verified $*.lastout >$@ ; then \
	        $(PRINTF) " OK\n"; \
	    else \
	        $(PRINTF) " FAILED\n";\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(PRINTF) " Diff for $$name (verified '<', current '>')\n" ;\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(CAT) $@ ;\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(PRINTF) " [ If new output is correct: mv $*.lastout $*.verified ]\n" ;\
	        $(PRINTF) "------------------------------------------------------\n" ;\
	        $(PRINTF) "\n" ;\
	    fi ;\
	    $(RM) $@ ;\
	else \
	    $(PRINTF) " DONE\n"; \
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(PRINTF) " Output of $$name :\n" ;\
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(CAT) $*.lastout ;\
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(PRINTF) " WARNING: verified output of $$name is missing. Please generate\n" ;\
	    $(PRINTF) " a file $*.verified that contains the verfied output of $$name.\n"; \
	    $(PRINTF) " [ If output is correct: mv $*.lastout $*.verified ]\n" ;\
	    $(PRINTF) "------------------------------------------------------\n" ;\
	    $(PRINTF) "\n" ;\
	fi

clean-local:
	$(RM) -rf $(TESTSOURCES:.cc=.test) $(TESTSOURCES:.cc=.o) \
	    $(TESTSOURCES:.cc=.lastout) $(TESTSOURCES:.cc=.diff) \
	    $(TESTSCRIPTS:.sh=.lastout) $(TESTSCRIPTS:.sh=.diff) \
	    $(TESTSOURCES:.cc=.exe) *.stackdump

#################################################################
