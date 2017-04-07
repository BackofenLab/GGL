#ifndef bool
#include <iostream>
#endif

#include <ggl/chem/SMILES_grammar.hh>
#include <ggl/chem/GS_SMILES.hh>
#include <ggl/chem/Reaction.hh>
#include <ggl/chem/MR_Reactions.hh>
#include <ggl/chem/Molecule.hh>
#include <ggl/chem/ChemRule.hh>
#include <ggl/chem/ChemRuleGraph.hh>
#include <ggl/Rule_GML_grammar.hh>
#include <ggl/Graph_GML_grammar.hh>
#include <ggl/MR_ApplyRule.hh>
#include <bin/toyChemUtil.hh>
#include <ggl/chem/AP_NSPDK.hh>
#include <ggl/Rule_GMLparser.hh>

extern "C" {
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
}
#ifdef bool
#undef bool
#include <iostream>
#endif


#define AROMATICITY_MODEL "Marvin:general:2013"


#define OUT_SMILES 2
#define OUT_RXNS 1

void
correctInputMolecules( const SMILES_container & targetSmiles, 
               SMILES_container & producedSmiles,
               ggl::chem::AromaticityPerception & aromaticity_withH);

void
correctMolecule(const ggl::chem::Molecule & m1,
               ggl::chem::Molecule & m2,
               ggl::chem::AromaticityPerception & aromaticity_withH);


SV* toyChem(SV* substrates, SV* ruleref, int outmode) {
  int i;

  if (!SvROK(substrates) || SvTYPE(SvRV(substrates)) != SVt_PVAV) {
    croak("substrates: expected ARRAY ref");
  }

  if (!SvROK(ruleref) || SvTYPE(SvRV(ruleref)) != SVt_PVAV) {
    croak("rules: expected ARRAY ref");
  }

  AV* subs_ary = (AV*) SvRV(substrates);
  AV* rules_ary = (AV*) SvRV(ruleref);
  
  SMILES_container parsedSmiles;
  SMILES_container substrateSmiles;
  SMILES_container strangeSmiles;

  ggl::chem::AP_NSPDK * aromaticity_withH = new ggl::chem::AP_NSPDK(AROMATICITY_MODEL);

  for (i=0; i<=av_len(subs_ary); i++) {
    SV** t_subs = av_fetch(subs_ary, i, 0);
    std::pair<ggl::chem::Molecule,int> result = 
      ggl::chem::SMILES_grammar::parseSMILES(SvPV_nolen(*t_subs));
      //GMLgrammar::parseGraph(SvPV_nolen(*t_subs));
    if (result.second != -1) {
      croak("Invalid input string: %s", SvPV_nolen(*t_subs));
    }

    size_t checkResult = ggl::chem::MoleculeUtil::isConsistent(result.first);
    if (checkResult != ggl::chem::MoleculeUtil::C_Consistent) {
      croak("Invalid molecule: %s", SvPV_nolen(*t_subs));
    }

    const std::string canSMILES =
      ggl::chem::SMILESwriter::getSMILES( result.first );
    parsedSmiles[canSMILES] = new ggl::chem::Molecule(result.first);
  }
  correctInputMolecules(parsedSmiles, substrateSmiles, *aromaticity_withH);

  std::vector<ggl::chem::ChemRule> rules;
  RulePatternMap rulePattern;
  for (i=0; i<=av_len(rules_ary); i++) {
    SV** t_rule = av_fetch(rules_ary, i, 0);
    std::string rule = SvPV_nolen(*t_rule);
    std::pair<ggl::chem::ChemRule,int> result =
      ggl::Rule_GMLparser::parseRule(rule);
      
    if (result.second != -1) {
      croak("Invalid rule string:\n%s", rule.c_str());
    }

    rules.push_back(result.first);
  }

  for (size_t r=0; r<rules.size(); r++) {
    ggl::chem::LeftSidePattern* pattern = new ggl::chem::LeftSidePattern(rules[r]);
    // store in the pattern list according to the component number
    size_t compNumber = pattern->getFirstOfEachComponent().size();
    rulePattern[compNumber].push_back( pattern );
  }


  SMILES_container toFill;
  ggl::chem::MR_Reactions::Reaction_Container producedReactions;


  applyRules(rulePattern, strangeSmiles, substrateSmiles, toFill, producedReactions, NULL, true, *aromaticity_withH);

  AV* ret = newAV();
  if (outmode == OUT_RXNS) {
    for ( ggl::chem::MR_Reactions::Reaction_Container::const_iterator
        r = producedReactions.begin();
        r != producedReactions.end(); r++) {

      HV* rxn = newHV();
      AV* rxn_substrates = newAV();
      AV* rxn_products = newAV();
      if ((*r).metabolites.size() == 0 || (*r).products.size() == 0) continue;
      ggl::chem::Reaction::Metabolite_Container::const_iterator s = (*r).metabolites.begin();
      ggl::chem::Reaction::Product_Container::const_iterator p = (*r).products.begin();
      for(; s != (*r).metabolites.end(); s++) 
        av_push(rxn_substrates, newSVpv((*s).c_str(), 0));
      for(; p != (*r).products.end(); p++)
        av_push(rxn_products, newSVpv((*p).c_str(), 0));

      hv_store(rxn, "substrates", 10, newRV((SV*)rxn_substrates), 0);
      hv_store(rxn, "products", 8, newRV((SV*)rxn_products), 0);
      av_push(ret, newRV((SV*)rxn));
    }
  }
  else {
    for (SMILES_container::const_iterator it=toFill.begin();
                                    it!=toFill.end() ; it++) {
      //std::cout << (*it).first << "\n";
      av_push(ret, newSVpv((*it).first.c_str(), 0));
    }
  }
  return newRV((SV*)ret);

}

SV* canonicalSmiles(SV* in_str, int noProtons) {

  using namespace ggl;
  using namespace ggl::chem;

  AP_NSPDK * aromaticity_withH = new AP_NSPDK(AROMATICITY_MODEL);
  std::pair<Molecule,int> result = 
    //GMLgrammar::parseGraph(SvPV_nolen(gml));
    SMILES_grammar::parseSMILES(SvPV_nolen(in_str));
  if (result.second != -1) {
    croak("Invalid input string: %s", SvPV_nolen(in_str));
  }

  Molecule corrected;
  correctMolecule(result.first, corrected, *aromaticity_withH);

  if (noProtons!=0) {
    MoleculeUtil::removeProtons(corrected);
  }

  const std::string canSMILES =
    SMILESwriter::getSMILES( corrected );

  return newSVpv(canSMILES.c_str(), 0);
}

SV* checkRules(SV* ruleref) {
  if (!SvROK(ruleref) || SvTYPE(SvRV(ruleref)) != SVt_PVAV) {
    croak("rules: expected ARRAY ref");
  }
  AV* rules_ary = (AV*) SvRV(ruleref);
  AV* ret = newAV();

  for(int i=0; i<=av_len(rules_ary); i++) {
    HV* h_res = newHV();
    std::string msg("");
    int res=1;
    SV** t_rule = av_fetch(rules_ary, i, 0);
    std::string rule = SvPV_nolen(*t_rule);
    std::pair<ggl::chem::ChemRule, int> result =
      ggl::Rule_GMLparser::parseRule(rule);

    if (result.second != -1) {
      res = 0;
      msg = "Invalid rule string";
    }
    size_t conStatus = result.first.isConsistent();
    if( conStatus != ggl::chem::ChemRule::C_Consistent ) {
      res = 0;
      std::ostringstream errorStream;
      // decode consistency status
      ggl::chem::MoleculeUtil::decodeConsistencyStatus(conStatus,errorStream);
      // get final error string
      msg = "consistency check failed for rule ID "+result.first.getID()+": "+errorStream.str();
    }
    hv_store(h_res, "res", 3, newSVnv(res), 0);
    hv_store(h_res, "msg", 3, newSVpv(msg.c_str(), 0), 0);
    av_push(ret, newRV((SV*)h_res));
  }
  return newRV((SV*)ret);
}


void correctInputMolecules( const SMILES_container & targetSmiles,
          SMILES_container & producedSmiles,
          ggl::chem::AromaticityPerception & aromaticity_withH ) {

  using namespace ggl;
  using namespace ggl::chem;

  Molecule m1;
  std::string SMILES;
  SMILESwriter writer;

  for(SMILES_container::const_iterator it=targetSmiles.begin(); it!=
        targetSmiles.end(); ++it) {
    
      // fill missing protons
    m1 = *(it->second);
    MoleculeUtil::fillProtons( m1 );
    // correct aromaticity
    aromaticity_withH.correctAromaticity( m1, false );
    // create SMILES from compressed molecule representation (protons
    // are removed if inferable)
    SMILES = SMILESwriter::getSMILES( m1, true );

      // store corrected molecule if not already present
    if (SMILES.size() > 0 && producedSmiles.find(SMILES) == producedSmiles.end()) {
      producedSmiles[SMILES] = new Molecule( m1 );
    }


  }
}

void correctMolecule(const ggl::chem::Molecule & hm, 
          ggl::chem::Molecule & m2,
          ggl::chem::AromaticityPerception & aromaticity_withH) {

  using namespace ggl;
  using namespace ggl::chem;

  Molecule m1 = hm;
  MoleculeUtil::fillProtons(m1);
  aromaticity_withH.correctAromaticity(m1,false);
  MoleculeUtil::compressHnodes(m1);
}



MODULE = Chemistry::GGL          	PACKAGE = Chemistry::GGL	

PROTOTYPES: DISABLE

SV *
toyChem(substrates, ruleref, outmode)
	SV *	substrates
	SV *	ruleref
	int	outmode
    
SV *
canonicalSmiles(gml, noProtons)
	SV *	gml
  int noProtons
    
SV * 
checkRules(ruleref) 
  SV *  ruleref

