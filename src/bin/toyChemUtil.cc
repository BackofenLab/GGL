
#include "toyChemUtil.hh"

#include "sgm/MR_restoreNodeLabel.hh"

#include "ggl/chem/ReactionRateCalculation.hh"
#include "ggl/Graph_GML_writer.hh"

#include <ggl/GS_stream.hh>
#include <ggl/Rule_GML_error.hh>

#include <ggl/chem/ChemRule.hh>
#include <ggl/chem/SMILESparser.hh>
#include <ggl/chem/GS_MolCheck.hh>
#include "ggl/chem/GS_SMILES.hh"
#include <ggl/chem/Molecule_Graph_noClass.hh>

#include <boost/lexical_cast.hpp>

#include <fstream>
#include <algorithm>




//////////////////////////////////////////////////////////////////////////

/*
 * Cuts a ':' separated list of files and stores all non-empty list elements
 * in the returned vector.
 * 
 * @param list the list of file names
 * @return the list of non-empty list elements
 */
std::vector< std::string >
splitToFileNames( const std::string & list ) {
	
	const char LIST_SEP = ':';
	
	size_t start = 0, end = 0;
	std::vector< std::string > files;
	while ( start != std::string::npos ) {
		end = list.find_first_of( LIST_SEP, start );
		if (end != std::string::npos) { // handle internal file name
			if (start != end) { // handle consecutive separators
				files.push_back( list.substr( start, end-start ) );
			}
			start = end + 1;
		} else {  // handle last file name
			if ( start < list.size() ) {
				files.push_back( list.substr( start ) );
			}
			start = end;
		}
	}
	
	return files;
}

//////////////////////////////////////////////////////////////////////////


#include <ggl/Rule_GMLparser.hh>

/*
 * Parses the given input stream for rules and adds each rule to the given
 * container.
 * 
 * @param in the input stream to read from
 * @param toFill the container to add the rules to
 * @param linesRead the number of lines already read from input (needed for
 *                  error reporting)
 *
 */
void
parseRules(	std::istream & in
			, std::vector<ggl::chem::ChemRule> & toFill
			, const size_t linesRead
			, const ggl::chem::GroupMap & groups
			)  throw(std::exception)
{
	const std::string whiteSpaces(" \t\n");
	std::string line;
	size_t lineNumber = linesRead;
	std::string ruleString = std::string("");
	  // parse whole input stream
	while (in.good()) {
		  // read next line
		std::getline( in, line );
		lineNumber++;

		  // skip comment lines
		if (	(!line.empty()) // no empty line
				&& line.at(0) == '#')  // leading comment sign
		{
			continue;
		} else {
			  // else append to rule string
			ruleString.append(line);
		}

		  // flag to mark end of stream
		bool endOfStream = in.eof();

		  // check if next rule header or end of file was reached
		size_t rulePos = line.find("rule",0);
		while (	rulePos != std::string::npos || endOfStream )
		{
			if (	endOfStream
					&& rulePos == std::string::npos
					)
			{
				if (ruleString.find("rule") != std::string::npos) {
					try {
						std::pair<ggl::Rule, int>
								ret = ggl::Rule_GMLparser::parseRule( ruleString );
						  // check for parsing error
						if (ret.second != -1) {
							std::ostringstream oss;
							oss	<<" Rule parsing error in rule specification :\n'"
								<<ruleString
								<<"'\n at rule string position " <<ret.second;
							throw ArgException(oss.str());
						}
						  // replace group label
						ggl::chem::ChemRule::CoreGraph newCore(ret.first.getCore());
						ggl::chem::ChemRule::insertGroups( newCore, groups );
						  // add to container
						toFill.push_back(ggl::chem::ChemRule(  newCore
								, ret.first.getID()
								, ret.first.getMatchConstraints()
								, ret.first.getCopyAndPasteOperations()));
					} catch (ggl::Rule_GML_error &ex) {
						  // catch GML parsing errors
						std::ostringstream oss;
						oss	<<" Rule parsing error in rule specification :\n'"
							<<ruleString
							<<"'\n\n\t" <<ex.what();
						throw ArgException( oss.str() );
					}
				}
				break;
			} else if (line.find("ruleID",rulePos) != rulePos) {
				  // calculate rulePos position in rule string
				size_t subEnd = ruleString.size() + rulePos - line.size();

				std::string curRuleString = ruleString.substr( 0, subEnd);
				if (curRuleString.find("rule") != std::string::npos) {
					try {
						std::pair<ggl::chem::ChemRule, int>
								ret = ggl::Rule_GMLparser::parseRule( curRuleString );
						  // check for parsing error
						if (ret.second != -1) {
							std::ostringstream oss;
							oss	<<" Rule parsing error in rule specification :\n'"
								<<curRuleString
								<<"'\n at rule string position " <<ret.second;
							throw ArgException(oss.str());
						}
						  // replace group label
						ggl::chem::ChemRule::CoreGraph newCore(ret.first.getCore());
						ggl::chem::ChemRule::insertGroups( newCore, groups );
						  // add to container
						toFill.push_back(ggl::chem::ChemRule(  newCore
								, ret.first.getID()
								, ret.first.getMatchConstraints()
								, ret.first.getCopyAndPasteOperations()));
					} catch (ggl::Rule_GML_error &ex) {
						  // catch GML parsing errors
						std::ostringstream oss;
						oss	<<" Rule parsing error in rule specification :\n'"
							<<curRuleString
							<<"'\n\n\t" <<ex.what();
						throw ArgException( oss.str() );
					}
				}
				  // shorten line and rule string : cut handled leading substr
				line = line.substr(rulePos);
				ruleString = ruleString.substr(subEnd);
				  // reset rulePos to start of shortened line
				rulePos = 0;
			}
			rulePos = line.find("rule",rulePos+1);
		}
	}

	// everything parsed .. hopefully .. ;)
}



//////////////////////////////////////////////////////////////////////////

/*
 * Parses the given input stream for rules in compacted GML and adds each rule
 * to the given container.
 * 
 * Rules are numbered consecutively.
 *
 * @param in the input stream to read from
 * @param toFill the container to add the rules to
 * @param linesRead the number of lines already read from input (needed for 
 *                  error reporting)
 * 
 */
void
parseCompactRules(	std::istream & in
					, std::vector<ggl::chem::ChemRule> & toFill
					, const size_t linesRead
					, const ggl::chem::GroupMap & groups
					)  throw(std::exception)
{
	const std::string whiteSpaces(" \t\n");
	std::string line;
	size_t lineNumber = linesRead;
	std::string ruleString = std::string("");
	size_t ruleID = 1;
	  // parse whole input stream
	while (in.good()) {
		  // read next line
		std::getline( in, line );
		lineNumber++;
		
		  // skip comment lines
		if (	(!line.empty()) // no empty line
				&& line.at(0) == '#')  // leading comment sign 
		{
			continue;
		} else {
			  // else append to rule string
			ruleString.append(line);
		}
		
		  // flag to mark end of stream
		bool endOfStream = in.eof();
		
		  // check if next rule header or end of file was reached
		size_t rulePos = line.find("graph",0);
		while (	rulePos != std::string::npos || endOfStream ) 
		{
			if (	endOfStream
					&& rulePos == std::string::npos
					)
			{
				if (ruleString.find("graph") != std::string::npos) {
					try {
						std::pair<ggl::Rule, int>
								ret = ggl::Rule_GMLparser::parseCompactRule( ruleString );
						  // check for parsing error
						if (ret.second != -1) {
							std::ostringstream oss;
							oss	<<" Rule parsing error in rule specification :\n'"
								<<ruleString
								<<"'\n at rule string position " <<ret.second;
							throw ArgException(oss.str());
						}
						  // set ID
						ret.first.setID("rule "+boost::lexical_cast<std::string>(ruleID++));
						  // replace group label
						ggl::chem::ChemRule::CoreGraph newCore(ret.first.getCore());
						ggl::chem::ChemRule::insertGroups( newCore, groups );
						  // add to container
						toFill.push_back(ggl::chem::ChemRule(  newCore
								, ret.first.getID()
								, ret.first.getMatchConstraints()
								, ret.first.getCopyAndPasteOperations()));
					} catch (ggl::Rule_GML_error &ex) {
						  // catch GML parsing errors
						std::ostringstream oss;
						oss	<<" Rule parsing error in rule specification :\n'"
							<<ruleString
							<<"'\n\n\t" <<ex.what();
						throw ArgException( oss.str() );
					}
				}
				break;
			} else {
				  // calculate rulePos position in rule string
				size_t subEnd = ruleString.size() + rulePos - line.size();
				
				std::string curRuleString = ruleString.substr( 0, subEnd);
				if (curRuleString.find("graph") != std::string::npos) {
					try {
						std::pair<ggl::chem::ChemRule, int>
								ret = ggl::Rule_GMLparser::parseCompactRule( curRuleString );
						  // check for parsing error
						if (ret.second != -1) {
							std::ostringstream oss;
							oss	<<" Rule parsing error in rule specification :\n'"
								<<curRuleString
								<<"'\n at rule string position " <<ret.second;
							throw ArgException(oss.str());
						}
						  // set ID
						ret.first.setID("rule "+boost::lexical_cast<std::string>(ruleID++));
						  // replace group label
						ggl::chem::ChemRule::CoreGraph newCore(ret.first.getCore());
						ggl::chem::ChemRule::insertGroups( newCore, groups );
						  // add to container
						toFill.push_back(ggl::chem::ChemRule(  newCore
								, ret.first.getID()
								, ret.first.getMatchConstraints()
								, ret.first.getCopyAndPasteOperations()));
					} catch (ggl::Rule_GML_error &ex) {
						  // catch GML parsing errors
						std::ostringstream oss;
						oss	<<" Rule parsing error in rule specification :\n'"
							<<curRuleString
							<<"'\n\n\t" <<ex.what();
						throw ArgException( oss.str() );
					}
				}
				  // shorten line and rule string : cut handled leading substr
				line = line.substr(rulePos);
				ruleString = ruleString.substr(subEnd);
				  // reset rulePos to start of shortened line
				rulePos = 0;
			}
			rulePos = line.find("graph",rulePos+1);
		}
	}
	
	// everything parsed .. hopefully .. ;)
}


//////////////////////////////////////////////////////////////////////////


void
parseRules(	const std::string & inSource
			, std::vector<ggl::chem::ChemRule> & toFill
			, const ggl::chem::GroupMap & groups
			, const bool compactGML ) throw(std::exception)
{
	if (inSource.compare("STDIN") == 0) {
		size_t linesRead = 0;
		if (compactGML) {
			parseCompactRules( std::cin, toFill, linesRead, groups );
		} else {
			parseRules( std::cin, toFill, linesRead, groups );
		}
		return;
	}
	
	  // split inSource into file name list if any
	std::vector< std::string > files = splitToFileNames( inSource );
	
	  // read all files
	size_t i=0;
	for (	std::vector< std::string >::const_iterator file=files.begin(); 
			file != files.end(); ++file)
	{
		i++;
		  // ensure there is no standard stream within the list
		if (file->compare("STDIN") == 0) {
			std::ostringstream oss;
			oss	<<"RULES input list : the " <<i
				<<". entry is 'STDIN', but it has to be a file list!";
			throw ArgException(oss.str());
		}
		 // open current file
		std::ifstream inFile( file->c_str(), std::ifstream::in );
		if (!inFile.is_open()) {
			std::ostringstream oss;
			oss	<<"cannot open " <<i <<". RULES input file '" <<(*file) <<"'";
			throw ArgException(oss.str());
		}
		 // parse current file
		size_t linesRead = 0;
		if (compactGML) {
			parseCompactRules( inFile, toFill, linesRead, groups );
		} else {
			parseRules( inFile, toFill, linesRead, groups );
		}
		 // close file
		inFile.close();
	}
}



//////////////////////////////////////////////////////////////////////////

#include <ggl/chem/SMILESparser.hh>
#include <ggl/chem/SMILESwriter.hh>

//////////////////////////////////////////////////////////////////////////

/*
 * Parses SMILES strings from stream. Each line should contain only ONE SMILES
 * string. Leading and tailing whitespaces are ignored.
 * 
 * @param in the stream to read from
 * @param toFill the inserter to add the found SMILES to
 * 
 */
size_t
parseSMILES(	std::istream & in
				, SMILES_container & toFill
				, const size_t linesRead
				, const ggl::chem::GroupMap & groups
				, const size_t setNextAtomClass
				) throw(std::exception)
{
	// get next atom class that is non-redundant
	size_t nextClassLabel = setNextAtomClass;

	std::string line;
	std::string smiles;
	size_t lineNumber = linesRead;
	const std::string WHITESPACES = std::string(" \t\n"); 
	  // parse whole input stream
	while (in.good()) {
		  // read next line
		std::getline( in, line );
		lineNumber++;
		
		  // check if line is not empty
		if (line.size() > 0) {
			size_t start = line.find_first_not_of(WHITESPACES);
			  // check if line is not only filled with blanks
			if (start != std::string::npos ) {
				size_t end = line.find_last_not_of(WHITESPACES);
				  // get SMILES string 
				smiles = line.substr( start, end - start + 1 );
				  // check if WHITESPACE within SMILES
				if (smiles.find(WHITESPACES) != std::string::npos) {
					std::ostringstream oss;
					oss	<<"parsing error in line " <<(linesRead+1)
						<<" in SMILES string '"
						<<smiles
						<<"' contains whitespaces";
					throw ArgException(oss.str());
				}
				  // check if '.' within SMILES, currently split not supported
				if (smiles.find('.') != std::string::npos) {
					std::ostringstream oss;
					oss	<<"parsing error in line " <<(linesRead+1)
						<<" in SMILES string '"
						<<smiles
						<<"' contains '.' connector, please split connected SMILES into single lines";
					throw ArgException(oss.str());
				}
				  // parse SMILES to graph
				std::pair<ggl::chem::Molecule,int> result 
					= ggl::chem::SMILESparser::parseSMILES(smiles,groups);
				  // check parsing result
				if (result.second != -1) {
					std::ostringstream oss;
					oss	<<"parsing error in line " <<(linesRead+1)
						<<" in SMILES string '"
						<<smiles
						<<"' at position " <<result.second;
					throw ArgException(oss.str());
				}

				// add atom class labels if needed
				if (nextClassLabel > 0) {
					// set class labels and overwrite existing labels
					nextClassLabel = setAtomClass( result.first, nextClassLabel, true );
				}

				try {
					  // derive canonical SMILES with own writer
					const std::string canSMILES =
						ggl::chem::SMILESwriter::getSMILES( result.first );

					  // store only if not already known
					if ( toFill.find(canSMILES) == toFill.end() ) {
						  // store parsed Molecule graph
						toFill[canSMILES] = new ggl::chem::Molecule(result.first);
					}
				} catch (std::runtime_error & e) {
					std::cerr <<"\n ERROR during SMILES parsing of '" <<smiles <<"' : " <<e.what() <<std::endl;
				}
			}
		}
	}

	return nextClassLabel;
}



//////////////////////////////////////////////////////////////////////////


size_t
parseSMILES(	const std::string & inSource
				, SMILES_container & toFill
				, const ggl::chem::GroupMap & groups
				, const size_t setNextAtomClass
				) throw(std::exception)
{
	if (inSource.compare("STDIN") == 0) {
		size_t linesRead = 0;
		return parseSMILES( std::cin, toFill, linesRead, groups, setNextAtomClass );
	}
	
	  // split inSource into file name list if any
	std::vector< std::string > files = splitToFileNames( inSource );
	
	  // read all files
	size_t i=0;
	size_t nextAtomClass = setNextAtomClass;
	for (	std::vector< std::string >::const_iterator file=files.begin(); 
			file != files.end(); ++file)
	{
		i++;
		  // ensure there is no standard stream within the list
		if (file->compare("STDIN") == 0) {
			std::ostringstream oss;
			oss	<<"SMILES input list : the " <<i
				<<". entry is 'STDIN', but it has to be a file list!";
			throw ArgException(oss.str());
		}
		 // open current file
		std::ifstream inFile( file->c_str(), std::ifstream::in );
		if (!inFile.is_open()) {
			std::ostringstream oss;
			oss	<<"cannot open " <<i <<". SMILES input file '" <<(*file) <<"'";
			throw ArgException(oss.str());
		}
		 // parse current file
		size_t linesRead = 0;
		nextAtomClass = parseSMILES( inFile, toFill, linesRead, groups, nextAtomClass );
		 // close file
		inFile.close();
	}

	return nextAtomClass;
}



//////////////////////////////////////////////////////////////////////////

#include <sgm/SubGraph.hh>
#include <ggl/chem/SMILESwriter.hh>
#include <ggl/Graph_GMLparser.hh>

////////////////////////////////////////////////////////////////////////////

size_t
parseMolGML(	std::istream & in
				, SMILES_container & toFill
				, const size_t linesRead
				, const ggl::chem::GroupMap & groups
				, const size_t setNextAtomClass
				, const bool pruneProtons
				, const bool reportToStream
				, std::ostream *error
				, const bool correctNodeByBonds
				)  throw(std::exception)
{
	using namespace ggl;
	using namespace ggl::chem;

	// get next atom class that is non-redundant
	size_t nextClassLabel = setNextAtomClass;

	typedef 
		ggl::Graph_GMLparser GraphParser;
	
	std::string line;
	size_t lineNumber = linesRead;
	std::string graphString = std::string("");
	std::string smiles = std::string("");
	  // parse whole input stream
	while (in.good()) {
		  // read next line
		std::getline( in, line );
		lineNumber++;
		
		  // check if next rule header or end of file was reached
		if (line.find("graph") != std::string::npos || in.eof()) {

			do  {
				  // append rest to previous graphString
				const size_t cutPos = line.find("graph");
				graphString.append( line.substr(0, cutPos) );
				line.erase( 0, cutPos );

				bool noParsingError = true;

				  // convert graph string to molecule object
				if (graphString.size() > 0) {
					std::pair<ggl::chem::Molecule, int>
							ret = GraphParser::parseGraph( graphString );
					  // check for parsing error
					if (ret.second != -1) {
						noParsingError = false;
						std::ostringstream oss;
						oss	<<" Graph parsing error in graph specification BEFORE line "
							<<line
							<<" at graph string position " <<ret.second;
						if (error != NULL){
							*error <<oss.str();
						} else {
							throw ArgException(oss.str());
						}
					}
					  // replace molecule component placeholders with
					  // according subgraphs from 'groups'
					MoleculeUtil::insertGroups(ret.first, groups);

					// add atom class labels if needed
					if (nextClassLabel > 0) {
						// set class labels and overwrite existing labels
						nextClassLabel = setAtomClass( ret.first, nextClassLabel, true );
					}

					  // check if final processing possible
					if (noParsingError) {

						  // check for individual connected component
						ggl::chem::Molecule_Graph parsedMols(ret.first);
						sgm::Graph_Interface::CompLabel compIDs;
						const size_t numOfMols = sgm::Graph_Interface::connectedComponents( parsedMols, compIDs );

						for (size_t curMolID = 0; curMolID < numOfMols; ++curMolID ) {

							  // get vector of all node ids for this component
							sgm::SubGraph::NodeList curCompIDs(compIDs.size());
							size_t curCompSize = 0;
							for (size_t x=0; x<compIDs.size(); ++x) {
								  // check if this node if part of the component
								if (compIDs.at(x) == curMolID) {
									  // add to component's node list
									curCompIDs[curCompSize] = x;
									++curCompSize;
								}
							}
							  // shrink container to final size
							curCompIDs.resize(curCompSize);

							  // get a molecule for the current component
							ggl::chem::Molecule g2;
							ggl::chem::MoleculeUtil::copy(sgm::SubGraph(parsedMols,curCompIDs),g2);

							  // get proton shrinked version of the molecule
	//						if (pruneProtons) {
	//							ggl::chem::MoleculeUtil::removeProtons(g2);
	//						} else {
	//							ggl::chem::MoleculeUtil::compressHnodes(g2);
	//						}
							  // correct aromaticity of node labels if needed
							if (correctNodeByBonds) {
								boost::property_map< ggl::chem::Molecule, ggl::PropNodeLabel >::type nodeLabel = get( ggl::PropNodeLabel(), g2);
								boost::property_map< ggl::chem::Molecule, ggl::PropEdgeLabel >::type edgeLabel = get( ggl::PropEdgeLabel(), g2);
								Molecule::vertex_iterator     vi, vi_end;
								Molecule::out_edge_iterator  e, e_end;

								using namespace boost;
								using namespace ggl::chem;

								for(boost::tie(vi, vi_end)=boost::vertices(g2); vi!=vi_end; ++vi) {
									  // check if node is currently aromatic
									bool nodeAromatic = MoleculeUtil::getAtomData( nodeLabel[*vi] )->isAromatic != 0;
									  // count aromatic edges
									size_t aromaticEdges = 0;
									for (boost::tie(e, e_end) = boost::out_edges(*vi, g2); e != e_end; ++e){
										  // check if aromatic edge
										if ( MoleculeUtil::getBondData( edgeLabel[*e] )->isAromatic != 0 ) {
											aromaticEdges++;
										}
									}
									  // error if only one aromatic edge
									if (aromaticEdges == 1) {
										noParsingError = false;
										if (error != NULL) {
											*error <<"\n ERROR : node with only one aromatic bond found!\n";
										} else {
											throw ArgException("ERROR : node with only one aromatic bond found!");
										}
									}
									  // do relabeling if needed
									if ( (aromaticEdges > 1 && !nodeAromatic)
										|| (aromaticEdges == 0 && nodeAromatic))
									{
										const std::string label = nodeLabel[*vi];
										nodeLabel[*vi] = MoleculeUtil::getComplexAtomLabel(
													*MoleculeUtil::getAromaticPendant( MoleculeUtil::getAtom( label ) )
													, MoleculeUtil::getProtons(label)
													, MoleculeUtil::getCharge(label)
													, MoleculeUtil::getClass(label)
													, false
												);
									}
								}
							}
							  // report conversion result
							if (reportToStream)  {
								ggl::GS_stream out(std::cout);

								out.add( ret.first );
								out.add( g2 );
							}

							try {
								  // create SMILES string for the molecule
								smiles = ggl::chem::SMILESwriter::getSMILES(g2,pruneProtons);

								  // store if not already known
								if (toFill.find(smiles) == toFill.end()) {
									  // store parsed Molecule graph
									toFill[smiles] = new ggl::chem::Molecule(ret.first);
								}
							} catch (std::runtime_error & e) {
								std::cerr <<"\n ERROR during GML parsing of '" <<graphString <<"' : " <<e.what() <<std::endl;
							}

						}  // for all connected components, ie. molecules in GML
					} // no parsing error
				}
				  // clear processed graph and initialize with rest of this line
				graphString = line;
			} while( line.find("graph",1) != std::string::npos );
		} else {
			  // append line content to graph string
			graphString.append(line);
		}
	}

	return nextClassLabel;
}

//////////////////////////////////////////////////////////////////////////


size_t
parseMolGML(	const std::string & inSource
				, SMILES_container & toFill
				, const ggl::chem::GroupMap & groups
				, const size_t setNextAtomClass
				) throw(std::exception)
{
	if (inSource.compare("STDIN") == 0) {
		size_t linesRead = 0;
		return parseMolGML( std::cin, toFill, linesRead, groups, setNextAtomClass );
	}
	
	  // split inSource into file name list if any
	std::vector< std::string > files = splitToFileNames( inSource );
	
	  // read all files
	size_t i=0;
	size_t nextAtomClass = setNextAtomClass;
	for (	std::vector< std::string >::const_iterator file=files.begin(); 
			file != files.end(); ++file)
	{
		i++;
		  // ensure there is no standard stream within the list
		if (file->compare("STDIN") == 0) {
			std::ostringstream oss;
			oss	<<"GML MOLECULES input list : the " <<i
				<<". entry is 'STDIN', but it has to be a file list!";
			throw ArgException(oss.str());
		}
		 // open current file
		std::ifstream inFile( file->c_str(), std::ifstream::in );
		if (!inFile.is_open()) {
			std::ostringstream oss;
			oss	<<"cannot open " <<i <<". GML MOLECULES input file '" 
				<<(*file) <<"'";
			throw ArgException(oss.str());
		}
		 // parse current file
		size_t linesRead = 0;
		nextAtomClass = parseMolGML( inFile, toFill, linesRead, groups, nextAtomClass );
		 // close file
		inFile.close();
	}
	return nextAtomClass;
}



//////////////////////////////////////////////////////////////////////////

void
parseGroups(	const std::string & inSource
				, ggl::chem::GroupMap & toFill ) throw(std::exception)
{
	if (inSource.compare("STDIN") == 0) {
		parseGroups( std::cin, toFill );
		return;
	}

	  // split inSource into file name list if any
	std::vector< std::string > files = splitToFileNames( inSource );

	  // read all files
	size_t i=0;
	for (	std::vector< std::string >::const_iterator file=files.begin();
			file != files.end(); ++file)
	{
		i++;
		  // ensure there is no standard stream within the list
		if (file->compare("STDIN") == 0) {
			std::ostringstream oss;
			oss	<<"SMILES input list : the " <<i
				<<". entry is 'STDIN', but it has to be a file list!";
			throw ArgException(oss.str());
		}
		 // open current file
		std::ifstream inFile( file->c_str(), std::ifstream::in );
		if (!inFile.is_open()) {
			std::ostringstream oss;
			oss	<<"cannot open " <<i <<". SMILES input file '" <<(*file) <<"'";
			throw ArgException(oss.str());
		}
		 // parse current file
		parseGroups( inFile, toFill );
		 // close file
		inFile.close();
	}
}



//////////////////////////////////////////////////////////////////////////

#include <ggl/chem/MoleculeComponent_GMLparser.hh>

/*
 * Parses SMILES strings from stream. Each line should contain only ONE SMILES
 * string. Leading and tailing whitespaces are ignored.
 *
 * @param in the stream to read from
 * @param toFill the inserter to add the found SMILES to
 *
 */
void
parseGroups(	std::istream & in
				, ggl::chem::GroupMap & toFill ) throw(std::exception)
{
	using namespace ggl;
	using namespace ggl::chem;

	typedef
		MoleculeComponent_GMLparser GroupParser;

	std::string line;
	std::string groupString = std::string("");
	const std::string * abbreviation = NULL;
	  // parse whole input stream
	while (in.good()) {
		  // read next line
		std::getline( in, line );

		  // check if next rule header or end of file was reached
		if (line.find("molcomp") != std::string::npos || in.eof()) {

			do  {
				  // append rest to previous groupString
				const size_t cutPos = line.find("molcomp");
				groupString.append( line.substr(0, cutPos) );
				line.erase( 0, cutPos );

				  // convert graph string to molecule object
				if (groupString.size() > 0) {
					  // do parsing
					std::pair<ggl::chem::MoleculeComponent, int>
							ret = GroupParser::parseGML( groupString );
					  // check for parsing error
					if (ret.second != -1) {
						std::ostringstream oss;
						oss	<<" Molecule group parsing error in graph specification BEFORE line "
							<<line
							<<" at graph string position " <<ret.second;
						throw ArgException(oss.str());
					}
					  // access to description
					abbreviation = &(ret.first.description);
					  // check if abbreviation is conform
					if (abbreviation->size()<2 || *(abbreviation->begin())!='{' || *(abbreviation->rbegin()) != '}') {
						std::ostringstream oss;
						oss	<<" Molecule Group parsing error : description has to be enclosed in brackets '{..}' : "
							<<"currently'" <<*abbreviation<<"'";
						throw ArgException(oss.str());
					}
					  // check if abbreviation unknown so far
					if (toFill.find(*abbreviation) != toFill.end()) {
						std::ostringstream oss;
						oss	<<" Molecule Group parsing error : description '"
								<<*abbreviation<<"' is already known or not unique";
						throw ArgException(oss.str());
					}
					  // check if only one compID given
					if (ret.first.compIDs.size() != 1) {
						std::ostringstream oss;
						oss	<<" Molecule Group parsing error : only one compID node id allowed : molcomp '"
								<<*abbreviation<<"'";
						throw ArgException(oss.str());
					}
					  // check if all nodes show a standard atom label
					boost::property_map<ggl::chem::MoleculeComponent::PatternGraph, ggl::PropNodeLabel>::type compNodeLabel
						= boost::get( ggl::PropNodeLabel(), ret.first.pattern);
					MoleculeComponent::PatternGraph::vertex_iterator  vi, vi_end;
					for(boost::tie(vi, vi_end)=boost::vertices(ret.first.pattern); vi!=vi_end; ++vi) {
						if (!ggl::chem::MoleculeUtil::isValidAtomLabel(compNodeLabel[*vi])) {
							std::ostringstream oss;
							oss	<<" Molecule Group parsing error : at least one node with invalid atom label : molcomp '"
									<<*abbreviation<<"'";
							throw ArgException(oss.str());
						}
					}
					  // store parsed molecule group
					toFill[*abbreviation] = ret.first;
				}
				  // clear processed graph and initialize with rest of this line
				groupString = line;
			} while( line.find("molcomp",1) != std::string::npos );
		} else {
			  // append line content to graph string
			groupString.append(line);
		}
	}
}


//////////////////////////////////////////////////////////////////////////

void
printRules( std::ostream& out, std::vector<ggl::chem::ChemRule> & rules )
{
	for (size_t i=0; i<rules.size(); ++i ) {
		printRule(out, rules[i]);
		out <<"\n";
	}

}

//////////////////////////////////////////////////////////////////////////

void
printRule( std::ostream& out, ggl::chem::ChemRule & rule )
{
	out <<"\n RULE : " <<rule.getID() <<"\n" <<std::endl; 
	{
		out <<"\n == LEFT_SIDE_PATTERN ==\n" <<std::endl;
		ggl::chem::LeftSidePattern g(rule);
		for (size_t i=0; i<g.getNodeNumber(); ++i ) {
			out <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
			for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i),
					eEnd = g.getOutEdgesEnd(i); e != eEnd; ++e)
			{
				out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
			}
			out <<" |\n";
		}
	}
	printConstraints(out, rule);
	{
		out <<"\n == RIGHT_SIDE_PATTERN ==\n" <<std::endl;
		ggl::RightSidePattern g(rule);
		for (size_t i=0; i<g.getNodeNumber(); ++i ) {
			out <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
			for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i),
					eEnd = g.getOutEdgesEnd(i); e != eEnd; ++e)
			{
				out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
			}
			out <<" |\n";
		}
	}
	printCopyAndPaste(out,rule);
	out <<"\n";
}

//////////////////////////////////////////////////////////////////////////


void printCopyAndPaste( std::ostream & out, const ggl::chem::ChemRule & rule )
{
	out <<"\n number of nodes with copy-and-paste operations = " <<rule.getCopyAndPasteOperations().size() <<"\n"<<std::endl;
	for ( ggl::Rule::CopyAndPasteOperations::const_iterator n=rule.getCopyAndPasteOperations().begin(); n!=rule.getCopyAndPasteOperations().end(); ++n) {
		for (ggl::Rule::CopyAndPasteOperations::mapped_type::const_iterator cnp=n->second.begin(); cnp!=n->second.end(); ++cnp) {
			out <<" copy-and-paste : source("<<cnp->source <<") pasteID("<<cnp->pasteID<<") target(";
			if (cnp->target != (size_t)INT_MAX)
				out <<cnp->target;
			out <<") edgeLabels (";
			for (ggl::Rule::RuleCnP::LabelSet::const_iterator l=cnp->edgeLabels.begin(); l!=cnp->edgeLabels.end(); ++l)
				out <<" '"<<*l<<"'";
			out <<")"<<std::endl;
		}
	}
}


//////////////////////////////////////////////////////////////////////////
void printConstraints( std::ostream & out, const ggl::chem::ChemRule & rule )
{
	out <<"\n number of constraints = " <<rule.getLeftSide().constraints.size() <<std::endl;

	typedef std::vector< sgm::Pattern_Interface::Match_Constraint* > CV;
	const CV & cv = rule.getLeftSide().constraints;
	for (CV::const_iterator c = cv.begin(); c != cv.end(); ++c ) {
		const sgm::MC_NodeLabel * cur = dynamic_cast< const sgm::MC_NodeLabel * >(*c);
		{
		if (cur != NULL) {
			out <<" + MC_NodeLabel "<<std::endl;
			out <<"   id " <<cur->constrainedNodeID <<std::endl;
			out <<"   op " <<(cur->compareType==sgm::MC_NodeLabel::ALLOWED?"=":"!") <<std::endl;
			out <<"   nodeLabels = " ;
			for (sgm::MC_NodeLabel::LabelSet::const_iterator l=cur->nodeLabels.begin(); l!=cur->nodeLabels.end(); ++l)
				out <<" '"<<*l<<"'" ;
			out <<std::endl;
			continue;
		}}{
		const sgm::MC_EdgeLabel * cur = dynamic_cast< const sgm::MC_EdgeLabel * >(*c);
		if (cur != NULL) {
			out <<" + MC_EdgeLabel "<<std::endl;
			out <<"   from " <<cur->constrainedFromID <<std::endl;
			out <<"   to   " <<cur->constrainedToID <<std::endl;
			out <<"   op " <<(cur->compareType==sgm::MC_EdgeLabel::ALLOWED?"=":"!") <<std::endl;
			out <<"   edgeLabels = " ;
			for (sgm::MC_EdgeLabel::LabelSet::const_iterator l=cur->edgeLabels.begin(); l!=cur->edgeLabels.end(); ++l)
				out <<" '"<<*l<<"'" ;
			out <<std::endl;
			continue;
		}}{
		const sgm::MC_NoEdge * cur = dynamic_cast< const sgm::MC_NoEdge * >(*c);
		if (cur != NULL) {
			out <<" + MC_NoEdge "<<std::endl;
			out <<"   from " <<cur->constrainedFromID <<std::endl;
			out <<"   to   " <<cur->constrainedToID <<std::endl;
			continue;
		}}{
		const sgm::MC_NodeAdjacency * cur = dynamic_cast< const sgm::MC_NodeAdjacency * >(*c);
		if (cur != NULL) {
			out <<" + MC_NodeAdjacency "<<std::endl;
			out <<"   id " <<cur->constrainedNodeID <<std::endl;
			out <<"   op " ;
			switch(cur->op) {
			case sgm::MC_NodeAdjacency::MC_EQ : out <<"="; break;
			case sgm::MC_NodeAdjacency::MC_G : out <<">"; break;
			case sgm::MC_NodeAdjacency::MC_GQ : out <<">="; break;
			case sgm::MC_NodeAdjacency::MC_L : out <<"<"; break;
			case sgm::MC_NodeAdjacency::MC_LQ : out <<"<="; break;
			default : out <<"unknown operator type";
			}
			out <<std::endl;
			out <<"   count " <<cur->count <<std::endl;
			out <<"   nodeLabels = " ;
			for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=cur->nodeLabels.begin(); l!=cur->nodeLabels.end(); ++l)
				out <<" '"<<*l<<"'" ;
			out <<std::endl;
			out <<"   edgeLabels = " ;
			for (sgm::MC_NodeAdjacency::LabelSet::const_iterator l=cur->edgeLabels.begin(); l!=cur->edgeLabels.end(); ++l)
				out <<" '"<<*l<<"'" ;
			out <<std::endl;

			continue;
		}}
		out <<" !!! unknown constraint type !!! " <<std::endl;
	}

}

//////////////////////////////////////////////////////////////////////////

void
printSMILES( std::ostream& out, SMILES_container & smiles )
{
	CompareString compare;
	  // create sorted SMILES list
	std::vector<std::string> toPrint(smiles.size(),"");
	size_t p=0;
	for (SMILES_container::const_iterator it = smiles.begin(); it != smiles.end(); ++it,++p ) {
		  // get proton compressed molecule
		try {
			toPrint[p] = ggl::chem::SMILESwriter::getSMILES(*(it->second), true);
		} catch (std::exception & ex) {
			std::cerr <<"ERROR in printSMILES for molecule '" <<it->first <<"' : " <<ex.what() <<std::endl;
		}
	}
	  // print sorted SMILES list
	std::sort(toPrint.begin(),toPrint.end(), compare);
	for (std::vector<std::string>::const_iterator s=toPrint.begin(); s!= toPrint.end(); ++s) {
		out <<" " <<*s <<"\n";
	}
}


//////////////////////////////////////////////////////////////////////////

void
printGraph( std::ostream& out, const sgm::Graph_Interface& g ) 
{
	for (size_t i=0; i<g.getNodeNumber(); ++i ) {
		out <<std::setw(6) <<i <<" ("<<g.getNodeLabel(i) <<")  --> ";
		for (sgm::Graph_Interface::OutEdge_iterator e = g.getOutEdgesBegin(i),
				eEnd = g.getOutEdgesEnd(i); e != eEnd; ++e)
		{
			out <<" | " <<e->getToIndex() <<" (" <<e->getEdgeLabel() <<")";
		}
		out <<" |\n";
	}

}

//////////////////////////////////////////////////////////////////////////

void
printGroups( std::ostream& out, const ggl::chem::GroupMap & groups )
{
	for (ggl::chem::GroupMap::const_iterator group = groups.begin();
			group != groups.end(); ++group)
	{
		  // print identifier
		out <<"group '" <<group->first <<"'\n";
		  // print graph
		ggl::Graph_GML_writer::write(out, group->second.pattern, true);
		out <<"\n";
	}

}

//////////////////////////////////////////////////////////////////////////

#include <sgm/SubGraph.hh>
#include <sgm/Graph_boostV_p.hh>
#include <sgm/MR_SymmBreak.hh>
#include <sgm/MR_Counting.hh>

//////////////////////////////////////////////////////////////////////////

bool 
singleRuleApplication(	sgm::SubGraphMatching& sgm
						, const ggl::chem::LeftSidePattern& rulePattern
						, const SMILES_container & allTargets
						, sgm::Match_Reporter& mrApplyRule
						, const sgm::Pattern_Automorphism* ruleSymmetry
						, const bool ignoreAtomClassLabel
						, const bool noRedundantMolecules
						)
{
	  // target graph holding compNbr different targets
	std::vector< const ggl::chem::Molecule* > 
			curTargets(rulePattern.getFirstOfEachComponent().size(),NULL);

	  // forward job to recursive version
	return singleRuleApplicationRec(	sgm
										, rulePattern
										, 0
										, allTargets
										, curTargets
										, mrApplyRule
										, ruleSymmetry
										, ignoreAtomClassLabel
										, noRedundantMolecules
										);
}


//////////////////////////////////////////////////////////////////////////

bool 
singleRuleApplicationRec(	sgm::SubGraphMatching& sgm
							, const ggl::chem::LeftSidePattern& rulePattern
							, const size_t ruleComponent
							, const SMILES_container & allTargets
							, std::vector< const ggl::chem::Molecule* >& curTargets
							, sgm::Match_Reporter& mrApplyRule
							, const sgm::Pattern_Automorphism* ruleSymmetry
							, const bool ignoreAtomClassLabel
							, const bool noRedundantMolecules
							)
{
	assert(ruleComponent <= rulePattern.getFirstOfEachComponent().size() /* ruleComponent exeeds number of rule components */);
	
	  // check for recursion abort
	  // --> all component targets have been determined and allTargets is full
	if (ruleComponent == rulePattern.getFirstOfEachComponent().size()) {

//		std::cout <<"## final target : noRed("<<(noRedundantMolecules?"true":"false") <<") ";
//		for (size_t i=0; i<curTargets.size(); i++) {
//			std::cout <<curTargets.at(i)<<" ";
//		}
//		std::cout <<std::endl;

		 // represent all target molecules as one graph to search
		sgm::Graph_boostV_p< ggl::chem::Molecule > targets(curTargets);
		 // create final target graph to search
		sgm::Graph_Interface* finalTargets = ignoreAtomClassLabel
						? (sgm::Graph_Interface*)new ggl::chem::Molecule_Graph_noClass( targets )
						: (sgm::Graph_Interface*)&targets;
		 // create final match reported for rule application
		sgm::Match_Reporter* mrApplyRuleFinal = ignoreAtomClassLabel
						? new sgm::MR_restoreNodeLabel( mrApplyRule )
						: &mrApplyRule;
		if (ruleSymmetry != NULL) {
			  // set up symmetry breaking interface
			sgm::MR_SymmBreak mrSymmBreak( *ruleSymmetry, *mrApplyRuleFinal);
			  // find all matches and apply rule
			sgm.findMatches( rulePattern, *finalTargets, mrSymmBreak, UINT_MAX );
		} else {
			  // find all matches and apply rule
			sgm.findMatches( rulePattern, *finalTargets, *mrApplyRuleFinal, UINT_MAX );
		}
		  // cleanup
		if (ignoreAtomClassLabel) { delete finalTargets; delete mrApplyRuleFinal; }
		  // all done --> end recursion here
		return true;
	}
	
	assert(curTargets.size() > ruleComponent /* not enough room to add target for current component */);
	
	  // get the label of the 'ruleComponent's component of the rule
	ggl::chem::LeftSidePattern::IndexSet::const_iterator firstIt 
			= rulePattern.getFirstOfEachComponent().begin();
	for (size_t c=0; c<ruleComponent; ++c) {
		firstIt++;
	}
	  // get label of the 'ruleComponent's component
	const int curCompLabel = rulePattern.getComponentLabeling().at(*firstIt);
	
	  // set up subgraph of the current rule component for checking its presence
	  // in target graphs
	sgm::SubGraphPattern component(	rulePattern
									, rulePattern.getComponentLabeling()
									, curCompLabel
									, rulePattern.getConstraints()
									, *(rulePattern.getUsedWildcard()) );
	
	  // assume this component is not present until not proven
	bool thisComponentPresent = false;
	  // dummy for matching
	sgm::MR_Counting mrCounting;
	  // for all targets
	for (	SMILES_container::const_iterator it=allTargets.begin(); 
			it!=allTargets.end() ; ++it ) 
	{
		// check if molecule is not already present
		if (noRedundantMolecules) {
			bool curMolNotPresent = true;
			for (size_t p=0; curMolNotPresent && p<ruleComponent; p++)  {
				// check if current molecule is not already present as a component
				curMolNotPresent = (curTargets.at(p) != it->second);
			}
			// check if molecule is already present
			if (!curMolNotPresent) {
				// skip this molecule and continue with next
				continue;
			}
		}
		
		  // set up wrapper for current target for matching
		sgm::Graph_boost< ggl::chem::Molecule > curTarget(*(it->second));
		 // create final target graph to search
		sgm::Graph_Interface* finalTarget = ignoreAtomClassLabel
						? (sgm::Graph_Interface*)new ggl::chem::Molecule_Graph_noClass( curTarget )
						: (sgm::Graph_Interface*)&curTarget;
		  // check if current target contains this component at least once
		const size_t numOfMatches = sgm.findMatches( component, *finalTarget, mrCounting,1 );
		  // cleanup
		if (ignoreAtomClassLabel) { delete finalTarget; }
		  // check if component was found
		if (numOfMatches == 0) {
			// this component is not contained in current target
			// --> go to next target
			continue;
		} else {
			thisComponentPresent = true;
		}
		  // add current target to the overall graph to apply this rule to
		curTargets[ruleComponent] = it->second;
		  // recursive call to handle next component or to start matching
		const bool remMatched = singleRuleApplicationRec(	sgm
															, rulePattern
															, (ruleComponent+1)
															, allTargets
															, curTargets
															, mrApplyRule
															, ruleSymmetry
															, ignoreAtomClassLabel
															, noRedundantMolecules
															);
		if (!remMatched)
			return false;
	}	
	
	  // return if this component or one of the other was present to allow for
	  // early abortion if one of the component cannot be matched
	return thisComponentPresent;
}


//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////

bool 
singleRuleApplication(	sgm::SubGraphMatching& sgm
						, const ggl::chem::LeftSidePattern& rulePattern
						, const SMILES_container & newTargets
						, const SMILES_container & oldTargets
						, sgm::Match_Reporter& mrApplyRule
						, const sgm::Pattern_Automorphism* ruleSymmetry
						, const bool ignoreAtomClassLabel
						, const bool noRedundantMolecules
						)
{
	  // target graph holding compNbr different targets
	std::vector< const ggl::chem::Molecule* > 
			curTargets(rulePattern.getFirstOfEachComponent().size(),NULL);

	  // forward job to recursive version
	return singleRuleApplicationRec(	sgm
										, rulePattern
										, 0
										, newTargets
										, oldTargets
										, false
										, curTargets
										, mrApplyRule
										, ruleSymmetry
										, ignoreAtomClassLabel
										, noRedundantMolecules
										);
}


//////////////////////////////////////////////////////////////////////////

int 
singleRuleApplicationRec(	sgm::SubGraphMatching& sgm
							, const ggl::chem::LeftSidePattern& rulePattern
							, const size_t ruleComponent
							, const SMILES_container & newTargets
							, const SMILES_container & oldTargets
							, const bool atLeastOneNewAdded
							, std::vector< const ggl::chem::Molecule* >& curTargets
							, sgm::Match_Reporter& mrApplyRule
							, const sgm::Pattern_Automorphism* ruleSymmetry
							, const bool ignoreAtomClassLabel
							, const bool noRedundantMolecules
							)
{
	assert(ruleComponent <= rulePattern.getFirstOfEachComponent().size() /* ruleComponent exeeds number of rule components */);
	
	  // check for recursion abort
	  // --> all component targets have been determined and allTargets is full
	if (ruleComponent == rulePattern.getFirstOfEachComponent().size()) {
		 // represent all target molecules as one graph to search
		sgm::Graph_boostV_p< ggl::chem::Molecule > targets(curTargets);
		 // create final target graph to search
		sgm::Graph_Interface* finalTargets = ignoreAtomClassLabel
						? (sgm::Graph_Interface*)new ggl::chem::Molecule_Graph_noClass( targets )
						: (sgm::Graph_Interface*)&targets;
		 // create final match reported for rule application
		sgm::Match_Reporter* mrApplyRuleFinal = ignoreAtomClassLabel
						? new sgm::MR_restoreNodeLabel( mrApplyRule )
						: &mrApplyRule;

		if (ruleSymmetry != NULL) {
			  // set up symmetry breaking interface
			sgm::MR_SymmBreak mrSymmBreak( *ruleSymmetry, *mrApplyRuleFinal);
			  // find all matches and apply rule
			sgm.findMatches( rulePattern, *finalTargets, mrSymmBreak, UINT_MAX );
		} else {
			  // find all matches and apply rule
			sgm.findMatches( rulePattern, *finalTargets, *mrApplyRuleFinal, UINT_MAX );
		}

		// cleanup if necessary
		if (ignoreAtomClassLabel) {
			delete finalTargets;
			delete mrApplyRuleFinal;
		}

		  // all done --> end recursion here
		return 0;
	}
	
	assert(curTargets.size() > ruleComponent /* not enough room to add target for current component */);
	
	  // get the label of the 'ruleComponent's component of the rule
	ggl::chem::LeftSidePattern::IndexSet::const_iterator firstIt 
			= rulePattern.getFirstOfEachComponent().begin();
	for (size_t c=0; c<ruleComponent; ++c) {
		firstIt++;
	}
	  // get label of the 'ruleComponent's component
	const int curCompLabel = rulePattern.getComponentLabeling().at(*firstIt);
	
	  // set up subgraph of the current rule component for checking its presence
	  // in target graphs
	sgm::SubGraphPattern component(	rulePattern
									, rulePattern.getComponentLabeling()
									, curCompLabel
									, rulePattern.getConstraints()
									, *(rulePattern.getUsedWildcard()) );
	
	  // for the presence check of the current component
	sgm::MR_Counting mrCounting;
	
	  // assume this component is not present until not proven
	bool thisComponentPresent = false;
	  // assume remaining components are present (set by the recursion)
	bool otherComponentsPresent = true;
	
	  // for new targets
	for (	SMILES_container::const_iterator it=newTargets.begin(); 
			otherComponentsPresent && it!=newTargets.end() ; ++it ) 
	{
		// check if molecule is not already present
		if (noRedundantMolecules) {
			bool curMolNotPresent = true;
			for (size_t p=0; curMolNotPresent && p<ruleComponent; p++)  {
				// check if current molecule is not already present as a component
				curMolNotPresent = (curTargets.at(p) != it->second);
			}
			// check if molecule is already present
			if (!curMolNotPresent) {
				// skip this molecule and continue with next
				continue;
			}
		}
		
		  // set up wrapper for current target for matching
		sgm::Graph_boost< ggl::chem::Molecule > curTarget(*(it->second));
		 // create final target graph to search
		sgm::Graph_Interface* finalTarget = ignoreAtomClassLabel
						? (sgm::Graph_Interface*)new ggl::chem::Molecule_Graph_noClass( curTarget )
						: (sgm::Graph_Interface*)&curTarget;
		  // check if current target contains this component at least once
		mrCounting.resetHits();
		sgm.findMatches( component, *finalTarget, mrCounting, 1 );
		  // cleanup
		if (ignoreAtomClassLabel) { delete finalTarget; }
		  // check if component was found
		if (mrCounting.getHits() == 0) {
			// this component is not contained in current target
			// --> go to next target
			continue;
		} else {
			thisComponentPresent = true;
		}

		  // add current target to the overall graph to apply this rule to
		curTargets[ruleComponent] = it->second;
		  // recursive call to handle next component or to start matching
		int recCall = singleRuleApplicationRec(	sgm
															, rulePattern
															, (ruleComponent+1)
															, newTargets
															, oldTargets
															, true
															, curTargets
															, mrApplyRule
															, ruleSymmetry
															, ignoreAtomClassLabel
															, noRedundantMolecules
															);
		otherComponentsPresent = (recCall == 0);
	}	
	
	  // avoid rule application if none of the newTargets was mapped onto one
	  // of the components
	if (	(ruleComponent+1) == rulePattern.getFirstOfEachComponent().size()
			&& !atLeastOneNewAdded && !thisComponentPresent ) 
	{
		return 2;
	}
	
	bool oldMolsPresent = true;
	
	  // for old targets
	for (	SMILES_container::const_iterator it=oldTargets.begin(); 
			oldMolsPresent && it!=oldTargets.end() ; ++it ) 
	{
		// check if molecule is not already present
		if (noRedundantMolecules) {
			bool curMolNotPresent = true;
			for (size_t p=0; curMolNotPresent && p<ruleComponent; p++)  {
				// check if current molecule is not already present as a component
				curMolNotPresent = (curTargets.at(p) != it->second);
			}
			// check if molecule is already present
			if (!curMolNotPresent) {
				// skip this molecule and continue with next
				continue;
			}
		}
		
		  // set up wrapper for current target for matching
		sgm::Graph_boost< ggl::chem::Molecule > curTarget(*(it->second));
		 // create final target graph to search
		sgm::Graph_Interface* finalTarget = ignoreAtomClassLabel
						? (sgm::Graph_Interface*)new ggl::chem::Molecule_Graph_noClass( curTarget )
						: (sgm::Graph_Interface*)&curTarget;
		 // check if current target contains this component at least once
		mrCounting.resetHits();
		sgm.findMatches( component, *finalTarget, mrCounting, 1 );
		  // cleanup
		if (ignoreAtomClassLabel) { delete finalTarget; }
		  // check if component was found
		if (mrCounting.getHits() == 0) {
			// this component is not contained in current target
			// --> go to next target
			continue;
		} else {
			thisComponentPresent = true;
		}
		  // add current target to the overall graph to apply this rule to
		curTargets[ruleComponent] = it->second;
		  // recursive call to handle next component or to start matching
		int recCall = singleRuleApplicationRec(	sgm
															, rulePattern
															, (ruleComponent+1)
															, newTargets
															, oldTargets
															, atLeastOneNewAdded
															, curTargets
															, mrApplyRule
															, ruleSymmetry
															, ignoreAtomClassLabel
															, noRedundantMolecules
															);
		switch (recCall) {
		case 0 : oldMolsPresent = true;  break;
		case 1 : oldMolsPresent = false; break;
		case 2 : return 2;
		}
	}	
	
	if (thisComponentPresent && (otherComponentsPresent || oldMolsPresent))
		return 0;
	else
		return 1;
	
//	  // return if this component or one of the other was present to allow for
//	  // early abortion if one of the components cannot be matched
//	return thisComponentPresent && otherComponentsPresent;
}


//////////////////////////////////////////////////////////////////////////

#include <sgm/SGM_vf2.hh>

//////////////////////////////////////////////////////////////////////////

void
applyRules( const RulePatternMap & rules
			, const SMILES_container & initialMolecules
			, SMILES_container & producedMolecules
			, ggl::chem::MR_Reactions::Reaction_Container & producedReactions
			, const ggl::chem::ReactionRateCalculation * rateCalc
			, const ggl::chem::AromaticityPerception & aromaticity
			, const bool ignoreAtomClassLabel
			, const bool noRedundantMolecules
		)
{
	  // set up graph matcher
	sgm::SGM_vf2 matcher;

	  // FOR EACH set of rules with equal number of connected components
	  //          in leftsidepattern :
	  // find matches and apply rules
	for (RulePatternMap::const_iterator pat = rules.begin();
			pat != rules.end(); ++pat)
	{
		  // set up Rule applyer 
		ggl::chem::MR_Reactions mr( initialMolecules, producedMolecules
						, producedReactions
//						, pat->first
						, true
						, rateCalc
						, aromaticity
						, ignoreAtomClassLabel );
		

		// for all rules with given number of components
		for (size_t curRule=0; curRule<pat->second.size(); ++curRule) {
			  // set up symmetry breaking conditions for current rule
			sgm::PA_OrderCheck ga(pat->second.at(curRule)->getGraphAutomorphism());
			  // rule application
			singleRuleApplication(  matcher
									, *(pat->second.at(curRule))
									, initialMolecules
									, mr
									, &ga
									, ignoreAtomClassLabel
									, noRedundantMolecules
									);
		}

	}

}
	
//////////////////////////////////////////////////////////////////////////

void
applyRules( const RulePatternMap & rules
			, const SMILES_container & oldMolecules
			, const SMILES_container & newMolecules
			, SMILES_container & producedMolecules
			, ggl::chem::MR_Reactions::Reaction_Container& producedReactions
			, const ggl::chem::ReactionRateCalculation * rateCalc
			, const bool allowAllIntra
			, const ggl::chem::AromaticityPerception & aromaticity
			, const bool ignoreAtomClassLabel
			, const bool noRedundantMolecules
		)
{
	  // set up graph matcher
	sgm::SGM_vf2 matcher;

	  // FOR EACH set of rules with equal number of connected components
	  //          in leftsidepattern :
	  // find matches and apply rules
	for (RulePatternMap::const_iterator pat = rules.begin();
			pat != rules.end(); ++pat)
	{
		  // check if single component rules are currently handled
		if (pat->first == 1) {
			
			// --> apply rules to newMolecules only
			
			  // set up Rule applyer 
			ggl::chem::MR_Reactions mr( newMolecules, oldMolecules, producedMolecules
							, producedReactions
							, !allowAllIntra
//							, pat->first, false
							, rateCalc
							, aromaticity
							, ignoreAtomClassLabel
							);
			

			// for all rules with given number of components
			for (size_t curRule=0; curRule<pat->second.size(); ++curRule) {
				  // set up symmetry breaking conditions for current rule
				sgm::PA_OrderCheck ga(pat->second.at(curRule)->getGraphAutomorphism());
				  // rule application
				singleRuleApplication(  matcher
										, *(pat->second.at(curRule))
										, newMolecules
										, mr
										, &ga
										, ignoreAtomClassLabel
										, noRedundantMolecules
										);
			}

		} else {
			
			// --> multi component rule --> at least one component has to be
			//     mapped on a newMolecules entry
			
			  // set up Rule applyer 
			ggl::chem::MR_Reactions mr( oldMolecules, newMolecules, producedMolecules
							, producedReactions
							, !allowAllIntra
//							, (allowAllIntra ? 1 : pat->first) 
//							, false
							, rateCalc
							, aromaticity
							, ignoreAtomClassLabel
							);
			
			// for all rules with given number of components
			for (size_t curRule=0; curRule<pat->second.size(); ++curRule) {
				  // set up symmetry breaking conditions for current rule
				sgm::PA_OrderCheck ga(pat->second.at(curRule)->getGraphAutomorphism());
				  // rule application
				singleRuleApplication(  matcher
										, *(pat->second.at(curRule))
										, newMolecules
										, oldMolecules
										, mr
										, &ga
										, ignoreAtomClassLabel
										, noRedundantMolecules
										);
			}
		}
	}
	
}

	
//////////////////////////////////////////////////////////////////////////

#include <ggl/MR_ApplyRule.hh>
#include <ggl/chem/GS_SMILES.hh>

//////////////////////////////////////////////////////////////////////////

///* Applies all rules onto a set of initial molecules. The resulting 
// * molecules from the applications as well as the reaction information is
// * written to provided containers.
// * 
// * @param rules (IN) the left side pattern of the rules to apply 
// * @param oldMolecules (IN) molecules to that the rules have been applied
// *        already (e.g. in last iteration), such that no rule application 
// *        to this set of rules only is done
// * @param newMolecules (IN) molecules the rules were not applied on 
// *        already such that at least one of these molecules will be within
// *        a single rule application target
// * @param producedMolecules (OUT) the container where the molecules 
// *         produced by the rule application are added
// */
//void
//applyRules( const RulePatternMap & rules
//			, const SMILES_container & oldMolecules
//			, const SMILES_container & newMolecules
//			, SMILES_container & producedMolecules
//		)
//{
//	  // set up graph matcher
//	sgm::SGM_vf2 matcher;
//	
//	// stores produced smiles if not existing in one of these storages
//	ggl::chem::GS_SMILES_MOLp<SMILES_container> gs(producedMolecules, oldMolecules, newMolecules);
//
//	  // FOR EACH set of rules with equal number of connected components
//	  //          in leftsidepattern :
//	  // find matches and apply rules
//	for (RulePatternMap::const_iterator pat = rules.begin();
//			pat != rules.end(); ++pat)
//	{
//		  // set up Rule applyer 
//		ggl::MR_ApplyRule mr( gs, pat->first, false);
//		  // check if single component rules are currently handled
//		if (pat->first == 1) {
//			
//			// --> apply rules to newMolecules only
//			
//			// for all rules with given number of components
//			for (size_t curRule=0; curRule<pat->second.size(); ++curRule) {
//				  // set up symmetry breaking conditions for current rule
//				sgm::PA_OrderCheck ga(pat->second.at(curRule)->getGraphAutomorphism());
//				  // rule application
//				singleRuleApplication(  matcher
//										, *(pat->second.at(curRule))
//										, newMolecules
//										, mr
//										, &ga );
//			}
//		} else {
//			
//			// --> multi component rule --> at least one component has to be
//			//     mapped on a newMolecules entry
//			
//			// for all rules with given number of components
//			for (size_t curRule=0; curRule<pat->second.size(); ++curRule) {
//				  // set up symmetry breaking conditions for current rule
//				sgm::PA_OrderCheck ga(pat->second.at(curRule)->getGraphAutomorphism());
//				  // rule application
//				singleRuleApplication(  matcher
//										, *(pat->second.at(curRule))
//										, newMolecules
//										, oldMolecules
//										, mr
//										, &ga );
//			}
//		}
//	}
//	
//}


//////////////////////////////////////////////////////////////////////////


size_t
correctInputMolecules( const SMILES_container & targetSmiles
						, SMILES_container & producedSmiles
						, ggl::chem::AromaticityPerception * aromaticity_full
						, const bool protonFilling
						, const size_t setNextAtomClass
						)
{
	using namespace ggl;
	using namespace ggl::chem;

	size_t nextAtomClass = setNextAtomClass;

	GS_SMILES_MOLp<SMILES_container> storage(producedSmiles);
	GS_MolCheck *molChecker = NULL;

	  // ensure predictor for proton-pruned molecules is present
	if(aromaticity_full != NULL) {
		molChecker = new GS_MolCheck( storage, *aromaticity_full, false );
	}


	Molecule mol;
	for (SMILES_container::const_iterator it=targetSmiles.begin(); it!=targetSmiles.end() ; ++it ) {

		try {
			  // copy orginal molecule
			mol = Molecule(*(it->second));

			  // fill protons
			if (protonFilling) {
				  // fill
				MoleculeUtil::fillProtons( mol );
				if (nextAtomClass > 0) {
					 // add atom class if needed
					nextAtomClass = setAtomClass( mol, nextAtomClass, false );
				}
			}

			if (molChecker != NULL) {
				  // correct aromaticity etc. and store if unknown
				molChecker->addMolecule( mol );
			} else {
				  // store if unknown
				storage.addMolecule( mol );
			}

		} catch (std::exception & ex) {
			std::cerr <<"\n\n ERROR for input molecule '"
					<<it->first
					<<"': " <<ex.what() <<"\n"
					<<"--> will be ignored !"<<std::endl;
		}
	}

	if (molChecker != NULL) {
		delete molChecker; molChecker = NULL;
	}

	return nextAtomClass;
}


//////////////////////////////////////////////////////////////////////////

size_t
setAtomClass( ggl::chem::Molecule & mol
		, const size_t nextClassLabelOption
		, const bool overwriteExistingClass )
{
	// get next atom class that is non-redundant
	size_t nextClassLabel = nextClassLabelOption;

	// access to node labels of graph
	boost::property_map< ggl::chem::Molecule, ggl::PropNodeLabel >::type
		nodeLabel = boost::get( ggl::PropNodeLabel(), mol );

	  // update property map
	typename ggl::chem::Molecule::vertex_iterator vi, v_end;
	boost::tie(vi,v_end) = boost::vertices(mol);

	// relabel all nodes not already with class label
	// i.e. add nextClassLabel suffix
	while (vi != v_end) {
		// ensure the class label was not known
		const size_t classLabel = ggl::chem::MoleculeUtil::getClass( nodeLabel[*vi] );
		// check if not labeled yet
		if ( classLabel == 0 ) {
			// set class label
			nodeLabel[*vi] += ":" + boost::lexical_cast<std::string>( nextClassLabel++ );
		} else if (overwriteExistingClass) {
			std::cerr <<"# overwriting atom class label "<<nodeLabel[*vi] <<" with "<<nextClassLabel<<std::endl;
			// replace class label
			nodeLabel[*vi]
			        = nodeLabel[*vi].substr(0,nodeLabel[*vi].find(":"))
					+ ":" + boost::lexical_cast<std::string>( nextClassLabel++ );
		}
		vi++;
	}
	return nextClassLabel;
}

//////////////////////////////////////////////////////////////////////////

