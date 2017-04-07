
#include <exception>
#include <algorithm>
#include <stdexcept>

#include "ggl/chem/MoleculeUtil.hh"

#include "sgm/SGM_vf2.hh"
#include "sgm/Pattern.hh"
#include "sgm/MC_Node.hh"
#include "sgm/MR_Storing.hh"

#include <boost/assign/list_of.hpp> // for 'list_of()'

#ifndef NDEBUG
#include "ggl/chem/SMILESwriter.hh"
#endif


namespace ggl {
	namespace chem {
	
//##############################################################################

	MoleculeUtil::AtomDataMap
	MoleculeUtil::atomData =	MoleculeUtil::AtomDataMap();

	MoleculeUtil::AromaticSwapMap
	MoleculeUtil::aromaticSwapMap =	MoleculeUtil::AromaticSwapMap();

	MoleculeUtil::BondDataMap
	MoleculeUtil::bondData =	MoleculeUtil::BondDataMap();
 
	const size_t 
	MoleculeUtil::C_Consistent = 0;

	const size_t 
	MoleculeUtil::C_AtomLabelInvalid = 5;

	const size_t 
	MoleculeUtil::C_BondLabelInvalid = 7;

	const size_t
	MoleculeUtil::C_AtomLabelWildcard = 11;

	const size_t
	MoleculeUtil::C_AtomComplexWithH = 13;

	const size_t
	MoleculeUtil::C_AtomValence = 17;

	const size_t
	MoleculeUtil::C_BondLoop = 19;

	const size_t
	MoleculeUtil::C_NonConnected = 23;

	const std::string 
	MoleculeUtil::DIGIT = std::string("1234567890");

	const std::string
	MoleculeUtil::WHITESPACE = std::string(" \t\n");
 	
	const std::string 
	MoleculeUtil::AtomLabelWildcard = std::string("*");

//##############################################################################

	// atomic orbitals : s, p,  d,  f
	// orbital sizes   : 2, 6, 10, 14  = shellSize(i)-shellSize(i-1)
	// (as number of electrons within)

	// Madelung rule / Klechkowski rule http://en.wikipedia.org/wiki/Aufbau_principle
	// "Aufbau principle" 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, ...
	// orbit electrons :  +2  +2  +6  +2  +6  +2 +10  +6  +2 +10  +6  +2 +14 +10  +6  +2
	// electron sum :      2,  4, 10, 12, 18, 20, 30, 36, 38, 48, 54, 56, 60, 70, 76, 78, ...
	// outmost shell :     1   2   2   3   3   4   4   4   5   5   5   6   6   6   6   7
	// outmost shell size  2   2   8   2   8   2   2   8   2   2   8   2   2   2   8   2

	// check http://www.talete.mi.it/help/dproperties_help/h_depleted_structures.htm


	// shell sizes : 2, 8, 18, 32, 50   = 2*n^2
	std::vector<MoleculeUtil::OneByte>
	MoleculeUtil::AtomLabelData::shellSizeSums = boost::assign::list_of(2)(10)(28)(60)(110);


//##############################################################################



	const
	MoleculeUtil
	::AtomDataMap&
	MoleculeUtil
	::getAtomData( void )
	{
		/////////////////////////////////////////////////////////////////
		// NOTE : EVERY SYMBOL ADDED HERE HAS TO BE ADDED TO THE SMILES
		//        PARSER
		/////////////////////////////////////////////////////////////////

		  // check if first call and data has to be set first
		if (atomData.size() == 0) {
			typedef AtomDataMap::value_type VT;
			typedef AromaticSwapMap::value_type MT;
			/*
			  prefered oxidation states are taken from
			  http://environmentalchemistry.com/yogi/periodic/FOO.html
			*/
			atomData.insert( VT(AtomLabelWildcard
										, AtomLabelData(  0,0,0,  0.000000, 0 )));
			atomData.insert( VT("H" 	, AtomLabelData(  1,1,0,  1.007940, 1 )));
			atomData.insert( VT("He" 	, AtomLabelData(  2,0,0,  4.002602, 1 )));
			atomData.insert( VT("Li" 	, AtomLabelData(  3,1,0,  6.941000, 1 )));
			atomData.insert( VT("Be" 	, AtomLabelData(  4,2,0,  9.012182, 1 )));
			atomData.insert( VT("B" 	, AtomLabelData(  5,3,0, 10.811000, 1 )));
			atomData.insert( VT("b" 	, AtomLabelData(  5,3,1, 10.811000, 1 )));
				aromaticSwapMap.insert( MT("b", "B" ));
				aromaticSwapMap.insert( MT("B", "b" ));
			atomData.insert( VT("C" 	, AtomLabelData(  6,4,0, 12.010700, 1 )));
			atomData.insert( VT("c" 	, AtomLabelData(  6,4,1, 12.010700, 1 )));
				aromaticSwapMap.insert( MT("c", "C" ));
				aromaticSwapMap.insert( MT("C", "c" ));
			atomData.insert( VT("N" 	, AtomLabelData(  7,3,0, 14.006740, 1, boost::assign::list_of(3)(5)(7) )));
			atomData.insert( VT("n" 	, AtomLabelData(  7,3,1, 14.006740, 1, boost::assign::list_of(3)(5)(7) )));
				aromaticSwapMap.insert( MT("n", "N" ));
				aromaticSwapMap.insert( MT("N", "n" ));
			atomData.insert( VT("O" 	, AtomLabelData(  8,2,0, 15.999400, 1 )));
			atomData.insert( VT("o" 	, AtomLabelData(  8,2,1, 15.999400, 1 )));
				aromaticSwapMap.insert( MT("o", "O" ));
				aromaticSwapMap.insert( MT("O", "o" ));
			atomData.insert( VT("F" 	, AtomLabelData(  9,1,0, 18.998402, 1 )));
			atomData.insert( VT("Ne" 	, AtomLabelData( 10,0,0, 20.179700, 1 )));
			atomData.insert( VT("Na" 	, AtomLabelData( 11,1,0, 22.989770, 1 )));
			atomData.insert( VT("Mg" 	, AtomLabelData( 12,2,0, 24.305000, 0 )));
			atomData.insert( VT("Al" 	, AtomLabelData( 13,3,0, 26.981538, 0 )));
			atomData.insert( VT("Si" 	, AtomLabelData( 14,4,0, 28.085500, 0 )));
			atomData.insert( VT("P" 	, AtomLabelData( 15,5,0, 30.973761, 1, boost::assign::list_of(3)(5)(7) )));
			atomData.insert( VT("p" 	, AtomLabelData( 15,5,1, 30.973761, 1, boost::assign::list_of(3)(5)(7) )));
				aromaticSwapMap.insert( MT("p", "P" ));
				aromaticSwapMap.insert( MT("P", "p" ));
			atomData.insert( VT("S" 	, AtomLabelData( 16,2,0, 32.066000, 1, boost::assign::list_of(2)(4)(6) )));
			atomData.insert( VT("s" 	, AtomLabelData( 16,3,1, 32.066000, 1, boost::assign::list_of(3)(5) )));  // "additional electron" for aromaticity check
				aromaticSwapMap.insert( MT("s", "S" ));
				aromaticSwapMap.insert( MT("S", "s" ));
			atomData.insert( VT("Cl" 	, AtomLabelData( 17,1,0, 35.452700, 1 )));
			atomData.insert( VT("Ar" 	, AtomLabelData( 18,0,0, 39.948000, 0 )));
			atomData.insert( VT("K" 	, AtomLabelData( 19,1,0, 39.098300, 0 )));
			atomData.insert( VT("Ca" 	, AtomLabelData( 20,2,0, 40.078000, 0 )));
			atomData.insert( VT("Sc" 	, AtomLabelData( 21,3,0, 44.955910, 0 )));
			atomData.insert( VT("Ti" 	, AtomLabelData( 22,4,0, 47.867000, 0 )));
			atomData.insert( VT("V" 	, AtomLabelData( 23,5,0, 50.941500, 0 )));
			atomData.insert( VT("Cr" 	, AtomLabelData( 24,3,0, 51.996100, 0 )));
			atomData.insert( VT("Mn" 	, AtomLabelData( 25,2,0, 54.938049, 0 )));
			atomData.insert( VT("Fe" 	, AtomLabelData( 26,3,0, 55.845000, 0 )));
			atomData.insert( VT("Co" 	, AtomLabelData( 27,2,0, 58.933200, 0 )));
			atomData.insert( VT("Ni" 	, AtomLabelData( 28,2,0, 58.693400, 0 )));
			atomData.insert( VT("Cu" 	, AtomLabelData( 29,2,0, 63.546000, 0 )));
			atomData.insert( VT("Zn" 	, AtomLabelData( 30,2,0, 65.390000, 0 )));
			atomData.insert( VT("Ga" 	, AtomLabelData( 31,3,0, 69.723000, 0 )));
			atomData.insert( VT("Ge" 	, AtomLabelData( 32,4,0, 72.610000, 0 )));
			atomData.insert( VT("As" 	, AtomLabelData( 33,3,0, 74.921600, 1, boost::assign::list_of(3)(5)(7) )));
			atomData.insert( VT("as" 	, AtomLabelData( 33,3,1, 74.921600, 1, boost::assign::list_of(3)(5)(7) )));
				aromaticSwapMap.insert( MT("as", "As" ));
				aromaticSwapMap.insert( MT("As", "as" ));
//			atomData.insert( VT("Se" 	, AtomLabelData( 34,4,0, 78.960000, 1, boost::assign::list_of(2)(4)(6) )));
//			atomData.insert( VT("se" 	, AtomLabelData( 34,4,1, 78.960000, 1, boost::assign::list_of(3)(5) )));  // "additional electron" for aromaticity check
				// for hydrogen filling: "standard valence" changed to 2 instead of 4 as given by http://environmentalchemistry.com/yogi/periodic/FOO.html
			atomData.insert( VT("Se" 	, AtomLabelData( 34,2,0, 78.960000, 1, boost::assign::list_of(2)(4)(6) )));
			atomData.insert( VT("se" 	, AtomLabelData( 34,3,1, 78.960000, 1, boost::assign::list_of(3)(5) )));  // "additional electron" for aromaticity check
				aromaticSwapMap.insert( MT("se", "Se" ));
				aromaticSwapMap.insert( MT("Se", "se" ));
			atomData.insert( VT("Br" 	, AtomLabelData( 35,1,0, 79.904000, 1 )));
			atomData.insert( VT("Kr" 	, AtomLabelData( 36,0,0, 83.800000, 0 )));
			atomData.insert( VT("Rb" 	, AtomLabelData( 37,1,0, 85.467800, 0 )));
			atomData.insert( VT("Sr" 	, AtomLabelData( 38,3,0, 87.620000, 0 )));
			atomData.insert( VT("Y" 	, AtomLabelData( 39,3,0, 88.905850, 0 )));
			atomData.insert( VT("Zr" 	, AtomLabelData( 40,4,0, 91.224000, 0 )));
			atomData.insert( VT("Nb" 	, AtomLabelData( 41,5,0, 92.906380, 0 )));
			atomData.insert( VT("Mo" 	, AtomLabelData( 42,6,0, 95.940000, 0 )));
			atomData.insert( VT("Tc" 	, AtomLabelData( 43,7,0, 98.906300, 0 )));
			atomData.insert( VT("Ru" 	, AtomLabelData( 44,4,0,101.070000, 0 )));
			atomData.insert( VT("Rh" 	, AtomLabelData( 45,3,0,102.905500, 0 )));
			atomData.insert( VT("Pd" 	, AtomLabelData( 46,2,0,106.420000, 0 )));
			atomData.insert( VT("Ag" 	, AtomLabelData( 47,1,0,107.868200, 0 )));
			atomData.insert( VT("Cd" 	, AtomLabelData( 48,2,0,112.411000, 0 )));
			atomData.insert( VT("In" 	, AtomLabelData( 49,3,0,114.818000, 0 )));
			atomData.insert( VT("Sn" 	, AtomLabelData( 50,4,0,118.710000, 0 )));
			atomData.insert( VT("Sb" 	, AtomLabelData( 51,3,0,121.760000, 0 )));
			atomData.insert( VT("Te" 	, AtomLabelData( 52,4,0,127.600000, 0 )));
			atomData.insert( VT("I" 	, AtomLabelData( 53,1,0,126.904470, 0 )));
			atomData.insert( VT("Xe" 	, AtomLabelData( 54,0,0,131.290000, 0 )));
			atomData.insert( VT("Cs" 	, AtomLabelData( 55,1,0,132.905450, 0 )));
			atomData.insert( VT("Ba" 	, AtomLabelData( 56,2,0,137.327000, 0 )));
			// 57-70 lanthanides
			atomData.insert( VT("Lu" 	, AtomLabelData( 71,3,0,174.967000, 0 )));
			atomData.insert( VT("Hf" 	, AtomLabelData( 72,4,0,178.490000, 0 )));
			atomData.insert( VT("Ta" 	, AtomLabelData( 73,5,0,180.947900, 0 )));
			atomData.insert( VT("W" 	, AtomLabelData( 74,6,0,183.840000, 0 )));
			atomData.insert( VT("Re" 	, AtomLabelData( 75,6,0,186.207000, 0 )));
			atomData.insert( VT("Os" 	, AtomLabelData( 76,4,0,190.230000, 0 )));
			atomData.insert( VT("Ir" 	, AtomLabelData( 77,4,0,192.217000, 0 )));
			atomData.insert( VT("Pt" 	, AtomLabelData( 78,4,0,195.078000, 0 )));
			atomData.insert( VT("Au" 	, AtomLabelData( 79,3,0,196.966550, 0 )));
			atomData.insert( VT("Hg" 	, AtomLabelData( 80,2,0,200.590000, 0 )));
			atomData.insert( VT("Tl" 	, AtomLabelData( 81,1,0,204.383300, 0 )));
			atomData.insert( VT("Pb" 	, AtomLabelData( 82,2,0,207.200000, 0 )));
			atomData.insert( VT("Bi" 	, AtomLabelData( 83,3,0,208.980380, 0 )));
			atomData.insert( VT("Po" 	, AtomLabelData( 84,4,0,208.982400, 0 )));
			atomData.insert( VT("At" 	, AtomLabelData( 85,1,0,209.987100, 0 )));
			atomData.insert( VT("Rn" 	, AtomLabelData( 86,0,0,222.017600, 0 )));
			atomData.insert( VT("Fr" 	, AtomLabelData( 87,1,0,223.019700, 0 )));
			atomData.insert( VT("Ra" 	, AtomLabelData( 88,2,0,226.025400, 0 )));
			// 89-102 actinides
			atomData.insert( VT("Lr" 	, AtomLabelData(103,3,0,262.110000, 0 )));
			atomData.insert( VT("Rf" 	, AtomLabelData(104,4,0,261.108900, 0 )));
			atomData.insert( VT("Db" 	, AtomLabelData(105,4,0,262.114400, 0 )));
			atomData.insert( VT("Sg" 	, AtomLabelData(106,2,0,263.118600, 0 )));
			atomData.insert( VT("Bh" 	, AtomLabelData(107,2,0,264.120000, 0 )));
			atomData.insert( VT("Hs" 	, AtomLabelData(108,2,0,265.130600, 0 )));
			atomData.insert( VT("Mt" 	, AtomLabelData(109,2,0,268.000000, 0 )));
			atomData.insert( VT("Uun" 	, AtomLabelData(110,2,0,269.000000, 0 )));
			atomData.insert( VT("Uuu" 	, AtomLabelData(111,2,0,272.000000, 0 )));
			atomData.insert( VT("Uub" 	, AtomLabelData(112,2,0,277.000000, 0 )));
			// Lanthanide Series
			atomData.insert( VT("La" 	, AtomLabelData( 57,3,0,138.905500, 0 )));
			atomData.insert( VT("Ce" 	, AtomLabelData( 58,3,0,140.116000, 0 )));
			atomData.insert( VT("Pr" 	, AtomLabelData( 59,3,0,140.907650, 0 )));
			atomData.insert( VT("Nd" 	, AtomLabelData( 60,3,0,144.240000, 0 )));
			atomData.insert( VT("Pm" 	, AtomLabelData( 61,3,0,144.912700, 0 )));
			atomData.insert( VT("Sm" 	, AtomLabelData( 62,3,0,150.360000, 0 )));
			atomData.insert( VT("Eu" 	, AtomLabelData( 63,3,0,151.964000, 0 )));
			atomData.insert( VT("Gd" 	, AtomLabelData( 64,3,0,157.250000, 0 )));
			atomData.insert( VT("Tb" 	, AtomLabelData( 65,3,0,158.925340, 0 )));
			atomData.insert( VT("Dy" 	, AtomLabelData( 66,3,0,162.500000, 0 )));
			atomData.insert( VT("Ho" 	, AtomLabelData( 67,3,0,164.930320, 0 )));
			atomData.insert( VT("Er" 	, AtomLabelData( 68,3,0,167.260000, 0 )));
			atomData.insert( VT("Tm" 	, AtomLabelData( 69,3,0,168.934210, 0 )));
			atomData.insert( VT("Yb" 	, AtomLabelData( 70,3,0,173.040000, 0 )));
			 // Actinide Series
			atomData.insert( VT("Ac" 	, AtomLabelData( 89,3,0,227.027700, 0 )));
			atomData.insert( VT("Th" 	, AtomLabelData( 90,4,0,232.038100, 0 )));
			atomData.insert( VT("Pa" 	, AtomLabelData( 91,5,0,231.035880, 0 )));
			atomData.insert( VT("U" 	, AtomLabelData( 92,6,0,238.028900, 0 )));
			atomData.insert( VT("Np" 	, AtomLabelData( 93,5,0,237.048200, 0 )));
			atomData.insert( VT("Pu" 	, AtomLabelData( 94,4,0,244.064200, 0 )));
			atomData.insert( VT("Am" 	, AtomLabelData( 95,3,0,243.061400, 0 )));
			atomData.insert( VT("Cm" 	, AtomLabelData( 96,3,0,247.070300, 0 )));
			atomData.insert( VT("Bk" 	, AtomLabelData( 97,3,0,247.070300, 0 )));
			atomData.insert( VT("Cf" 	, AtomLabelData( 98,3,0,251.079600, 0 )));
			atomData.insert( VT("Es" 	, AtomLabelData( 99,3,0,252.083000, 0 )));
			atomData.insert( VT("Fm" 	, AtomLabelData(100,3,0,257.095100, 0 )));
			atomData.insert( VT("Md" 	, AtomLabelData(101,3,0,258.098400, 0 )));
			atomData.insert( VT("No" 	, AtomLabelData(102,3,0,259.101100, 0 )));
		}
		return atomData;
	}


//##############################################################################


	const MoleculeUtil::AtomLabelData * const
	MoleculeUtil
	::getAtomData( const std::string& label )
	{
		  // try to locate directly
		AtomDataMap::const_iterator i = getAtomData().find(label);
		  // if found enable pointer access
		if (i!=getAtomData().end()) {
			return &(i->second);
		}
		  // try extraction of atom identifier only
		i = getAtomData().find(getAtom(label));
		  // if found enable pointer access
		if (i!=getAtomData().end()) {
			return &(i->second);
		}

		  // not found ...
		return NULL;
	}


//##############################################################################


	const std::string * const
	MoleculeUtil
	::getAromaticPendant( const std::string& atomLabel )
	{
		  // dummy call to ensure that aromaticSwapMap is initialized
		getAtomData();
		getBondData();

		  // try to locate directly
		AromaticSwapMap::const_iterator i = aromaticSwapMap.find(atomLabel);

		  // if found enable pointer access
		if (i!=aromaticSwapMap.end()) {
			return &(i->second);
		}

		  // not found ...
		return NULL;
	}


//##############################################################################


	std::string
	MoleculeUtil
	::getAtom( const std::string& label )
	{
		  // ensure that no charge or proton information is present
		if ( getAtomData().find(label) != getAtomData().end() ) {
			return label;
		} else {

			  // find first additional information delimiter
			const size_t delim = label.find_first_of("H+-:",1);

			  // check if no delimiter was found -> return no label
			if ( delim == std::string::npos ) {
				return "";
			}

			  // return atom label
			return label.substr(0,delim);
		}
	}


//##############################################################################




	size_t
	MoleculeUtil
	::getProtons( const std::string & label )
	{
		size_t protons = 0;

		  // handle complex atom label
		if ( label.size()>=2 && label.find_first_of("H",1)!=std::string::npos ) {
			const size_t delim = label.find_first_of('H',1);
			  // at least one proton present
			protons = 1;
			  // if something else is present
			if ((delim+1) != label.size()) {
				  // check if a digit is following
				size_t digEnd = label.find_first_not_of(DIGIT,delim+1);
				if ( digEnd != (delim+1) ) {
					  // parse number
					protons = (size_t)atoi( label.substr(delim+1,digEnd-delim+1).c_str() );
					assert( protons > 0 /*otherwise parsing error*/);
				}
			}
		} else if (label == "H" || (label.at(0) == 'H' && label.find_first_of("+-:",1)==1)) {
			protons = 1;
		}

		return protons;
	}

//##############################################################################




	int
	MoleculeUtil
	::getCharge( const std::string & label )
	{
		int charge = 0;

		  // handle complex atom label
		if ( label.size()>=2 && label.find_first_of("+-",1)!=std::string::npos ) {
			const size_t delim = label.find_first_of("+-",1);
			  // at least sign information present
			charge = (label.at(delim) == '+' ? +1 : -1 );
			  // if present
			if (delim != std::string::npos && (delim+1) != label.size()) {
				  // check if a digit is following
				const size_t digEnd = label.find_first_not_of(DIGIT,delim+1);
				if ( digEnd != (delim+1) ) {
					  // parse number
					charge *= atoi( label.substr(delim+1,digEnd-delim+1).c_str() );
				}
			}
		}

		return charge;
	}

//##############################################################################




	int
	MoleculeUtil
	::getClass( const std::string & label )
	{
		int classID = 0;

		  // handle complex atom label
		if ( label.size()>=2 && label.find_first_of(":",1)!=std::string::npos ) {
			const size_t delim = label.find(':',1);
			  // if present
			if (delim != std::string::npos && (delim+1) != label.size()) {
				  // check if a digit is following
				const size_t digEnd = label.find_first_not_of(DIGIT,delim+1);
				if ( digEnd != (delim+1) ) {
					  // parse number
					classID = atoi( label.substr(delim+1,digEnd-delim+1).c_str() );
				}
			}
		}

		return classID;
	}


//##############################################################################



	std::string
	MoleculeUtil
	::getComplexAtomLabel( const std::string& atom
					, const size_t protons
					, const int charge
					, const int classID
					, const bool explicitChargeValue )
	{
		std::stringstream label;

		assert( atom.size() > 0 /*no atom label given*/);
		assert( MoleculeUtil::getAtomData(atom) != NULL /*unknown atom label*/);

		  // set atom label
		label <<atom;

		  // set protons
		if (protons > 0 && atom != "H") {
			label <<"H";
			if (protons > 1) {
				label <<protons;
			}
		}

		  // set charge
		if (charge != 0) {
			int val = abs(charge);
			label <<( charge < 0 ? "-" : "+");
			if (explicitChargeValue || val > 1) {
				label <<val;
			}
		}

		  // set class
		if (classID != 0) {
			label <<":" <<classID;
		}

		return label.str();
	}

//##############################################################################



	const
	MoleculeUtil
	::BondDataMap&
	MoleculeUtil
	::getBondData( void )
	{
		/////////////////////////////////////////////////////////////////
		// NOTE : EVERY SYMBOL ADDED HERE HAS TO BE ADDED TO THE SMILES
		//        PARSER
		/////////////////////////////////////////////////////////////////

		  // check if first call and data has to be set first
		if (bondData.size() == 0) {
			typedef BondDataMap::value_type VT;
			typedef AromaticSwapMap::value_type MT;
			bondData.insert( VT("^" ,  BondLabelData(0,0)));	// new special label for coordinative bonds
			bondData.insert( VT("-" ,  BondLabelData(1,0)));
				aromaticSwapMap.insert( MT("-", ":" ));
			bondData.insert( VT("=" ,  BondLabelData(2,0)));
				aromaticSwapMap.insert( MT("=", ":" ));
			bondData.insert( VT("#" ,  BondLabelData(3,0)));
			bondData.insert( VT("$" ,  BondLabelData(4,0)));
			bondData.insert( VT(":" ,  BondLabelData(1,1)));
//			bondData.insert( VT("/" ,  BondLabelData(1,0)));
//			bondData.insert( VT("\\" , BondLabelData(1,0)));
		}
		return bondData;
	}


//##############################################################################




	const MoleculeUtil::BondLabelData * const
	MoleculeUtil
	::getBondData( const std::string& label )
	{
		  // try to locate
		BondDataMap::const_iterator i = getBondData().find(label);

		  // if found enable pointer access
		if ( i != getBondData().end() ) {
			return &(i->second);
		}

		  // not found ...
		return NULL;
	}


//##############################################################################




	bool
	MoleculeUtil
	::isValidAtomLabel( const std::string & label_ )
	{
		const AtomDataMap& lm = getAtomData();
		  // check if label is found in the available label map
		if ( lm.find( label_ ) != lm.end() ) {

			  // atom is known within registered atom label
			return true;

		/////////////////////////////////////////////
		} else { // atom identifier was NOT found
		/////////////////////////////////////////////
			// not found --> decompose for further checks

			  // ensure no whitespace present
			if (label_.find(WHITESPACE) != std::string::npos) {
//std::cerr <<" whitespace\n";
				return false;
			}

			  // assure complex label
			if (label_.find_first_of("H+-:",1) == std::string::npos) {
//std::cerr <<" no complex label\n";
				return false;
			}

			  // make a copy of the label and attach a delimiter to distinguish if additional information is present or not
			const std::string label = label_ + "]";

			  //////////////////////////////////////////////////////////
			  // find first H that is not the atom label itself
			  //////////////////////////////////////////////////////////
			const size_t hPos = label.find('H',1);
			  // get hydrogen information end
			const size_t hEnd = ( hPos == std::string::npos )
								? hPos
								: label.find_first_not_of(DIGIT,hPos+1);
			if ( hPos != std::string::npos ) {
				  // check for [HH] --> not allowed
				if ( label.substr(0,hPos) == "H" ) {
	//std::cerr <<" HH\n";
					return false;
				}
				  // ensure it is the only additional H position
				if (	(hPos+1) < label.size()
						&& label.find('H',hPos+1) != std::string::npos )
				{
	//std::cerr <<" multiple H\n";
					return false;
				}
				  // check if hEnd >= hPos
				  // and only max ONE digits present
				if (	hEnd != std::string::npos
						&&  (hEnd - hPos) > 2 ) {
	//std::cerr <<" H end < pos\n";
					return false;
				}
			}

			  //////////////////////////////////////////////////////////
			  // find charge information
			  //////////////////////////////////////////////////////////
			const size_t signPos = label.find_first_of("+-",1);
			  // get charge end
			const size_t chargeEnd = ( signPos == std::string::npos )
									? signPos
									: label.find_first_not_of(DIGIT,signPos+1);
			if ( signPos != std::string::npos ) {
				  // ensure signPos > hPos if hPos present
				if (	hPos != std::string::npos
						&& signPos < hPos ) {
	//std::cerr <<" signPos < hPos\n";
					return false;
				}
				  // ensure signPos succeeds hPos, i.e. signPos == hEnd
				if (	hEnd != std::string::npos
						&& signPos != hEnd) {
	//std::cerr <<" signPos != hPos\n";
					return false;
				}
				  // ensure it is the only charge information
				if (	(signPos+1) < label.size()
						&& label.find_first_of("+-",signPos+1) != std::string::npos)
				{
	//std::cerr <<" multiple sign pos\n";
					return false;
				}
				  // check if chargeEnd > signPos
				  // and only max TWO digits present
				if (	chargeEnd != std::string::npos  // <- no digit given
						&&  (chargeEnd - signPos) > 3 ) {
	//std::cerr <<" charge length > 2\n";
					return false;
				}
			}

			  //////////////////////////////////////////////////////////
			  // find class information ":"
			  //////////////////////////////////////////////////////////
			const size_t classPos = label.find(':',1);
			  // get class end
			const size_t classEnd = ( classPos == std::string::npos )
									? classPos
									: label.find_first_not_of(DIGIT,classPos+1);
			if ( classPos != std::string::npos ) {
				  // ensure classPos > hPos
				if (	hPos != std::string::npos
						&& classPos < hPos ) {
	//std::cerr <<" class < hPos\n";
					return false;
				}
				  // ensure classPos succeeds hPos if no charge present,
				  // i.e. classPos == hEnd
				if (	signPos == std::string::npos
						&& hEnd != std::string::npos
						&& classPos != hEnd) {
	//std::cerr <<" class == hEnd if no charge \n";
					return false;
				}
				  // ensure classPos > signPos
				if (	signPos != std::string::npos
						&& classPos < signPos ) {
	//std::cerr <<" class < signPos \n";
					return false;
				}
				  // ensure classPos succeeds signPos, i.e. classPos == signEnd
				if (	chargeEnd != std::string::npos
						&& classPos != chargeEnd) {
	//std::cerr <<" classPos != chargeEnd if charge present\n";
					return false;
				}
				  // check if chargeEnd > classPos
				if (	classEnd != std::string::npos
						&& classEnd == (classPos+1) ) {
	//std::cerr <<" classPos == classEnd \n";
					return false;
				}
			}

			  //////////////////////////////////////////////////////////
			  // check atom label
			  //////////////////////////////////////////////////////////
			size_t atomEnd = std::min( hPos, std::min( signPos, classPos ) );
			  // if no H, sign, charge information found --> label is invalid
			if (atomEnd == std::string::npos) {
//std::cerr <<" atomEnd == npos \n";
//std::cerr <<" atom label unknown but no H, sign, charge information found
				return false;
			}
			assert(atomEnd != std::string::npos);

			  // try to find atom label without brackets
			if ( lm.find( label.substr(0,atomEnd) ) == lm.end() ) {
//std::cerr <<" label " <<label.substr(1,atomEnd-1) <<" not found \n";
				return false;
			}

			  //////////////////////////////////////////////////////////
			  // ensure nothing strange left at the end
			  //////////////////////////////////////////////////////////
			  // get end of handled string
			size_t end = (hEnd != std::string::npos) ? hEnd : 0;
			end = std::max( end, ((chargeEnd != std::string::npos) ? chargeEnd : 0));
			end = std::max( end, ((classEnd != std::string::npos) ? classEnd : 0));
			if ( end == 0 ) {
				end = atomEnd;
			}
			  // ensure end is final
			if ( end != (label.size()-1) ) {
//std::cerr <<" end " <<end <<" != label.size()-2 \n";
				return false;
			}

			  // else it is a valid complex atom label
			return true;

		}
		return true;
	}


//##############################################################################




	size_t
	MoleculeUtil
	::isConsistent( const Molecule & mol )
	{

		size_t retVal = 1;

		bool atomLabelsOK = true, bondLabelsOK = true;

		  // perform atom label tests
		if ( ! checkAtomLabel(mol) ) {
			retVal *= C_AtomLabelInvalid;
			atomLabelsOK = false;
		}

		if ( atomLabelsOK && ! checkAtomComplexWithH(mol) ) {
			retVal *= C_AtomComplexWithH;
			atomLabelsOK = false;
		}

		if ( ! checkAtomLabelWildcard(mol) ) {
			retVal *= C_AtomLabelWildcard;
			atomLabelsOK = false;
		}

		if ( ! checkBondLabel(mol) ) {
			retVal *= C_BondLabelInvalid;
			bondLabelsOK = false;
		}

		if ( atomLabelsOK && bondLabelsOK && ! checkAtomValence(mol) ) {
			retVal *= C_AtomValence;
		}

		if ( ! checkBondLoop(mol) ) {
			retVal *= C_BondLoop;
		}

		if ( ! checkNonConnected(mol) ) {
			retVal *= C_NonConnected;
		}

		if (retVal == 1) {
			retVal = C_Consistent;
		}

		return retVal;
	}


//##############################################################################




	bool
	MoleculeUtil
	::decodeConsistencyStatus(	const size_t consistencyCode
								, std::ostream& errorStream )
	{
		if (consistencyCode == C_Consistent)
			return true;

		  // temporary error code handling
		size_t errorCode = consistencyCode;
		errorStream	<<"The molecule is not consistent due to : \n";
		if ( errorCode % MoleculeUtil::C_AtomLabelInvalid == 0) {
			errorStream <<" + not supported atom labels\n";
			errorCode /= MoleculeUtil::C_AtomLabelInvalid;
		}
		if ( errorCode % MoleculeUtil::C_BondLabelInvalid == 0) {
			errorStream <<" + not supported bond labels\n";
			errorCode /= MoleculeUtil::C_BondLabelInvalid;
		}
		if ( errorCode % MoleculeUtil::C_AtomLabelWildcard == 0) {
			errorStream <<" + wildcards in atom labels (currently only allowed in rules not in molecules)\n";
			errorCode /= MoleculeUtil::C_AtomLabelWildcard;
		}
		if ( errorCode % MoleculeUtil::C_AtomComplexWithH == 0) {
			errorStream <<" + implicit protons within complex atom labels (currently not supported in the framework, make explicit)\n";
			errorCode /= MoleculeUtil::C_AtomComplexWithH;
		}
		if ( errorCode % MoleculeUtil::C_AtomValence == 0) {
			errorStream <<" + electron distribution of an atom is not fulfilling 0 <= [(bondValenceSum+aromAdd) - (valence+charge)] <= aromAdd, where aromAdd == 1 if the atom is part of an aromatic ring and 0 otherwise\n";
			errorCode /= MoleculeUtil::C_AtomValence;
		}
		if ( errorCode % MoleculeUtil::C_BondLoop == 0) {
			errorStream <<" + an atom forms a bond with itself (bond loop)\n";
			errorCode /= MoleculeUtil::C_BondLoop;
		}
		if ( errorCode % MoleculeUtil::C_NonConnected == 0) {
			errorStream <<" + not connected, i.e. more than one connected component within molecule graph\n";
			errorCode /= MoleculeUtil::C_NonConnected;
		}
		if ( errorCode != 1 ) {
			errorStream <<" + error code "<<errorCode<<" has no description, sorry...\n" <<std::endl;
		}
		return false;
	}



//##############################################################################




	bool
	MoleculeUtil
	::checkAtomLabel
	( const Molecule & mol )
	{
		Molecule::vertex_iterator vIt, vItEnd;
		boost::property_map<	Molecule, PropNodeLabel >::const_type
			nodeLabel = boost::get( PropNodeLabel(), mol );

		  // check all atom nodes
		for (boost::tie(vIt,vItEnd) = boost::vertices(mol); vIt != vItEnd; ++vIt) {
			  // check atom label
			if (!isValidAtomLabel( nodeLabel[*vIt] )) {
				return false;
			}

		}

		return true;
	}


//##############################################################################




	bool
	MoleculeUtil
	::checkAtomLabelWildcard
	( const Molecule & mol )
	{
		Molecule::vertex_iterator vIt, vItEnd;
		boost::property_map<	Molecule, PropNodeLabel >::const_type
			nodeLabel = boost::get( PropNodeLabel(), mol );

		  // check all atom nodes
		for (boost::tie(vIt,vItEnd) = boost::vertices(mol); vIt != vItEnd; ++vIt) {
			  // check if atom label is a wildcard
			if (nodeLabel[*vIt].compare(AtomLabelWildcard) == 0) {
				return false;
			}
		}

		return true;
	}


//##############################################################################




	bool
	MoleculeUtil
	::checkAtomComplexWithH
	( const Molecule & mol )
	{
		Molecule::vertex_iterator vIt, vItEnd;
		boost::property_map<	Molecule, PropNodeLabel >::const_type
			nodeLabel = boost::get( PropNodeLabel(), mol );

		  // check all atom nodes
		for (boost::tie(vIt,vItEnd) = boost::vertices(mol); vIt != vItEnd; ++vIt) {
			  // access shortcut
			const std::string& curNodeLabel = nodeLabel[*vIt];
			  // check if the node shows implicit protons
			if (	getAtom(curNodeLabel).compare("H") != 0
					&& getProtons(curNodeLabel) > 0 )
			{
				return false;
			}
		}

		return true;
	}


//##############################################################################




	bool
	MoleculeUtil
	::checkAtomValence
	( const Molecule & mol )
	{
		Molecule::vertex_iterator vIt, vItEnd;
		Molecule::out_edge_iterator ae, end;
		boost::property_map<	Molecule, PropNodeLabel >::const_type
			nodeLabel = boost::get( PropNodeLabel(), mol );
		boost::property_map<	Molecule, PropEdgeLabel >::const_type
			edgeLabel = boost::get( PropEdgeLabel(), mol );

		  // check all atom nodes
		for (boost::tie(vIt,vItEnd) = boost::vertices(mol); vIt != vItEnd; ++vIt) {

			  // access shortcut
			const std::string& curNodeLabel = nodeLabel[*vIt];

			  // get atom data
			const AtomLabelData * const atomData = getAtomData(curNodeLabel);
			assert( atomData != NULL );

			  // check if we have to do a valence check for this type of atom
			if (atomData->isToBeChecked == 0) {
				continue;
			}

			  // corrected charge value for positively charged atoms
			const int curCharge = getCharge( curNodeLabel ) * (getAtom(curNodeLabel).compare("H") == 0 ? -1 : +1);

			int bondValenceSum = 0;
			int aromAdd = 0;  // can be either 0 or 1
			for (boost::tie(ae,end) = out_edges( *vIt, mol ); ae!=end; ++ae) {
				const BondLabelData * const bondData = getBondData( edgeLabel[*ae] );
				assert( bondData != NULL);
				  // if aromatic edge adjacent -> increase counter and stop
				if (aromAdd == 0 && bondData->isAromatic != 0) {
					aromAdd++;
				}
				  // count edge valence
				bondValenceSum += (int)bondData->valence;
			}

			  // compute final valence to test for
			const int curValence = (bondValenceSum) - (curCharge);


			  // get available valences
			const AtomLabelData::OneByteVec & valences = atomData->valenceAlternatives;
			  // screen for a fitting one
			bool noFittingValence = true;
			for (size_t i=0; i<valences.size() && noFittingValence; ++i) {

				  // do next valence test with and without aromaticity factor
				noFittingValence = (int)valences.at(i) != curValence
								&& (int)valences.at(i) != (curValence+aromAdd);

//				  // check if the electron distribution is sane
//				int difference = (bondValenceSum + aromAdd) - (curValence + curCharge);
//				  // i.e. (0 <= difference <= aromAdd)
//				if ( difference < 0 || aromAdd < difference )
//				{
//					return false;
//				}
			}
			if (noFittingValence) {

//#ifndef NDEBUG
////					std::stringstream msg;
//					std::ostream & msg = std::cerr;
//					msg <<"unbalanced atom '" <<curNodeLabel <<"'"
//						 <<" valences ";
//					for (size_t i=0; i<valences.size() && noFittingValence; ++i) {
//						msg <<(int)valences.at(i) <<" ";
//					}
//					msg	 <<": charge " << curCharge
////						 <<", proton count " << bondNumProton
//						 <<", bond valences " << bondValenceSum
//						 <<", aromaticEdges " << aromAdd
//						 <<" : curValence = " <<curValence
//						 ;
////					throw std::runtime_error(msg.str());
//#endif

				return false;
			}
		}

		return true;
	}


//##############################################################################




	bool
	MoleculeUtil
	::checkBondLabel
	( const Molecule & mol )
	{
		Molecule::edge_iterator eIt, eItEnd;
		boost::property_map<	Molecule, PropEdgeLabel >::const_type
			edgeLabel = boost::get( PropEdgeLabel(), mol );

		  // check all bonds
		for (boost::tie(eIt,eItEnd) = boost::edges(mol); eIt != eItEnd; ++eIt) {
			  // check edge label
			if (!isValidBondLabel( edgeLabel[*eIt] )) {
				return false;
			}
		}

		return true;
	}


//##############################################################################




	bool
	MoleculeUtil
	::checkBondLoop
	( const Molecule & mol )
	{
		Molecule::edge_iterator eIt, eItEnd;

		  // check all bonds
		for (boost::tie(eIt,eItEnd) = boost::edges(mol); eIt != eItEnd; ++eIt) {
			  // check if ring bond, ie. connects an atom with itself
			if (boost::source(*eIt,mol) == boost::target(*eIt,mol)) {
				return false;
			}
		}

		return true;
	}


//##############################################################################


	void
	MoleculeUtil
	::compressHnodes( Molecule & mol )
	{

		// vertex and edge iterator types
		boost::graph_traits<Molecule>::vertex_iterator     vi, vi_end;
		boost::graph_traits<Molecule>::adjacency_iterator  i, i_end;
		typedef boost::graph_traits<Molecule>::vertex_descriptor M_t;
		std::map< M_t, size_t > hCount;
		boost::property_map<Molecule, PropNodeLabel>::type vname = boost::get(PropNodeLabel(), mol);

		  // access to label
		  // compress labels (remove H-nodes)
		for(boost::tie(vi, vi_end)=boost::vertices(mol); vi!=vi_end; ++vi) {
			if (vname[*vi] == "H") {
				  // check if degree is one
				if (boost::out_degree(*vi,mol) != 1) {
					continue;
				}
				  // update H-count of neighbored nodes
				for (boost::tie(i, i_end) = boost::adjacent_vertices(*vi, mol); i != i_end; ++i){
					  // ensure it is no H-H connection
					if (vname[*i] == "H") {
						continue;
					}
					  // update H-count
					if (hCount.find(*i) == hCount.end()) {
						hCount[*i] = 1;
					} else {
						hCount[*i] += 1;
					}
				}
			}
		}
		  // update all node labels with adjacent H-nodes
		for (std::map< M_t, size_t >::const_iterator u = hCount.begin(); u!=hCount.end(); ++u) {
			std::string name = vname[u->first];
			assert(name.size() > 0 /* no node label information available */);
			  //  check if H available
			assert(name.find('H') == std::string::npos /* H information already available.. update necessary but not implemented */);
			  // find insertion position if any
			const size_t insertPos = name.find_first_of("+-:");
			  // compile new label
			std::stringstream sstr;
			  // add atom label
			sstr << name.substr(0,insertPos);
			  // add H information
			sstr <<'H';
			if (u->second > 1) {
				sstr <<u->second;
			}
			  // add remaining label
			if (insertPos != std::string::npos) {
				 // check if there was already a complex label ending with ']'
				sstr <<name.substr(insertPos,(name.size()-insertPos));
			}
			  // update node label
			vname[u->first] = sstr.str();
		}
		  // remove H-nodes
		bool tryToRemove = (hCount.size() > 0);
		while (tryToRemove) {
			tryToRemove = false;
			for(boost::tie(vi, vi_end)=boost::vertices(mol); vi!=vi_end; ++vi) {
				if (vname[*vi] == "H") {
					  // check if degree is one
					if (boost::out_degree(*vi,mol) != 1) {
						continue;
					}
					  // ensure it is no H-H connection
					if (vname[*(boost::adjacent_vertices(*vi, mol).first)] == "H") {
						continue;
					}
					  // update H-count
					boost::clear_vertex(*vi,mol);
					boost::remove_vertex(*vi,mol);
					tryToRemove = true;
					break;
				}
			}
		}

	}


//##############################################################################

//##############################################################################


	void
	MoleculeUtil
	::removeProtons( Molecule & mol )
	{

		// vertex and edge iterator types
		boost::graph_traits<Molecule>::vertex_iterator     vi, vi_end;
		boost::graph_traits<Molecule>::adjacency_iterator  i, i_end;
		boost::property_map<Molecule, PropNodeLabel>::type vname = boost::get(PropNodeLabel(), mol);
		boost::property_map<Molecule, PropEdgeLabel>::type ename = boost::get(PropEdgeLabel(), mol);

		const std::string protonInfo = "§H§";

		  // remove all H-nodes
		bool oneNodeRemoved = false;
		do {
			oneNodeRemoved = false;
			for(boost::tie(vi, vi_end)=boost::vertices(mol); vi!=vi_end; ++vi) {
				if (vname[*vi] == "H") {
					  // check if degree is one
					if (boost::out_degree(*vi,mol) != 1) {
						continue;
					}
					  // ensure it is no H-H connection
					if (vname[*(boost::adjacent_vertices(*vi, mol).first)] == "H") {
						continue;
					} else {
						// add proton removal information to adjacent atom label
						std::string oldLabel = vname[*(boost::adjacent_vertices(*vi, mol).first)];
						vname[*(boost::adjacent_vertices(*vi, mol).first)] = oldLabel + protonInfo;
					}
					  // update graph
					boost::clear_vertex(*vi,mol);
					boost::remove_vertex(*vi,mol);
					oneNodeRemoved = true;
					break;
				}
			}
		} while (oneNodeRemoved);

		  // remove all additional proton information from atom labels
		for(boost::tie(vi, vi_end)=boost::vertices(mol); vi!=vi_end; ++vi) {
			  // check non-proton labels for implicit proton information
			if (vname[*vi] != "H") {
				std::string oldAtomLabel = vname[*vi];
				int adjacentProtons = 0;
				  // check for the number of protons removed
				if (oldAtomLabel.find(protonInfo) != std::string::npos) {
					  // get number of removed protons
					const size_t infoPos = oldAtomLabel.find(protonInfo);
					adjacentProtons += ((int)(oldAtomLabel.size()-infoPos) / (int)protonInfo.size());
					  // remove proton removal information
					oldAtomLabel = oldAtomLabel.substr(0,infoPos);
					vname[*vi] = oldAtomLabel;
				}
				  // ignore group labels
				if (isGroupLabel(vname[*vi])) {
					continue;
				}

				  // get number of implicit protons within label
				int implicitProtons = (int)getProtons(oldAtomLabel);

				  // check if any protons there -> if not ignore
				if ((implicitProtons+adjacentProtons) == 0) {
					continue;
				}

				  // check for adjacent aromatic edges
				Molecule::out_edge_iterator ae, end;
				int aromAdd = 0;
				int bondValences = 0;
				const BondLabelData * bondData = NULL;
				for (boost::tie(ae,end) = out_edges( *vi, mol ); ae!=end; ++ae) {
					bondData = getBondData( ename[*ae] );
					assert( bondData != NULL);
					  // if aromatic edge adjacent -> increase counter and stop
					if (aromAdd == 0 && bondData->isAromatic != 0) {
						aromAdd++;
					}
					  // count edge valence
					bondValences += (int)bondData->valence;
				}

				  // get access to atom specific data
				const AtomLabelData * atomData = getAtomData( oldAtomLabel );
				assert(atomData != NULL);

				  // check if protons can be savely removed without loosing information
				  // ie. proton number can be inferred later
				if ((adjacentProtons + implicitProtons + bondValences + aromAdd) == ( atomData->valence + getCharge( oldAtomLabel ))) {
					  // check if we have to remove implicit proton information from the atom label
					if (implicitProtons > 0) {
						// additional proton (H) information present that has to be removed
						const size_t posH = oldAtomLabel.find("H",1);
						const size_t endH = oldAtomLabel.find_first_not_of("H1234567890",posH);
						  // check if nothing else within label behind proton information
						if (endH == std::string::npos) {
							// replace atom label with pure atom name without H-information
							vname[*vi] = oldAtomLabel.substr(0,posH);
						} else {
							// replace atom label with anything but the H-information
							vname[*vi] = oldAtomLabel.substr(0,posH) + oldAtomLabel.substr(endH);
						}
					}
				} else {
					  // proton information is important and cannot be inferred
					  // thus, we have to create an according complex atom label
					vname[*vi] = getComplexAtomLabel( getAtom(oldAtomLabel)
							, (adjacentProtons+implicitProtons)
							, getCharge( oldAtomLabel )
							, getClass( oldAtomLabel) );
				}
			}
		}

	}



//##############################################################################

	size_t
	MoleculeUtil
	::getProtonsToAdd( const std::string & atomLabel
					, const int atomCharge
					, const size_t bondValenceSum
					, const size_t bondNum
					, const size_t bondNumAromatic
					, const size_t bondNumProton
					)
	{

		  // get atom information
		const AtomLabelData * atomData = getAtomData( atomLabel );
		assert(atomData != NULL /*atom label unknown*/);

		  // get overall allowed valence
		int remainingValence = (int)atomData->valence;

		  // charge
		remainingValence += (atomLabel=="H"?-1:+1) * atomCharge;

		  // bonds
		remainingValence -= (int)bondValenceSum;

		  // handle aromatic edges
		remainingValence -= (int)floor( (double)bondNumAromatic * 0.5 );

//		// special handlings
//		if (atomData->atomicNumber == 1) { // "H"
//
//		} else if (atomData->atomicNumber == 7) { // "N"
//
////			if ((bondNum) == 3 && bondValenceSum == 4 && charge != 1) {
////				throw std::runtime_error("'N' atom with 3 bonds and bond valence sum of 4 should be positively charged!");
////			}
////
////			if (remainingValence < 0) {
////				if (bondNum == 3 && bondValenceSum == 5) {
////					remainingValence = 0;
////				}
////				if (bondNum == 5 && bondValenceSum == 7 && bondNumProton == 0) {
////					remainingValence = 0;
////				}
////			}
//
//		} else if (atomData->atomicNumber == 15) { // "P"
//
////			if ((bondNum) == 3 && bondValenceSum == 4 && charge != 1) {
////				throw std::runtime_error("'P' atom with 3 bonds and bond valence sum of 4 should be positively charged!");
////			}
//
//		} else if (atomData->atomicNumber == 16) { // "S"
//
//			  // set upper bound on attached protons
//			if (remainingValence > 0 && remainingValence > (int)atomData->maxProtons) {
//				remainingValence = (int)atomData->maxProtons;
//			}
//		}

		  // check for unbalanced valence
//#ifndef NDEBUG
//		if (atomData->atomicNumber<=20 && remainingValence < 0) {
//			std::stringstream msg;
//			msg <<"unbalanced atom '" <<atomLabel <<"'"
//				 <<" valence " <<(int)atomData->valence
//				 <<", charge " << atomCharge
//				 <<", proton count " << bondNumProton
//				 <<", bond valences " << bondValenceSum
//				 <<", aromaticEdges " << bondNumAromatic
//				 <<" : remaining valence = " <<remainingValence
//				 ;
//			throw std::runtime_error(msg.str());
//			assert( remainingValence >= 0 /*unbalanced atom*/);
//		}
//#endif

//		  // check if an explicit maximal proton number was given
		 // CURRENTLY ONLY FOR "S" ENABLED -> THUS HANDLED ABOVE
//		if (atomData->valence != atomData->maxProtons) {
//			remainingValence += (int)atomData->maxProtons - (int)atomData->valence;
//			assert( (int)atomData->maxProtons >= ((int)bondNumProton - atomCharge) /*too many protons attached*/ );
//		}


		  // truncate at 0
		return remainingValence >= 0 ? (size_t)remainingValence : 0;
	}

//##############################################################################


	void
	MoleculeUtil
	::fillProtons( Molecule & mol )
	{
		  // check for unbalanced valence
#ifndef NDEBUG
		const std::string initialSMILES = SMILESwriter::getSMILES(mol, true, true);
#endif

		// vertex and edge iterator types
		boost::graph_traits<Molecule>::vertex_iterator     vi, vi_end;
		boost::graph_traits<Molecule>::out_edge_iterator   ei, ei_end;
		boost::graph_traits<Molecule>::vertex_descriptor addedVertex;
		boost::graph_traits<Molecule>::edge_descriptor addedEdge;
		boost::property_map<Molecule, PropNodeLabel>::type vname = boost::get(PropNodeLabel(), mol);
		boost::property_map<Molecule, PropEdgeLabel>::type ename = boost::get(PropEdgeLabel(), mol);

		  // check all nodes
		for(boost::tie(vi,vi_end)=boost::vertices(mol); vi!=vi_end; ++vi) {
			  // get atom information
			const AtomLabelData * atomData = getAtomData( getAtom(vname[*vi]) );
			assert(atomData != NULL /*atom label unknown*/);

			  // get overall allowed valence
			int remainingValence = (int)atomData->valence;

			  // take charge into account (take care if current atom is a proton!)
			const int curCharge = (int)getCharge( vname[*vi] );
//			remainingValence += (getAtom(vname[*vi])=="H"?-1:+1) * curCharge;

			  // get implicit protons (take care if current atom is a proton!)
			const int curImplicitProtons = getProtons(vname[*vi]) - (getAtom(vname[*vi])=="H"?1:0);
			  // relabel vertex if necessary
			if (curImplicitProtons > 0) {
				vname[*vi] = getComplexAtomLabel( getAtom(vname[*vi])
									, 0  // no implicit protons anymore
									, getCharge(vname[*vi])
									, getClass(vname[*vi])
									, false );
			}
			  // create explicit protons
			for (int p=0; p < curImplicitProtons; ++p) {
				  // add proton
				addedVertex = boost::add_vertex( mol );
				vname[addedVertex] = "H";
				  // make adjacent
				addedEdge = boost::add_edge( *vi, addedVertex, mol ).first;
				ename[addedEdge] = "-";
			}

			  // get bond valence reduction
			boost::tie( ei, ei_end ) = boost::out_edges(*vi, mol);
			int aromaticEdges = 0;
			int protonCount = 0;
			int bondValences = 0;
			int bondNumber = 0;

			for (; ei!=ei_end; ++ei) {
				const BondLabelData * bondData = getBondData( ename[*ei] );
				assert(bondData != NULL /*bond label unknown*/);

//				remainingValence -= (int)bondData->valence;
				bondValences += (int)bondData->valence;
				bondNumber++;

				if (bondData->isAromatic != 0) {
					aromaticEdges++;
				}
				protonCount += getAtom(vname[boost::target(*ei,mol)])=="H" ? 1 : 0;
			}

			try {
				  // get number of protons still to add
				remainingValence = (int)getProtonsToAdd( getAtom(vname[*vi]), curCharge, bondValences, bondNumber, aromaticEdges, protonCount);

			} catch (std::exception & ex) {
#ifndef NDEBUG
				throw std::runtime_error("fillProtons( '"+initialSMILES+"' )\n\t'"
				+SMILESwriter::getSMILES(mol, false, true)+"' : \n\t"+
#else
				throw std::runtime_error(
#endif
						ex.what());
			}


//			  // handle aromatic edges
//			remainingValence -= (int)floor( aromaticEdges * 0.5 );
//
//			  // check if an explicit maximal proton number was given
//			if (atomData->valence != atomData->maxProtons) {
//				remainingValence += (int)atomData->maxProtons - (int)atomData->valence;
//				assert( (int)atomData->maxProtons >= (protonCount - curCharge) /*too many protons attached*/ );
//			}

			  // add protons
			for (; remainingValence > 0; remainingValence--) {
				  // add proton
				addedVertex = boost::add_vertex( mol );
				vname[addedVertex] = "H";
				  // make adjacent
				addedEdge = boost::add_edge( *vi, addedVertex, mol ).first;
				ename[addedEdge] = "-";
			}
		}
	}


//##############################################################################


	std::ostream&
	MoleculeUtil
	::convertCML( const Molecule& m, std::ostream& molCML )
	{
		  // copy to enable label compression
		Molecule m2;
		copy( m, m2 );
		  // compress protons into the node labels
		compressHnodes( m2 );

		  // constant access
		const Molecule& mol = m2;
		boost::property_map< Molecule, PropNodeLabel>::const_type
			mNodeLabel = (boost::get( PropNodeLabel(), mol ));
		boost::property_map< Molecule, PropEdgeLabel>::const_type
			mEdgeLabel = (boost::get( PropEdgeLabel(), mol ));
		boost::property_map< Molecule, PropNodeIndex>::const_type
			mNodeIndex = (boost::get( PropNodeIndex(), mol ));

		Molecule::vertex_iterator vi, v_end;
		Molecule::edge_iterator ei, e_end;

		  // gather atom data
		molCML <<"<molecule>\n"
				<<"<atomArray>\n";
		for (boost::tie(vi,v_end) = boost::vertices(mol); vi!=v_end; ++vi) {
			molCML <<" <atom id=\"a"<<mNodeIndex[*vi]
					<<"\" elementType=\""<<getAtom(mNodeLabel[*vi]);
			int charge = getCharge(mNodeLabel[*vi]);
			if (charge != 0) {
				molCML <<"\" charge=\""<<charge;
			}
			size_t hydrogenCount = getProtons(mNodeLabel[*vi]);
			if (hydrogenCount > 0) {
				molCML <<"\" hydrogenCount=\""<<hydrogenCount;
			}
			molCML <<"\"/>\n";
		}
		molCML <<"</atomArray>\n"
				<<"<bondArray>\n";
		for (boost::tie(ei,e_end) = boost::edges(mol); ei!=e_end; ++ei) {
			const BondLabelData* data = getBondData(mEdgeLabel[*ei]);
			molCML <<" <bond atomRefs2=\"a"<<mNodeIndex[boost::source(*ei,mol)]
					<<" a"<<mNodeIndex[boost::target(*ei,mol)]
					<<"\" order=\"";
			if (data->isAromatic==1) {
				molCML <<"A";
			} else {
				molCML <<((size_t)data->valence);
			}
			molCML <<"\"/>\n";

		}
		molCML <<"</bondArray>\n"
				<<"</molecule>\n";

		  // return string CML representation
		return molCML;
	}


//##############################################################################

	void
	MoleculeUtil
	::copy( const sgm::Graph_Interface& graph, Molecule & mol)
	{
		assert( boost::num_vertices(mol) == 0 );

		// vertex and edge iterator types
		boost::graph_traits<Molecule>::vertex_descriptor addedVertex;
		boost::graph_traits<Molecule>::edge_descriptor addedEdge;
		boost::property_map<Molecule, PropNodeLabel>::type vname = boost::get(PropNodeLabel(), mol);
		boost::property_map<Molecule, PropEdgeLabel>::type ename = boost::get(PropEdgeLabel(), mol);

		  // add all atoms
		const size_t atomNumber = graph.getNodeNumber();
		for (size_t i=0; i<atomNumber; ++i) {
			  // add atom to molecule
			addedVertex = boost::add_vertex( mol );
			  // set atom label
			vname[addedVertex] = graph.getNodeLabel(i);
		}

		  // copy edges
		sgm::Graph_Interface::OutEdge_iterator edge, edgesEnd;
		for (size_t i=0; i<atomNumber; ++i) {
			for(edge = graph.getOutEdgesBegin( i ),edgesEnd = graph.getOutEdgesEnd( i ); edge != edgesEnd; ++edge)
			{
				  // handle only ordered edges to avoid multiple edge insertions
				if (edge->getFromIndex() < edge->getToIndex()) {
					  // add bond
					addedEdge = boost::add_edge( boost::vertex( edge->getFromIndex(), mol )
												,boost::vertex( edge->getToIndex(), mol )
												,mol ).first;
					  // set bond label
					ename[addedEdge] = edge->getEdgeLabel();
				}
			}

		}

	}

//##############################################################################

	void
	MoleculeUtil::
	insertGroups( Molecule & mol, const GroupMap & groups ) throw(std::runtime_error)
	{

		using namespace boost;

		  // access to the property maps to fill
		boost::property_map< Molecule, PropNodeLabel >::type nodeLabel = get( PropNodeLabel(), mol);
		boost::property_map< Molecule, PropEdgeLabel >::type edgeLabel = get( PropEdgeLabel(), mol);
		Molecule::vertex_descriptor atom, newAtom;
		Molecule::edge_descriptor newEdge;

		  // iterate over all original molecule nodes (excluding new components)
		  // (avoids recursive group replacement)
		const size_t oldMolSize = boost::num_vertices(mol);
		for ( size_t atomID = 0; atomID < oldMolSize; ++atomID ) {

			  // access node
			atom = boost::vertex( atomID, mol );

			std::string &atomLabel = nodeLabel[atom];

			  // check if node is a group label proxy node
			if (! isGroupLabel(atomLabel) ) {
				continue;
			}

			  // derive groupID
			std::string groupID = atomLabel.substr(0,atomLabel.find('}')+1);
			  // check if groupID known
			if (groups.find(groupID) == groups.end()) {
				std::ostringstream oss;
				oss	<<"ggl::chem::MoleculeUtil::insertGroup : the node " <<atomID
					<<" shows the unknown group label '"<<groupID
					<<"' within its node label '"<<atomLabel<<"'!";
				throw std::runtime_error(oss.str());
			}

			  // access according group (presence ensured with test from above)
			const MoleculeComponent &group = groups.find(groupID)->second;
			const size_t groupProxy = *group.compIDs.begin();
			  // group graph access
			boost::property_map< MoleculeComponent::PatternGraph, PropNodeLabel >::const_type groupNodeLabel = get( PropNodeLabel(), group.pattern);
			boost::property_map< MoleculeComponent::PatternGraph, PropNodeIndex >::const_type groupNodeIndex = get( PropNodeIndex(), group.pattern);
			boost::property_map< MoleculeComponent::PatternGraph, PropEdgeLabel >::const_type groupEdgeLabel = get( PropEdgeLabel(), group.pattern);


			// relabel compID node (preserve charge/protons = rest of original atom label)
			nodeLabel[atom] = groupNodeLabel[boost::vertex(groupProxy, group.pattern)]
							 + atomLabel.substr(groupID.size(),std::string::npos);

			// insert remaining atoms/bonds of molecule component

			  // mapping of component atom indices to inserted atoms within mol
			std::map< MoleculeComponent::PatternGraph::vertex_descriptor
					, Molecule::vertex_descriptor> group2mol;
			  // proxy node exists already -> mapping known
			group2mol[boost::vertex(*group.compIDs.begin(), group.pattern)] = atom;

			  // add and set all group atoms (excluding proxy since already there)
			MoleculeComponent::PatternGraph::vertex_iterator  vi, vi_end;
			for(boost::tie(vi, vi_end)=boost::vertices(group.pattern); vi!=vi_end; ++vi) {
				if (groupNodeIndex[*vi] != groupProxy) {
					  // add new node
					newAtom = boost::add_vertex( mol );
					  // store mapping
					group2mol[*vi] = newAtom;
					  // set label
					nodeLabel[newAtom] = groupNodeLabel[*vi];
				}
			}
			  // add and set all group bonds
			MoleculeComponent::PatternGraph::edge_iterator  ei, ei_end;
			for(boost::tie(ei, ei_end)=boost::edges(group.pattern); ei!=ei_end; ++ei) {
				  // add new edge
				newEdge = boost::add_edge( group2mol[boost::source(*ei,group.pattern)]
										 , group2mol[boost::target(*ei,group.pattern)]
										 , mol).first;
				  // set label of new edge
				edgeLabel[newEdge] = groupEdgeLabel[*ei];
			}
		}


	}

//##############################################################################

//		struct CompareGT {
//			bool operator()(const std::pair<size_t, std::string> & p1, const std::pair<size_t, std::string> & p2 ) const {
//				return p1.first > p2.first
//						|| (p1.first == p2.first && p1.second < p2.second);
//			}
//		} compare;
	void
	MoleculeUtil::
	compressGroups( ggl::chem::Molecule & mol
				, const GroupMap & groups
				)
	{
		boost::property_map< Molecule, PropNodeLabel >::type nodeLabel = get( PropNodeLabel(), mol);
		Molecule::vertex_descriptor atom;

		  // create list of groups in decreasing size
		typedef std::pair<int, std::string> GT;
		std::vector< GT > orderedGroups( groups.size() );
		size_t i = 0;
		for ( GroupMap::const_iterator g = groups.begin(); g != groups.end(); ++g,++i) {
			orderedGroups[i] = GT( - boost::num_vertices( g->second.pattern ), g->first );
			assert(g->second.compIDs.size() == 1);
		}

		  // sort list with decreasing group size
		std::sort(orderedGroups.begin(),orderedGroups.end());

		  // check each group if present
		for (std::vector< GT >::const_iterator g = orderedGroups.begin(); g != orderedGroups.end(); ++g ) {

			  // skip all groups that are larger than the remaining molecule
			if (-(g->first) > (int)boost::num_vertices(mol) ) {
				continue;
			}

			  // get group
			const MoleculeComponent & group = groups.find(g->second)->second;

			// create list of degree constraints for all but the interface nodes
			sgm::Pattern_Interface::ConstraintVec constraints( boost::num_vertices(group.pattern) - group.compIDs.size() );
			for (size_t i = 0, n=0; n < boost::num_vertices(group.pattern); ++n) {
				  // ensure this is no compID node
				if (group.compIDs.find(n) == group.compIDs.end()) {
					 // create degree constraint
					constraints[i] = new sgm::MC_NodeAdjacency( n
							, sgm::MC_NodeAdjacency::MC_EQ
							, boost::degree( boost::vertex( n, group.pattern), group.pattern )
							, AtomLabelWildcard);
					++i;
				}
			}

			// create pattern
			sgm::Graph_boost<MoleculeComponent::PatternGraph> patternGraph(group.pattern);
			sgm::Pattern searchPattern( patternGraph, constraints, AtomLabelWildcard );

			// iterate find matches and compress until no matching
			bool matchFound = false;
			sgm::SGM_vf2 vf2;
			do {
				matchFound = false;
				  // skip all groups that are larger than the remaining molecule
				if (-(g->first) > (int)boost::num_vertices(mol) ) {
					break;
				}
				Molecule_Graph molGraph(mol);
				  // setup match storing
				std::vector< sgm::Match > match;
				sgm::MR_Storing matchReporter(match);
				  // run matching
				vf2.findMatches( searchPattern, molGraph, matchReporter, 1 );
				  // handle match
				if (!match.empty()) {
					matchFound = true;
					  // collect all nodes to remove and rename hub node
					std::set< size_t > toRemove;
					for (size_t i=0; i<match.begin()->size(); ++i) {
						if (group.compIDs.find(i) == group.compIDs.end()) {
							  // mark for deletion
							toRemove.insert(match.begin()->at(i));
						} else {
							  // replace compID node with group description
							nodeLabel[ boost::vertex( match.begin()->at(i), mol) ] = group.description;
						}
					}
					  // remove all remaining nodes in decreasing index order
					for (std::set<size_t>::const_reverse_iterator n = toRemove.rbegin(); n != toRemove.rend(); ++n) {
						atom = boost::vertex( *n, mol );
						  // clear adjacent edges
						boost::clear_vertex( atom, mol );
						  // remove vertex
						boost::remove_vertex( atom, mol );
					}
				}

			} while( matchFound );

			// garbage collection
			for( size_t i=0; i<constraints.size(); ++i) {
				delete constraints[i];
			}
			constraints.clear();
		}
	}

//##############################################################################

	bool
	MoleculeUtil::
	isGroupLabel( const std::string & nodeLabel )
	{
		size_t closePos = std::string::npos;
				// check for leading bracket
		return (*(nodeLabel.begin()) == '{')
				// check for closing bracket
				&& (closePos = nodeLabel.find('}',1) != std::string::npos)
				// check if at least of length 3
				&& (nodeLabel.size() > 2)
				;

	}

//##############################################################################

	std::string
	MoleculeUtil::
	getGroupLabel( const std::string & nodeLabel )
	{
		std::string groupLabel = "";

		  // check if group label and return accordingly
		if (isGroupLabel(nodeLabel)) {
			  // grep group label part
			groupLabel = nodeLabel.substr(0,nodeLabel.find('}',2)+1);
		}

		return groupLabel;
	}

//##############################################################################

 	
 } // namespace chem
} // namespace ggl

