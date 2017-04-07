/*
 * MoleculeComponent.cc
 *
 *  Created on: 09.11.2010
 *      Author: mmann
 */

#include "ggl/chem/MoleculeComponent.hh"

namespace ggl {
namespace chem {


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

MoleculeComponent
::RingFragment
::RingFragment( const RingFragmentType& type_
				, const RingFragmentList& fragment_)
: type(type_), fragment(fragment_)
{}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


MoleculeComponent
:: MoleculeComponent()
 :	description("")
	, priority(0)
	, freeEnergy(0.0)
	, pattern()
	, compIDs()
	, constraints()
	, ringFragments()
{
}


//////////////////////////////////////////////////////////////////////////////


MoleculeComponent
:: MoleculeComponent( const MoleculeComponent& toCopy )
 :	description( toCopy.description )
	, priority( toCopy.priority )
	, freeEnergy( toCopy.freeEnergy )
	, pattern( toCopy.pattern )
	, compIDs( toCopy.compIDs )
	, constraints( sgm::Pattern::copyConstraintVec(toCopy.constraints) )
	, ringFragments( toCopy.ringFragments )
{ }

//////////////////////////////////////////////////////////////////////////////


MoleculeComponent &
MoleculeComponent
:: operator = ( const MoleculeComponent& toCopy )
{
	  // flat copy data update
	this->description = toCopy.description ;
	this->priority = toCopy.priority ;
	this->freeEnergy = toCopy.freeEnergy ;
	this->pattern = toCopy.pattern ;
	this->compIDs = toCopy.compIDs ;
	this->ringFragments = toCopy.ringFragments ;
	  // manual cloning of constraints
	for (size_t i=0; i<constraints.size(); ++i ) {
		delete constraints[i];
	}
	constraints = sgm::Pattern::copyConstraintVec(toCopy.constraints);
	  // access to changed object
	return *this;
}


//////////////////////////////////////////////////////////////////////////////


MoleculeComponent
:: ~MoleculeComponent()
{
	  // clear heap memory
	for (size_t i=0; i<constraints.size(); ++i ) {
		delete (constraints[i]);
	}
	constraints.clear();
}


//////////////////////////////////////////////////////////////////////////////

}  // namespace chem
} // namespace ggl
