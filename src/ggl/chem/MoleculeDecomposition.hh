
#ifndef MOLECULEDECOMPOSITION_HH_
#define MOLECULEDECOMPOSITION_HH_

#include <stdexcept>

#include "sgm/Graph_Interface.hh"
#include "sgm/SubGraphMatching.hh"
#include "sgm/Match_Reporter.hh"
#include "sgm/RingReporter.hh"
#include "sgm/GM_vf2.hh"
#include "sgm/SGM_vf2.hh"

#include "ggl/Graph.hh"
#include "ggl/chem/Molecule.hh"
#include "ggl/chem/MoleculeComponent.hh"

namespace ggl {
 namespace chem {


	 /*! @brief Molecule energy estimation ala Jankowski et al.
	  *
	  * Molecule decomposition approach to analyze a molecules thermodynamics.
	  * The approach follows the paper by Jankowski et al. (2008)
	  *
	  * Group Contribution Method for Thermodynamic Analysis of Complex
	  *  Metabolic Networks
	  * Matthew D. Jankowski, Christopher S. Henry, Linda J. Broadbelt
	  *  and Vassily Hatzimanikatis
	  * Biophysical Journal, Volume 95, Issue 3, 1487-1499, 1 August 2008
	  * doi:10.1529/biophysj.107.124784
	  *
	  *
	  * @author Martin Mann - 2010 - http://www.bioinf.uni-freiburg.de/~mmann/
	  *
	  */
	class MoleculeDecomposition : public sgm::Match_Reporter, public sgm::RingReporter {
	public:

		  /*!
		   * Abstract class to report matched molecule components and
		   * accordingly the labeled match-graph.
		   */
		class DecompositionReporter {
		public:
			  //! construction
			DecompositionReporter(){};
			  //! destruction
			virtual ~DecompositionReporter(){};

			  /*!
			   * Reports a component that was matched and its according ID
			   * in the final match graph
			   * @param component the component that was matched
			   * @param matchCount the number of matches found
			   * @param matchID the ID of the component in the final match graph
			   */
			virtual
			void
			reportComponent(	const MoleculeComponent & component,
								const size_t matchCount,
								const size_t matchID ) = 0;
			  /*!
			   * Reports an interaction pattern that was matched
			   * @param interactionDescription the interaction that was matched
			   * @param matchCount the number of matches found
			   */
			virtual
			void
			reportInteraction(	const std::string & interactionDescription,
								const size_t matchCount ) = 0;

			  /*!
			   * Reports the matching of a full small molecule.
			   * @param smallMolecule the small molecule matched
			   */
			virtual
			void
			reportSmallMolecule( const MoleculeComponent & smallMolecule ) = 0;

			  /*!
			   * Reports the final match graph where each matched component ID
			   * is given as class information of the matched nodes.
			   *
			   * @param graph the matched graph
			   */
			virtual
			void
			reportMatchGraph(	const sgm::Graph_Interface & graph ) = 0;

			  /*!
			   * Reports whether or not the complete graph (all nodes) were
			   * successfully mapped or not.
			   *
			   * @param matchComplete true, if all nodes have been mapped;
			   *                      false otherwise.
			   */
			virtual
			void
			reportMatchComplete( const bool matchComplete ) = 0;

			  /*!
			   * Reports the matching of a polyphosphate chain. Only the number
			   * and match ID of the middle phosphates is given. The end
			   * phosphate is reported as independent component.
			   * @param numOfMiddlePhosphates the number of middle phosphate
			   *          matches found
			   * @param matchID the ID of the component in the final match graph
			   */
			virtual
			void
			reportPolyPhosphate( const size_t numOfMiddlePhosphates
								, const size_t matchID ) = 0;
		};


	protected:

		  /*!
		   * helper class to order MoleculeComponent pointers by the priority
		   * of the MoleculeComponents within an ordered STL container.
		   */
		class priority_order {
		public:

			  /*!
			   * operator to order MoleculeComponent pointers by the priority
			   * of the MoleculeComponents.
			   * @param c1 the first component
			   * @param c2 the second component
			   * @return true if the first component has higher priority compared
			   *         to the second; false otherwise
			   */
			bool operator()(const MoleculeComponent & c1,
							const MoleculeComponent & c2 ) const;
		};

		typedef std::set< MoleculeComponent, priority_order > ComponentContainer;

		  //! set of all groups that are matched exclusively starting with the
		  //! one with lowest priority value
		static ComponentContainer groups;

		  //! set of all interaction patterns to be screened for
		static ComponentContainer interactions;

		  //! set of all small molecules to be screened for
		static ComponentContainer smallMolecules;

		  //! correction term if a molecule is a hydrocarbon, i.e.
		  //! it contains only H and C
		static double correctionHydroCarbon;

		  //! correction term for each heteroaromatic ring within the molecule
		static double correctionPerHeteroaromaticRing;

		  //! correction term for each three-membered ring
		  //! independently of the atom labels
		static double correctionPerThreeMemberedRing;

		  //! correction term for each amide within the molecule; note,
		  //! for nitrogen atom involved in two amides only one occurrence is
		  //! considered
		static double correctionPerAmide;

		  //! correction term for each thioester within the molecule; note,
		  //! for sulfur atom involved in two amides only one occurrence is
		  //! considered
		static double correctionPerThioester;

		  //! correction term for vicinal chlorines within the molecule
		static double correctionVicinalChlorine;

		  //! energy contribution term for a middle phosphate within a phosphate
		  //! chain
		static double contributionMiddlePhosphate;

		  //! Access to the groups member. This function
		  //! initializes the container if not already done.
		  //! @return the groups container to be matched.
		  //! @thrown std::logic_error in case the group initialization failed
		  //!         which should never happen
		static
		const ComponentContainer &
		getGroups( void ) throw (std::logic_error);

		  //! Access to the interactions member. This function
		  //! initializes the container if not already done.
		  //! @return the interactions container to be matched.
		  //! @thrown std::logic_error in case the group initialization failed
		  //!         which should never happen
		static
		const ComponentContainer &
		getInteractions( void ) throw (std::logic_error);

		  //! Access to the smallMolecules member. This function
		  //! initializes the container if not already done.
		  //! @return the smallMolecules container to be matched.
		  //! @thrown std::logic_error in case the group initialization failed
		  //!         which should never happen
		static
		const ComponentContainer &
		getSmallMolecules( void ) throw (std::logic_error);

		  //! Computes the number of amide matches within the molecule where for
		  //! nitrogen atom involved in two amides only one occurrence is
		  //! counted
		  //! @param matcher the sub graph matcher to be used
		  //! @param molGraph the molecule graph to check
		  //! @return the number of amides within the molecule
		static
		size_t
		getAmideNumber( sgm::SubGraphMatching & matcher
						, sgm::Graph_Interface & molGraph );

		  //! Computes the number of thioester matches within the molecule
		  //! where for sulfur atoms involved in two thioester only one
		  //! occurrence is counted
		  //! @param matcher the sub graph matcher to be used
		  //! @param molGraph the molecule graph to check
		  //! @return the number of thioester within the molecule
		static
		size_t
		getThioesterNumber( sgm::SubGraphMatching & matcher
							, sgm::Graph_Interface & molGraph );

		  //! Corrects the mapping of cyclic phosphates. Due to their symmetric
		  //! symmetric pattern they are matched twice. Thus, one oxygen is
		  //! mapped too much.
		  //! @param curMatchedIDs the node sets for each match that has to be
		  //!         corrected
		static
		void
		correctCyclicPhosphates( std::set< std::set<size_t> > & curMatchedIDs );

		  //! Corrects the mapping of enclosed primary phosphates. Due to their
		  //! symmetric symmetric pattern they are matched all along the
		  //! phosphate fragment and not only at the ends.
		  //! @param curMatchedIDs the node sets for each match that has to be
		  //!         corrected
		static
		void
		correctEnclosedPhosphates( std::set< std::set<size_t> > & curMatchedIDs );

	protected:

		  //! the energy of the currently parsed molecule
		double energy;
		  //! the currently matched MoleculeComponent
		const MoleculeComponent* curComponent;
		  //! the matched set of node indices for current MoleculeComponent
		std::set< std::set<size_t> > curMatchedIDs;

		  //! set of all nodes participating in a ring within the current
		  //! molecule handled
		RingNodes curMolRingNodes;

		  //! descriptor of a ring within the molecule
		class RingDescriptor {
		public:
			typedef std::pair<size_t,size_t> Edge;
			typedef std::set<Edge> EdgeSet;
		public:
			  //! edge set of this ring
			EdgeSet edges;

			  //! constructor
			  //! @param edges the edge set to set
			RingDescriptor( const EdgeSet& edges )
			 : edges(edges)
			{}

			  //! constructor
			  //! @param ringList the ring to describe
			RingDescriptor( const sgm::RingReporter::RingList & ringList );

			  //! decides order based on ring set size
			bool
			operator <( const RingDescriptor& toCompare ) const;

			  /*!
			   * Checks whether or not this ring is contained within the other
			   * larger ring as a fuse with some other ring. Therefore, the
			   * edge set difference leaves only one edge.
			   * @param largerRing the larger ring this ring might be contained
			   * @return true, if the set difference of this edge set with the
			   *         edge set of the larger ring leaves only one edge
			   */
			bool
			isContained( const RingDescriptor& largerRing ) const;

		};
		  //! string encoding of all rings within the current molecule grouped
		  //! by their type
		std::map< MoleculeComponent::RingFragmentType, std::vector< RingDescriptor > > curMolRings;

		  //! an optional decomposition reporter
		DecompositionReporter* reporter;

		  //! the graph matching implementation to be used
		sgm::GraphMatching * fullMatcher;
		  //! the sub graph matching implementation to be used
		sgm::SubGraphMatching * matcher;

		  /*!
		   * Node constraint to ensure that a matched node is part of a ring
		   * based on the given ring node set information.
		   *
		   * This increases the probability that the ring fragment constraints,
		   * that are checked only for full matches, are fulfilled.
		   */
		class MC_MC_RingNode : public sgm::MC_Node {

		public:

			  //! the set of possible ring node ids the constrainedNodeID can
			  //! be matched on
			const RingNodes * ringNodes;

			  //! construction
			  //! @param constrainedNodeID the node ID to be constrained
			  //! @param ringNodes the set of possible ring node ids the
			  //! constrainedNodeID can be matched on
			MC_MC_RingNode( const size_t constrainedNodeID
							, const RingNodes & ringNodes)
			  :	sgm::MC_Node(constrainedNodeID), ringNodes(&ringNodes)
			{}

			  //! copy construction
			  //! @param toCopy the object to copy
			MC_MC_RingNode( const MC_MC_RingNode& toCopy )
			  :	sgm::MC_Node(toCopy.constrainedNodeID)
				, ringNodes(toCopy.ringNodes)
			{}

			  //! destruction
			virtual
			~MC_MC_RingNode()
			{}

			 /*!
			  * Creates a new heap object that equals the current object.
			  * NOTE: YOU have to delete it later on! There is no garbage
			  *       collection!
			  * @return a new allocated object that equals this
			  */
			virtual
			MC_Node*
			clone( void ) const {
				return new MC_MC_RingNode(*this);
			}

			 /*!
			  * Checks whether or not a given label is part of the constraint
			  * information. This check is needed by some parsers to verify the
			  * wildcard definition.
			  *
			  * @param label the label of interest
			  * @return false, since no label information is used/stored
			  */
			virtual
			bool
			isConstrainedLabel(const std::string& label) const {
				return false;
			}

			 /*!
			  * Creates a new heap object that equals the current
			  * object but uses the new indices given by old2newIndexMapping.
			  * NOTE: YOU have to delete it later on! There is no garbage
			  *       collection!
			  * @param old2newIndexMapping the index mapping to be used for the
			  *        remapping
			  * @param unmatchedIndex an optional specific index that marks
			  *        unmatched nodes within old2newIndexMapping. if this
			  *        constrains one of these nodes, no remapping is done and
			  *        NULL is returned
			  * @return a new allocated object
			  */
			virtual
			MC_Node*
			remap( const sgm::Match & old2newIndexMapping, const size_t unmatchedIndex = UINT_MAX )
			{
				assert(this->constrainedNodeID < old2newIndexMapping.size());
				  // check if this node is an unmatched node and thus to be ignored
				if (old2newIndexMapping.at(this->constrainedNodeID)==unmatchedIndex) {
					return NULL;
				}
				  // create copy
				MC_MC_RingNode* copy = new MC_MC_RingNode(*this);
				  // do remapping
				copy->constrainedNodeID = old2newIndexMapping.at(this->constrainedNodeID);
				  // return remapped copy
				return copy;
			}



			 /*!
			  * Checks whether or not a match on a given target fulfills the
			  * additional node constraint for the pattern matching.
			  *
			  * @param pattern the pattern graph that was matched
			  * @param target the target graph the pattern was matched on
			  * @param matchedTargetID the matched node index within
			  *        the target graph
			  * @return true if the match is valid; false if the constraint is
			  *         violated
			  */
			virtual
			bool
			isValidMatch(	const sgm::Pattern_Interface & pattern,
							const sgm::Graph_Interface & target,
							const size_t matchedTargetID ) const
			{
				assert(ringNodes != NULL);
				  // check whether or not the given target id is among the
				  // known ring nodes
				return ringNodes->find(matchedTargetID) != ringNodes->end();
			}



		};

	protected:

		  /*!
		   * Checks whether or not a given ring fragment is present in the
		   * molecule.
		   *
		   * @param ringFragment the ring fragment to check
		   * @param match the node matching that triggered the check
		   * @return true, if the ring fragment is matched; false otherwise
		   */
		bool
		isValidRingFragment( const MoleculeComponent::RingFragment & ringFragment
							, const sgm::Match & match);

		  //! node label map typedef for shorter parameter lists
		typedef boost::property_map< Molecule, PropNodeLabel>::type NodeLabelMap;
		  //! node index map typedef for shorter parameter lists
		typedef boost::property_map< Molecule, PropNodeIndex>::const_type NodeIndexMap;

		  /*!
		   * Extend the labeling of a phosphate chain starting from a given
		   * phosphorus node using a specified component label.
		   *
		   * @param m the molecule to label
		   * @param mNodeLabel the node label access for m
		   * @param mNodeIndex the node index access for m
		   * @param idxNextP the phosphorus node to start the labeling
		   * @param compID the label for the phosphate chain
		   * @return the length of the chain
		   */
		size_t
		extendPolyphosphateChain( Molecule & m
							, NodeLabelMap & mNodeLabel
							, NodeIndexMap & mNodeIndex
							, const size_t idxNextP
							, const size_t compID
							);

	public:

		  //! construction
		  //! @param fullMatcher the graph matching implementation to be used
		  //! @param matcher the sub graph matching implementation to be used
		MoleculeDecomposition( sgm::GraphMatching & fullMatcher
								, sgm::SubGraphMatching & matcher );

		  //! construction with a given reporter
		  //! @param fullMatcher the graph matching implementation to be used
		  //! @param matcher the sub graph matching implementation to be used
		  //! @param reporter the decomposition reporter to be used
		MoleculeDecomposition(  sgm::GraphMatching & fullMatcher
								, sgm::SubGraphMatching & matcher
								, DecompositionReporter & reporter );

		  //! destruction
		virtual ~MoleculeDecomposition();

		  /*!
		   * Decomposes the given molecule into MoleculeComponents and derives
		   * an estimate of its energy in kcal/mol.
		   *
		   * @param mol the molecule to be analyzed
		   * @return an estimate of its energy in kcal/mol
		   */
		double
		getEnergy( const Molecule & mol );

	public:

		  //! Handles the occurence of a MoleculeComponent within a decomposed
		  //! molecule graph.
		  //!
		  //! @param componentPattern the pattern of the MoleculeComponent to
		  //! be matched.
		  //! @param targetMol the molecule graph the pattern was found within
		  //! @param match contains the indices of the matched pattern nodes in
		  //! the target graph. match[i] corresponds to the mapping of the ith
		  //! vertex in the pattern graph.
		void
		reportHit (	const sgm::Pattern_Interface & componentPattern,
					const sgm::Graph_Interface & targetMol,
					const sgm::Match & match );

		  /*!
		   * Is called to report a ring. The ring is stored for lookup of ring
		   * constraints of the MoleculeComponents matched.
		   * @param graph the graph that contains the ring
		   * @param ringList the ring to report
		   */
		void
		reportRing( const sgm::Graph_Interface& graph, const RingList & ringList );

	};



 } // namespace chem
} // namespace ggl

#include "ggl/chem/MoleculeDecomposition.icc"

#endif /* MOLECULEDECOMPOSITION_HH_ */
