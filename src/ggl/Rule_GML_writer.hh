#ifndef GGL_RULE_GML_WRITER_HH_
#define GGL_RULE_GML_WRITER_HH_

#include <iostream>

#include "ggl/Rule.hh"


namespace ggl
{


  /*! @brief Rule GML writer
   *
   * Algorithm wrapper class to write a ggl rule in GGL GML format to stream.
   * 
   * @author Martin Mann (c) 2013 http://www.bioinf.uni-freiburg.de/~mmann/
   * 
   */
class Rule_GML_writer
{
public:
	  //! construction
	Rule_GML_writer();
	  //! destruction
	~Rule_GML_writer();
	
	
	  /*!
	   * Writes a boost rule in GGL GML format to stream.
	   * 
	   * NOTE: currently no constraints are printed...
	   *
	   * @param out the stream to write to
	   * @param rule the rule to print
	   * @param withSpaces if true newlines and white spaces are used to write a
	   *        user friendly and readable output, otherwise a one-line output
	   *        with minimal space requirement is produced. 
	   */
	static 
	void
	write(	std::ostream& out
			, const Rule & rule
			, const bool withSpaces = true );
	
	  /*!
	   * Writes a boost rule in compacted GGL GML format to stream.
	   *
	   * Therin, each node/label occurs only once. To indicate label changes or
	   * edge insertions/deletions, a label is in one of the following 4 forms:
	   *
	   * - "C" : the according node/edge has this label left and right = context
	   * - "L|" : the according edge is removed and the left side label L is given
	   * - "|R" : the according edge is inserted and the right side label R is given
	   * - "L|R" : the node label changes from left label L to right label R or
	   *           the edge changes from valence L to valence R during the reaction
	   *
	   * @param out the stream to write to
	   * @param rule the rule to print
	   * @param withSpaces if true newlines and white spaces are used to write a
	   *        user friendly and readable output, otherwise a one-line output
	   *        with minimal space requirement is produced.
	   */
	static
	void
	writeCompact(	std::ostream& out
			, const Rule & rule
			, const bool withSpaces = true );

private:

	  /*!
	   * Writes the GML part for the given context
	   *
	   * @param graph the rule core graph
	   * @param context the context to print
	   * @param whether to print with whitespaces or minimal
	   *
	   * @return the GML notation for the current context
	   */
	static
	std::string
	getContextGML( const Rule::CoreGraph & graph
				, const Rule::RuleContext & context
				, const bool withSpaces );


};


} // namespace ggl


  // include implementation
#include "ggl/Rule_GML_writer.icc"

#endif /*RULE_GML_WRITER_HH_*/

