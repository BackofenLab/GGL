/*
 * svmmodel.h
 *
 *  Created on: 12.05.2011
 *      Author: mmann
 */

#ifndef SGD_SVMMODEL_H_
#define SGD_SVMMODEL_H_

#include "sgd/vectors.h"

namespace sgd {

	 /*!
	  * A model to be trained by a linear SVM based on a sparse vector
	  * representation of the learned weights.
	  */
	class SvmModel {

	public: // members

		  //! the weights for each feature
		SVector w;

		  //! the weight scale for better precision
		double  wscale;

		  //! the bias that marks the treshold to reach to be classified as a
		  //! positive instance
		double bias;


	public: // functions

		 /*!
		  * Initializes the model to hold "dim" features and using the given
		  * initial weight scale.
		  *
		  * @param wscaleStart the initial value for the weight scaling
		  * @param bias
		  */
		SvmModel( const double wscaleStart = 1.0
					, const double biasStart = 0.0 )
		 : w(), wscale(wscaleStart), bias(biasStart)
		{}

		  /*!
		   * copy construction
		   * @param model the model to copy
		   */
		SvmModel( const SvmModel& model )
		 : w(model.w), wscale(model.wscale), bias(model.bias)
		{}

		 /*!
		  * Computes the prediction value
		  *
		  * ((dot( w, toTest ) * wscale) + bias)
		  *
		  * @param toTest the instance to test
		  * @return the predicted value
		  */
		double
		predictValue( const SVector& toTest ) const {

			  // get the accumulated weight for the test features
			return ((dot( w, toTest ) * wscale) + bias);

		}

		 /*!
		  * Predicts if an instance is positive, i.e.
		  *
		  * "predictValue(toTest) >= 0"
		  *
		  * @param toTest the instance to test
		  * @return true if the instance is positive; false otherwise
		  */
		bool
		predict( const SVector& toTest ) const {

			  // get the accumulated weight for the test features
			return predict(predictValue(toTest));

		}

		 /*!
		  * Predicts if an instance is positive, i.e.
		  *
		  * "predValue >= 0"
		  *
		  * @param predValue the predicted value to decide for
		  * @return true if the instance is positive; false otherwise
		  */
		bool
		predict( const double predValue ) const {

			  // decides based on the prediction value
			return predValue >= 0.0;

		}

	};

} // namespace sgd

#endif /* SVMMODEL_H_ */
