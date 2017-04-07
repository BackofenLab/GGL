// -*- C++ -*-
// SVM with stochastic gradient
// Copyright (C) 2007- Leon Bottou

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA



// $Id: svmsgd.h,v 1.2 2011/05/25 09:40:19 mmann Exp $

#ifndef SGM_SVMSGD_H
#define SGM_SVMSGD_H

//#include <iostream>
//#include <fstream>
//#include <iomanip>
//#include <string>
//#include <map>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <limits>

#include "sgd/vectors.h"
#include "sgd/timer.h"
#include "sgd/svmmodel.h"

namespace sgd {

using namespace std;


typedef vector<SVector> xvec_t;
typedef vector<double> yvec_t;


// Available losses
#define HINGELOSS 1
#define SMOOTHHINGELOSS 2
#define SQUAREDHINGELOSS 3
#define LOGLOSS 10
#define LOGLOSSMARGIN 11

// Select loss
#define LOSS HINGELOSS

// Zero when no bias
// One when bias term
#define BIAS 1


inline
double loss(double z)
{
#if LOSS == LOGLOSS
	if (z > 18)
		return exp(-z);
	if (z < -18)
		return -z;
	return log(1+exp(-z));
#elif LOSS == LOGLOSSMARGIN
	if (z > 18)
		return exp(1-z);
	if (z < -18)
		return 1-z;
	return log(1+exp(1-z));
#elif LOSS == SMOOTHHINGELOSS
	if (z < 0)
		return 0.5 - z;
	if (z < 1)
		return 0.5 * (1-z) * (1-z);
	return 0;
#elif LOSS == SQUAREDHINGELOSS
	if (z < 1)
		return 0.5 * (1 - z) * (1 - z);
	return 0;
#elif LOSS == HINGELOSS
	if (z < 1)
		return 1 - z;
	return 0;
#else
# error "Undefined loss"
#endif
}

inline
double dloss(double z)
{
#if LOSS == LOGLOSS
	if (z > 18)
		return exp(-z);
	if (z < -18)
		return 1;
	return 1 / (exp(z) + 1);
#elif LOSS == LOGLOSSMARGIN
	if (z > 18)
		return exp(1-z);
	if (z < -18)
		return 1;
	return 1 / (exp(z-1) + 1);
#elif LOSS == SMOOTHHINGELOSS
	if (z < 0)
		return 1;
	if (z < 1)
		return 1-z;
	return 0;
#elif LOSS == SQUAREDHINGELOSS
	if (z < 1)
		return (1 - z);
	return 0;
#else
	if (z < 1)
		return 1;
	return 0;
#endif
}


// -- stochastic gradient

	class SvmSgd
	{
	private: // members

		double  t;
		double  lambda;

		 //! the model to be trained and used for tests
		SvmModel model;

	public: // functions

		SvmSgd(double lambda);

		 /*!
		  * Initializes with a given model, e.g. to apply only tests.
		  *
		  * @param model the model to load
		  * @param lambda
		  */
		SvmSgd( const SvmModel& model, double lambda );


		 /*!
		  * Trains the SVM model based on the given data.
		  *
		  * @param trainFeatures the feature vectors to be used for training
		  * @param trainTarget the target values to be optimized, i.e. the true
		  *            classification
		  * @param imin the first index within trainFeatures to be used for
		  *            training
		  * @param imax the last index within trainFeatures to be used for
		  *            training
		  * @param epochs the number of training iterations to perform
		  *
		  */
		void train(	const xvec_t & trainFeatures
					, const yvec_t & trainTarget
					, size_t imin = 0
					, size_t imax = numeric_limits<size_t>::max()
//					, const char *prefix = ""
					, const size_t epochs = 1 );

		 /*!
		  * Tests the SVM model based on the given data.
		  *
		  * @param testFeatures the feature vectors to be used for testing
		  * @param testTarget the targeted values to check if the prediction
		  *            is correct
		  * @param imin the first index within testFeatures to be used for
		  *            testing
		  * @param imax the last index within testFeatures to be used for
		  *            testing
		  *
		  * @return first = missclassification in %; second = cost
		  */
		pair< double, double >
		test(	const xvec_t & testFeatures
					, const yvec_t & testTarget
					, size_t imin = 0
					, size_t imax = numeric_limits<size_t>::max()
//					, const char *prefix = ""
				);

		 /*!
		  * Access to the trained model.
		  * @return the trained model.
		  */
		const SvmModel & getModel() const;

	};






//// --- options
//
//string trainfile;
//string testfile;
//double lambda = 1e-4;
//int epochs = 5;
//int maxtrain = -1;
//
//void
//usage();
//
//void
//parse(int argc, const char **argv);
//
//
//// --- loading data
//
//int dim;
//xvec_t xtrain;
//yvec_t ytrain;
//xvec_t xtest;
//yvec_t ytest;
//
//void
//load(const char *fname, xvec_t &xp, yvec_t &yp);

} // namespace

#endif
