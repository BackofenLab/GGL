
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



// $Id: run_svmsgd.cc,v 1.2 2011/11/19 23:05:55 mmann Exp $


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <cmath>

#include "sgd/svmsgd.h"
#include "sgd/vectors.h"
#include "sgd/timer.h"

using namespace std;
using namespace sgd;


//typedef vector<SVector> xvec_t;
//typedef vector<double> yvec_t;


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




// --- options

string trainfile;
string testfile;
double lambda = 1e-4;
int epochs = 5;
int maxtrain = -1;

void
usage()
{
  cerr << "Usage: svmsgd [options] trainfile [testfile]" << endl
       << "Options:" << endl
       << " -lambda <lambda>" << endl
       << " -epochs <epochs>" << endl
       << " -maxtrain <n>" << endl
       << endl;
  exit(10);
}

void
parse(int argc, const char **argv)
{
  for (int i=1; i<argc; ++i )
    {
      const char *arg = argv[i];
      if (arg[0] != '-')
        {
          if (trainfile.empty())
            trainfile = arg;
          else if (testfile.empty())
            testfile = arg;
          else
            usage();
        }
      else
        {
          while (arg[0] == '-') arg += 1;
          string opt = arg;
          if (opt == "lambda" && i+1<argc)
            {
              lambda = atof(argv[++i]);
              cout << "Using lambda=" << lambda << "." << endl;
              assert(lambda>0 && lambda<1e4);
            }
          else if (opt == "epochs" && i+1<argc)
            {
              epochs = atoi(argv[++i]);
              cout << "Going for " << epochs << " epochs." << endl;
              assert(epochs>0 && epochs<1e6);
            }
          else if (opt == "maxtrain" && i+1<argc)
            {
              maxtrain = atoi(argv[++i]);
              assert(maxtrain > 0);
            }
          else
            usage();
        }
    }
  if (trainfile.empty())
    usage();
}


// --- loading data

int dim;
xvec_t xtrain;
yvec_t ytrain;
xvec_t xtest;
yvec_t ytest;

void
load(const char *fname, xvec_t &xp, yvec_t &yp)
{
  cout << "Loading " << fname << "." << endl;

  ifstream f;
  f.open(fname);
  if (! f.good())
    {
      cerr << "ERROR: cannot open " << fname << "." << endl;
      exit(10);
    }
  int pcount = 0;
  int ncount = 0;

  bool binary;
  string suffix = fname;
  if (suffix.size() >= 7)
    suffix = suffix.substr(suffix.size() - 4);
  if (suffix == ".dat")
    binary = false;
  else if (suffix == ".bin")
    binary = true;
  else
    {
      cerr << "ERROR: filename should end with .bin or .dat" << endl;
      exit(10);
    }

  while (f.good())
    {
      SVector x;
      double y;
      if (binary)
        {
          y = (f.get()) ? +1 : -1;
          x.load(f);
        }
      else
        {
          f >> y >> x;
        }
      if (f.good())
        {
          assert(y == +1 || y == -1);
          xp.push_back(x);
          yp.push_back(y);
          if (y > 0)
            pcount += 1;
          else
            ncount += 1;
          if (x.size() > dim)
            dim = x.size();
        }
    }
  cout << "Read " << pcount << "+" << ncount
       << "=" << pcount + ncount << " examples." << endl;
}



int
main(int argc, const char **argv)
{
  parse(argc, argv);

  // load training set
  load(trainfile.c_str(), xtrain, ytrain);
  cout << "Number of features " << dim << "." << endl;
  int imin = 0;
  int imax = xtrain.size() - 1;
  if (maxtrain > 0 && imax >= maxtrain)
    imax = imin + maxtrain -1;
  // prepare svm
  SvmSgd svm( lambda );
  Timer timer;

  // load testing set
  if (! testfile.empty())
    load(testfile.c_str(), xtest, ytest);
  int tmin = 0;
  int tmax = xtest.size() - 1;

      timer.start();
      cout << "train: " << "Training on [" << imin << ", " << imax << "] epochs = "<<epochs<<"." << endl;
      svm.train( xtrain, ytrain, imin, imax, epochs);
      double wnorm =  dot(svm.getModel().w,svm.getModel().w) * svm.getModel().wscale * svm.getModel().wscale;
      cout << "train: " << setprecision(6)
           << "Norm: " << wnorm << ", Bias: " << svm.getModel().bias << endl;

      timer.stop();
      cout << "Total training time " << setprecision(6)
           << timer.elapsed() << " secs." << endl;

	  cout << "train: " << "Testing on [" << imin << ", " << imax << "]." << endl;
      std::pair<double,double> res = svm.test( xtrain, ytrain, imin, imax );
        cout << "train: " << setprecision(4)
             << "Misclassification: " << res.first << "%." << endl;
        cout << "train: " << setprecision(12)
             << "Cost: " << res.second << "." << endl;
      if (tmax >= tmin) {
    	  cout << "test:  " << "Testing on [" << tmin << ", " << tmax << "]." << endl;
        std::pair<double,double> res = svm.test( xtest, ytest, tmin, tmax);
        cout << "test:  " << setprecision(4)
             << "Misclassification: " << res.first << "%." << endl;
        cout << "test:  " << setprecision(12)
             << "Cost: " << res.second << "." << endl;
      }

      cout <<"\nmodel:\n"
				<<" wscale = "<<svm.getModel().wscale <<"\n"
				<<" bias   = "<<svm.getModel().bias <<"\n"
				<<" w      = "<<svm.getModel().w
				<<endl;
}
