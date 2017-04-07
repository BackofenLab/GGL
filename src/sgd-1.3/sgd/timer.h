// -*- C++ -*-
// A simple timer.
// Copyright (C) 2007- Leon Bottou

// This library is free software; you can redistribute it and/or
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



// $Id: timer.h,v 1.2 2011/05/12 13:50:38 mmann Exp $


#ifndef SGD_TIMER_H
#define SGD_TIMER_H 1

namespace sgd {

	class Timer
	{
	public:
	  Timer();
	  void reset();
	  double start();
	  double stop();
	  double elapsed();
	private:
	  double a, s;
	  int    r;
	};

} // namespace

/* -------------------------------------------------------------
   Local Variables:
   c++-font-lock-extra-types: ("\\sw+_t" "[A-Z]\\sw*[a-z]\\sw*" "std::\\sw+_t")
   End:
   ------------------------------------------------------------- */


#endif
