#!/bin/sh
# 
# Run this before configure
#
# $Id: autogen.sh,v 1.1 2008/11/14 12:50:48 xtof Exp $

set -e
#test -d config || mkdir config
# Produce aclocal.m4, so autoconf gets the automake macros it needs
echo "Creating aclocal.m4..."
aclocal
#autoheader

# Produce all the `Makefile.in's, verbosely, and create neat missing things
# like `libtool', `install-sh', etc.
automake --add-missing --gnu --copy

# If there's a config.cache file, we may need to delete it.  
# If we have an existing configure script, save a copy for comparison.
if [ -f config.cache ] && [ -f configure ]; then
  cp configure configure.$$.tmp
fi

# Produce ./configure
echo "Creating configure..."
autoconf

echo ""
echo "You can run ./configure [--prefix=$HOME] now."
echo ""
