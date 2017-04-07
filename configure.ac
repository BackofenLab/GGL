#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
# 5 argument version only available with aclocal >= 2.64
AC_INIT([GGL], [4.1.2], [], [ggl], [http://www.bioinf.uni-freiburg.de/~mmann/])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_SRCDIR([src/bin/toyChem.cc])
AC_CANONICAL_HOST

m4_include([misc/m4/m4_ax_boost_base.m4])
m4_include([misc/m4/m4_ax_prog_doxygen.m4])

# backwards compatibility for autoconf >= 2.64
# PACKAGE_URL should be the fifth argument of AC_INIT
m4_ifndef([AC_PACKAGE_URL],
	[m4_define([AC_PACKAGE_URL], [http://www.bioinf.uni-freiburg.de/~mmann/])
	AC_DEFINE(PACKAGE_URL, "AC_PACKAGE_URL", [Package URL])
	AC_SUBST(PACKAGE_URL, AC_PACKAGE_URL)])

lt_enable_auto_import=""
case "$host_os" in 
	cygwin* | mingw* | cegcc*)
		AM_LDFLAGS="-Wl,--enable-auto-import $AM_LDFLAGS"
esac

# store current compiler flags to avoid default setup via AC_PROG_CXX and *_CC
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB # for VF2-library
# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS

# program checks
AC_CHECK_PROG([HAVE_PKG_CONFIG],[pkg-config],[yes],[no])
AC_CHECK_PROG([HAVE_PWD],[pwd],[yes],[no])

###############################################################################
# DOCUMENTATION BUILD CHECK
###############################################################################
dnl Find latex and friends

AC_MSG_CHECKING([whether to build manuals and tutorials])
AC_ARG_WITH(docs,
        AS_HELP_STRING([--without-docs],
                           [do not try to build manuals and tutorials]),
           ,with_docs=yes)
AC_MSG_RESULT($with_docs)

if test "$with_docs" = "yes"; then
        AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
        if test -z "$PDFLATEX"; then
        	AC_MSG_WARN([Cannot find 'pdflatex' in your path!])
        fi
        if test ! '(' -z "$PDFLATEX" ')'; then
        	BUILDDOC="yes"
        fi
else
        BUILDDOC="no"
fi

AM_CONDITIONAL(build_docs, [test "$BUILDDOC" = "yes"])


##########################################################################
# Checks for programs needed for tests
##########################################################################

# Check for 'diff' and get full path.
AC_ARG_VAR([DIFF],[the 'diff' program to use for test output comparison])
AC_PATH_PROG([DIFF],[diff],[])
if test "x$DIFF" = "x"; then
	AC_MSG_NOTICE([==> diff command not found!])
	AC_MSG_NOTICE([==> Set DIFF variable if present in non-standard path!])
	AC_MSG_ERROR([diff is mandatory to run the tests : will stop here!])
fi

# Check for 'cat' and get full path.
AC_ARG_VAR([CAT],[the 'cat' program used for printing test output files])
AC_PATH_PROG([CAT],[cat],[])
if test "x$CAT" = "x"; then
	AC_MSG_NOTICE([==> cat command not found!])
	AC_MSG_NOTICE([==> Set CAT variable if present in non-standard path!])
	AC_MSG_ERROR([cat is mandatory to run the tests : will stop here!])
fi

# Check for 'printf' and get full path.
AC_ARG_VAR([PRINTF],[the 'printf' program used to print test information])
AC_PATH_PROG([PRINTF],[printf],[])
if test "x$PRINTF" = "x"; then
	AC_MSG_NOTICE([==> printf command not found!])
	AC_MSG_NOTICE([==> Set PRINTF variable if present in non-standard path!])
	AC_MSG_ERROR([printf is mandatory to run the tests : will stop here!])
fi

##########################################################################


# setup of local VF2 library version
LOCAL_VF2_FOLDER="vflib2-080701"
AC_SUBST([LOCAL_VF2_FOLDER])

# setup of local SGD library version
LOCAL_SGD_FOLDER="sgd-1.3"
AC_SUBST([LOCAL_SGD_FOLDER])

# setup of local NSPDK library version
LOCAL_NSPDK_FOLDER="nspdk"
AC_SUBST([LOCAL_NSPDK_FOLDER])

# minimal required version of the OpenBabel library if depending module enabled
OPENBABEL_REQUIRED_VERSION=2.3.1
# minimal required version of the boost library (for bugfix of clear_vertex()) 
BOOST_REQUIRED_VERSION=1.50.0


###############################################################################
# PRELIMINARY CHECKS
###############################################################################


# Check for pkg-config manually first, as if its not installed the
# PKG_PROG_PKG_CONFIG macro won’t be defined.
if test x”$HAVE_PKG_CONFIG” == xno; then
	AC_MSG_ERROR([pkg-config is required for full dependency checks])
else
	# set default search pathes
	export PKG_CONFIG_PATH="usr/lib/pkgconfig/:usr/local/lib/pkgconfig/:$PKG_CONFIG_PATH"
	# call pgk-config setup 
	PKG_PROG_PKG_CONFIG 
fi

if test x”$HAVE_PWD” == xno; then
	AC_MSG_ERROR([pwd command is required for building the package])
fi


AM_INIT_AUTOMAKE([subdir-objects])

AC_LANG([C++])

###############################################################################
# BOOST CHECK
###############################################################################

AX_BOOST_BASE([$BOOST_REQUIRED_VERSION], [FOUND_BOOST=1;], [FOUND_BOOST=0;])

###############################################################################
# DOXYGEN SETUP
###############################################################################

DX_HTML_FEATURE([ON])
DX_CHM_FEATURE([OFF])
DX_CHI_FEATURE([OFF])
DX_MAN_FEATURE([OFF])
DX_RTF_FEATURE([OFF])
DX_XML_FEATURE([OFF])
DX_PDF_FEATURE([OFF])
DX_PS_FEATURE([OFF])
DX_INIT_DOXYGEN([$PACKAGE_NAME], [misc/doxygen/doxygen.config], [doc])

###############################################################################
# ANSI CONFORM CODE CHECK SETUP
###############################################################################

AC_MSG_CHECKING([whether to build strict ANSI C++ code])
ansi=no
AC_ARG_ENABLE([ANSI],
	[AS_HELP_STRING([--enable-ANSI],
	    [enable check for strict ANSI C++ source code (def=no)])],
  	[ansi="$enableval"])
AC_MSG_RESULT([$ansi])
if test x"$ansi" = x"yes"; then
	AM_CXXFLAGS="$AM_CXXFLAGS -ansi -pedantic"
fi

###############################################################################
# SETUP FOR COMPILATION WITH 2011 ISO C++ STANDARD 
###############################################################################

AC_MSG_CHECKING([whether to build with 2011 ISO C++ standard])
cpp11=no
AC_ARG_ENABLE([c++11],
	[AS_HELP_STRING([--enable-c++11],
	    [enable build with 2011 ISO C++ standard (Check if supported by your compiler! def=no)])],
  	[cpp11="$enableval"])
AC_MSG_RESULT([$cpp11])
if test x"$cpp11" = x"yes"; then
	AM_CXXFLAGS="$AM_CXXFLAGS -std=c++11"
fi


###############################################################################
# DEBUG SUPPORT SETUP
###############################################################################

AC_MSG_CHECKING([whether to build with debug information])
debuger=no
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
	    [enable debug data generation (def=no)])],
  	[debuger="$enableval"])
AC_MSG_RESULT([$debuger])
if test x"$debuger" = x"yes"; then
	AC_DEFINE([DEBUG])
	AM_CXXFLAGS="$AM_CXXFLAGS -g -O0 -Wall -Werror -Wno-uninitialized -Wno-deprecated"
else
	AC_DEFINE([NDEBUG])
	AM_CXXFLAGS="$AM_CXXFLAGS -O3 -fno-strict-aliasing -Wno-uninitialized -Wno-deprecated"
fi



###############################################################################
# GPROF PROFILER BUILT SUPPORT 
###############################################################################

AC_MSG_CHECKING([whether to build with '-pg' for gprof profiling])
profiling=no
AC_ARG_ENABLE([profiling],
	[AS_HELP_STRING([--enable-profiling],
	    [enable the built of gprof profiling information (-pg) (def=no)])],
  	[profiling="$enableval"])
AC_MSG_RESULT([$profiling])
if test x"$profiling" = x"yes"; then
	AM_CXXFLAGS="$AM_CXXFLAGS -pg"
fi



###############################################################################
# OPENBABEL-MODULE SETUP
###############################################################################

AC_ARG_ENABLE([OpenBabel],
    [AS_HELP_STRING([--enable-OpenBabel],
	    [enables the use of the OpenBabel-2 library. NOTE: THIS SUPPORT IS ONLY PRELIMINARY AND NOT FULLY TESTED! (def=no)])],
    [enable_OpenBabel=$enableval],
    [enable_OpenBabel=no])

AC_MSG_CHECKING([whether to use the OpenBabel-2 library])
AC_MSG_RESULT([$enable_OpenBabel])


###############################################################################
# STATIC LINKING SUPPORT
###############################################################################

AC_ARG_ENABLE([static-linking],
    [AS_HELP_STRING([--enable-static-linking],
	    [enables static linking of the executables (default=disabled)])],
    [enable_static_linking=$enableval],
    [enable_static_linking=no])

if test "$enable_static_linking" = "yes"; then
	AM_LDFLAGS="-static $AM_LDFLAGS"
else
	AS_CASE([$host_os],
			[cygwin*],
			[AC_MSG_NOTICE([NOTE: '-fPIC' g++ option not used for Windows Cygwin build])],
			[AM_CXXFLAGS="$AM_CXXFLAGS -fPIC"])
fi
AC_MSG_CHECKING([whether to link only with static libraries])
AC_MSG_RESULT([$enable_static_linking])


###############################################################################
# CHECKS FOR BEST AVAILABLE HASH OR MAP
###############################################################################


# check for unordered_map
AC_MSG_CHECKING([for std::unordered_map])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <unordered_map>]],
						[[std::unordered_map< int, int > myMap]])], 
	[AC_MSG_RESULT([yes]) 
		AC_SUBST([HAVE_UNORDERED_MAP],[1])],
	[AC_MSG_RESULT([no])
		AC_SUBST([HAVE_UNORDERED_MAP],[0])])
	
# check for tr1::unordered_map
AC_MSG_CHECKING([for std::tr1::unordered_map])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <tr1/unordered_map>]],
						[[std::tr1::unordered_map< int, int > myMap]])], 
	[AC_MSG_RESULT([yes]) 
		AC_SUBST([HAVE_TR1_UNORDERED_MAP],[1])],
	[AC_MSG_RESULT([no])
		AC_SUBST([HAVE_TR1_UNORDERED_MAP],[0])])
	
# check for gnu hash_map
AC_MSG_CHECKING([for GNU hash_map])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <ext/hash_map>]],
						[[__gnu_cxx::hash_map< int, int > myMap]])], 
	[AC_MSG_RESULT([yes]) 
		AC_SUBST([HAVE_GNU_HASH_MAP],[1])],
	[AC_MSG_RESULT([no])
		AC_SUBST([HAVE_GNU_HASH_MAP],[0])])



###############################################################################
# OPENBABEL LIBRARY PATH SUPPORT
###############################################################################

if test "$enable_OpenBabel" = "yes"; then
	# setup openbabel dependencies
	LIBS="-lm $LIBS"
	# add -lz if Z-library is present
	if test "$host_os" = "cygwin"; then
		AC_DEFINE([GGL_HOST_OS_CYGWIN])
		LIBS="-lz $LIBS"
	fi
#	AC_CHECK_LIB([z], [inflateEnd], [AC_DEFINE([HAVE_LIBZ],[1])],[])
	AC_CHECK_FUNCS([strcasecmp strncasecmp])
	
	# OpenBabel library path support, if not installed in usual directories
	AC_ARG_WITH([OpenBabel],
			[AS_HELP_STRING([--with-OpenBabel=prefix],
			    [try this for non-standard OpenBabel-library prefix install directory (optional)])], 
			[OPENBABELPATHSET=1],
			[OPENBABELPATHSET=0])
	if test  $OPENBABELPATHSET = 1 ; then
		export PKG_CONFIG_PATH="$with_OpenBabel/lib/pkgconfig/:$PKG_CONFIG_PATH"
	fi
fi








###############################################################################
# DEPENDENCY CHECKS
###############################################################################


###########################################
# check for stl headers
AC_CHECK_HEADERS([vector])

###########################################
# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T



###########################################
# check for OpenBabel library
if test "$enable_OpenBabel" = "yes"; then

	######### ensure no static linking is targeted ###########
	if test "$enable_static_linking" = "yes"; then
		AC_MSG_ERROR([OpenBabel is a dynamic library, thus no static linking possible : cannot compile and stop here !])
	fi

	AC_MSG_NOTICE([check if OpenBabel version >= $OPENBABEL_REQUIRED_VERSION])
	PKG_CHECK_MODULES([OpenBabel], 
				[openbabel-2.0 >= $OPENBABEL_REQUIRED_VERSION],
				[FOUND_OPENBABEL_2=1],
				[AC_MSG_NOTICE([$OpenBabel_PKG_ERRORS])
					FOUND_OPENBABEL_2=0])
	 # set compilation and linking flags
	AM_CXXFLAGS="$AM_CXXFLAGS $OpenBabel_CFLAGS"
	LIBS="$OpenBabel_LIBS $LIBS"
	 # set run path for linking the shared library
	if test  $OPENBABELPATHSET = 1 ; then
		AM_LDFLAGS="-Wl,-R,$with_OpenBabel/lib -Wl,-rpath,$with_OpenBabel/lib $AM_LDFLAGS"
	fi
else
   FOUND_OPENBABEL_2=0;
fi



###############################################################################
# COMPILER AND LINKER FLAG DISTRIBUTION
###############################################################################

AM_CXXFLAGS="$BOOST_CPPFLAGS $AM_CXXFLAGS -DBOOST_SYSTEM_NO_DEPRECATED"
AM_LDFLAGS="$BOOST_LDFLAGS $AM_LDFLAGS"

# distribute the configure based compiler and linker setups
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])
AC_SUBST([PREFIX_PATH], [${prefix}])


###############################################################################
# DEPENDENCY CHECK RESULT HANDLING
###############################################################################



if test "$enable_OpenBabel" = "yes"; then
	# error if OpenBabel path set but not found
	if test $FOUND_OPENBABEL_2 != 1; then
		AC_MSG_NOTICE([])
		AC_MSG_NOTICE([OpenBabel module:])
		AC_MSG_NOTICE([The following dependecies are not found:])
		if test $FOUND_OPENBABEL_2 != 1; then
			AC_MSG_NOTICE([ - OpenBabel-2 library version >= $OPENBABEL_REQUIRED_VERSION])
		fi
		if test  $OPENBABELPATHSET = 1 ; then
			AC_MSG_NOTICE([ The path '$with_OpenBabel' was set. Maybe wrong ???])
		else
			AC_MSG_NOTICE([ No non-default path was set. You might want to use '--with-OpenBabel=PREFIX' ?!?])
		fi
		AC_MSG_NOTICE([ ==> OpenBabel module is disabled])
		AC_MSG_NOTICE([])
	fi
fi

AM_CONDITIONAL([with_openbabel], [test "$FOUND_OPENBABEL_2" = "1"])

# FOUND_BOOST is only defined if want_boost is "yes"
if test $want_boost = "no" || test $FOUND_BOOST != 1; then
	AC_MSG_NOTICE([])
	AC_MSG_NOTICE([The Boost Graph Library (BGL) was not found!])
	AC_MSG_NOTICE([ You might want to use '--with-boost=PREFIX' ?!?])
	AC_MSG_NOTICE([])
	AC_MSG_ERROR([BGL is an essential dependency : have to stop here !])
fi

###############################################################################
# FINAL SETUP
###############################################################################


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/vflib2-080701/Makefile])
AC_CONFIG_FILES([src/sgd-1.3/Makefile])
AC_CONFIG_FILES([src/nspdk/Makefile])
AC_CONFIG_FILES([src/sgm/Makefile])
AC_CONFIG_FILES([src/sgm/HashMap.hh])
AC_CONFIG_FILES([src/ggl/Makefile])
AC_CONFIG_FILES([src/bin/Makefile])
AC_CONFIG_FILES([perl/Makefile])
AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([misc/Makefile])
AC_CONFIG_FILES([misc/doxygen/Makefile])
AC_CONFIG_FILES([src/doxygen.hh])
AC_CONFIG_FILES([misc/logo/Makefile])
AC_CONFIG_FILES([misc/literature/Makefile])
AC_CONFIG_FILES([misc/manual/Makefile])
AC_CONFIG_FILES([misc/manual/faq.tex])
AC_CONFIG_FILES([misc/manual/tutorial-rules.tex])
AC_CONFIG_FILES([misc/manual/tutorial-groups.tex])
AC_CONFIG_FILES([misc/manual/tutorial-viz.tex])
AC_CONFIG_FILES([data/Makefile])
AC_CONFIG_FILES([data/kegg/Makefile])
AC_CONFIG_FILES([data/aromModel/Makefile])
AC_CONFIG_FILES([ggl.pc])

# link manual source files to enable custom build path setup
AC_CONFIG_LINKS([misc/manual/mhchem.sty:misc/manual/mhchem.sty])

AC_OUTPUT


###############################################################################
# SETUP OUTPUT
###############################################################################

echo ""
echo "----------------------------------"
echo " Build Settings for ${PACKAGE_TARNAME} v${PACKAGE_VERSION}"
echo "----------------------------------"
echo "  host type            = $host_os"
echo "  install dir          = $prefix"
echo "  CC                   = $CC"
echo "  CFLAGS               = $AM_CPPFLAGS $CFLAGS"
echo "  CXX                  = $CXX"
echo "  CXXFLAGS             = $AM_CXXFLAGS $CXXFLAGS"
echo "  LDFLAGS              = $AM_LDFLAGS $LDFLAGS"
echo "  LIBS                 = $LIBS"
echo "----------------------------------"
echo " Compiled modules :"
echo "----------------------------------"
if test "$FOUND_OPENBABEL_2" = "1"; then
	echo "  OpenBabel module     : enabled"
else
	echo "  OpenBabel module     : disabled"
fi
echo "----------------------------------"
echo ""
echo " You can run ./make now."
echo ""
