# Copyright (C) 2001-2024 Quantum ESPRESSO group
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

AC_INIT(ESPRESSO, 7.4, , espresso)

AC_PREREQ(2.64)
AC_CONFIG_MACRO_DIR([m4/])

# Initialize variables, filling with what comes from the environment...
X_AC_QE_DEFAULT_ENV()

# configure for current directory by default
if test "$topdir" = "" ; then topdir="`pwd`" ; fi

# check system type (no cross-compilation for now)
AC_CANONICAL_BUILD

# Checking Architecture...
X_AC_QE_ARCH()

# Add all needed -D options to try_dflags
try_dflags=""

# Add needed include directories
try_iflags="-I. -I\$(TOPDIR)/include"

# Checking archiver...
X_AC_QE_AR()

# Checking OpenMP...
X_AC_QE_OPENMP()

AC_LANG_PUSH(Fortran) 
AC_FC_SRCEXT(f90) 

# Checking MPIF90...
X_AC_QE_MPIF90()

# Checking CC...
X_AC_QE_CC()

# Checking F90...
X_AC_QE_F90()

# Checking preprocessor...
X_AC_QE_CPP()

# Checking linker...
X_AC_QE_LD()

# Checking CUDA...
X_AC_QE_CUDA()

# Checking F90 rule...
X_AC_QE_F90RULE()

FFLAGS="$test_fflags"
LDFLAGS="$test_ldflags"

# Checking BLAS...
X_AC_QE_BLAS()

# Checking LAPACK...
X_AC_QE_LAPACK()

# Checking for FFT...
X_AC_QE_FFT()

# Setting libxc...
ACX_LIBXC()

# Setting FoX...
X_AC_QE_FOX()

# Checking for IBM MASS library...
X_AC_QE_MASS()

# check for MPI library...
X_AC_QE_MPI()

# Setting ScaLAPACK...
X_AC_QE_SCALAPACK()

# Setting ELPA...
X_AC_QE_ELPA()

# Setting HDF5...
X_AC_QE_HDF5()

# Checking environ...
X_AC_QE_ENVIRON()

# Checking for optional compilation macros
X_AC_QE_PLUGINS() 

# Checking SIGNAL...
X_AC_QE_SIGNAL()

# Checking STATUS...
X_AC_QE_STATUS()

# Checking for ranlib...
X_AC_QE_RANLIB()

# Checking wget or curl...
X_AC_QE_WGET()

# set up git
X_AC_QE_GIT()

# Show dflags before adding $(MANUAL_DFLAGS) and adapt to XLF (if needed)
if test "$dflags" = "" ; then dflags="$try_dflags" ; fi
echo setting DFLAGS... $try_dflags

# xlf compilers (AIX and powerpc) want comma-separated -D directives
X_AC_QE_AIX_DFLAGS()

if test "$iflags" = "" ; then iflags="$try_iflags" ; fi
echo setting IFLAGS... $iflags

# export additional settings to generated files
AC_SUBST(dflags)
AC_SUBST(fdflags)
AC_SUBST(iflags)
AC_SUBST(ld_libs)
AC_SUBST(topdir)
AC_SUBST(extlib_flags)

AC_CONFIG_FILES([include/configure.h:include/configure.h.in])

AC_CONFIG_FILES(make.inc)
AC_CONFIG_FILES(configure.msg)
AC_CONFIG_FILES(install/make_wannier90.inc)


AC_OUTPUT

# final messages
sed '/@delete@/d' configure.msg
echo configure: success
