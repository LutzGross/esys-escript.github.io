// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef IFPACK2_ETIHELPERMACROS_H_
#define IFPACK2_ETIHELPERMACROS_H_

#include <Ifpack2_ConfigDefs.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>

#define IFPACK2_INSTANTIATE_SL(INSTMACRO)\
	INSTMACRO( double , int )\
	INSTMACRO( std_complex0double0 , int )


#define IFPACK2_INSTANTIATE_SL_REAL(INSTMACRO)\
	INSTMACRO( double , int , int )


#define IFPACK2_INSTANTIATE_L(INSTMACRO)\
	INSTMACRO( int )


#define IFPACK2_INSTANTIATE_SLG(INSTMACRO)\
	INSTMACRO( double , int , int )\
	INSTMACRO( std_complex0double0 , int , int )


#define IFPACK2_INSTANTIATE_SLG_REAL(INSTMACRO)\
	INSTMACRO( double , int , int )


#define IFPACK2_INSTANTIATE_LG(INSTMACRO)\
	INSTMACRO( int , int )


#define IFPACK2_INSTANTIATE_SLGN(INSTMACRO)\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )\
	INSTMACRO( std_complex0double0 , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( std_complex0double0 , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define IFPACK2_INSTANTIATE_SLGN_REAL(INSTMACRO)\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define IFPACK2_INSTANTIATE_LGN(INSTMACRO)\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )\
	INSTMACRO( std_complex0double0 , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( std_complex0double0 , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define IFPACK2_INSTANTIATE_N(INSTMACRO)\
	INSTMACRO( Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define IFPACK2_ETI_MANGLING_TYPEDEFS()  \
	typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode Tpetra_KokkosCompat_KokkosSerialWrapperNode; \
	typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Tpetra_KokkosCompat_KokkosOpenMPWrapperNode; \
	typedef std::complex<double> std_complex0double0;

#endif // IFPACK2_ETIHELPERMACROS_H_
