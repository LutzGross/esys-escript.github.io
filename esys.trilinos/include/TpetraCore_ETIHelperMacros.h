#ifndef TPETRACORE_ETIHELPERMACROS_H
#define TPETRACORE_ETIHELPERMACROS_H

#include <Tpetra_ConfigDefs.hpp>

/* Tpetra provides official support for the following nodes */
#include "Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp"

/* Tpetra provides official support for dd_real and qd_real */
#if defined(HAVE_TPETRA_QD)
#include <qd/qd_real.h>
#endif

#define TPETRA_INSTANTIATE_SLGN( INSTMACRO ) \
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( INSTMACRO ) \
  TPETRA_INSTANTIATE_GLGN( INSTMACRO )

#define TPETRA_INSTANTIATE_VECTOR( INSTMACRO ) \
  TPETRA_INSTANTIATE_SLGN( INSTMACRO )

#define TPETRA_INSTANTIATE_MULTIVECTOR( INSTMACRO ) \
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( INSTMACRO ) \
  TPETRA_INSTANTIATE_GLGN_NO_INT( INSTMACRO ) \
  TPETRA_INSTANTIATE_LLGN( INSTMACRO )

#define TPETRA_INSTANTIATE_LLGN(INSTMACRO)\
	INSTMACRO( int , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( int , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define TPETRA_INSTANTIATE_GLGN_NO_INT(INSTMACRO)


#define TPETRA_INSTANTIATE_GLGN(INSTMACRO)\
	INSTMACRO( int , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( int , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define TPETRA_INSTANTIATE_GLG(INSTMACRO)\
	INSTMACRO( int , int , int )


#define TPETRA_INSTANTIATE_PLGN(INSTMACRO)\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )\
	INSTMACRO( std_complex0double0 , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( std_complex0double0 , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )\
	INSTMACRO( int , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( int , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define TPETRA_INSTANTIATE_LGN(INSTMACRO)\
	INSTMACRO( int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define TPETRA_INSTANTIATE_SLG( INSTMACRO ) \
  TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR( INSTMACRO ) \
  TPETRA_INSTANTIATE_GLG( INSTMACRO )

#define TPETRA_INSTANTIATE_LG(INSTMACRO)\
	INSTMACRO( int , int )


#define TPETRA_INSTANTIATE_SL(INSTMACRO)\
	INSTMACRO( double , int )\
	INSTMACRO( std_complex0double0 , int )\
	INSTMACRO( int , int )


#define TPETRA_INSTANTIATE_SN(INSTMACRO)\
	INSTMACRO( double , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( double , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )\
	INSTMACRO( std_complex0double0 , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( std_complex0double0 , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )\
	INSTMACRO( int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define TPETRA_INSTANTIATE_GN(INSTMACRO)\
	INSTMACRO( int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define TPETRA_INSTANTIATE_S(INSTMACRO)\
	INSTMACRO( double )\
	INSTMACRO( std_complex0double0 )\
	INSTMACRO( int )


#define TPETRA_INSTANTIATE_L(INSTMACRO)\
	INSTMACRO( int )


#define TPETRA_INSTANTIATE_G(INSTMACRO)\
	INSTMACRO( int )


#define TPETRA_INSTANTIATE_N(INSTMACRO)\
	INSTMACRO( Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define TPETRA_INSTANTIATE_TSLGN(INSTMACRO)


#define TPETRA_INSTANTIATE_TSLG(INSTMACRO)


#define TPETRA_INSTANTIATE_CONVERT(INSTMACRO)


#define TPETRA_INSTANTIATE_CONVERT_SSL(INSTMACRO)


#define TPETRA_INSTANTIATE_TESTMV( INSTMACRO ) \
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( INSTMACRO )

#define TPETRA_INSTANTIATE_DOUBLE_INT_INT_N(INSTMACRO)\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(INSTMACRO)\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )\
	INSTMACRO( std_complex0double0 , int , int , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( std_complex0double0 , int , int , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR(INSTMACRO)\
	INSTMACRO( double , int , int )\
	INSTMACRO( std_complex0double0 , int , int )


#define TPETRA_INSTANTIATE_SL_NO_ORDINAL_SCALAR(INSTMACRO)\
	INSTMACRO( double , int )\
	INSTMACRO( std_complex0double0 , int )


#define TPETRA_INSTANTIATE_SN_NO_ORDINAL_SCALAR(INSTMACRO)\
	INSTMACRO( double , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( double , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )\
	INSTMACRO( std_complex0double0 , Tpetra_KokkosCompat_KokkosSerialWrapperNode )\
	INSTMACRO( std_complex0double0 , Tpetra_KokkosCompat_KokkosOpenMPWrapperNode )


#define TPETRA_INSTANTIATE_S_NO_ORDINAL_SCALAR(INSTMACRO)\
	INSTMACRO( double )\
	INSTMACRO( std_complex0double0 )


#define TPETRA_ETI_MANGLING_TYPEDEFS()  \
	typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode Tpetra_KokkosCompat_KokkosSerialWrapperNode; \
	typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Tpetra_KokkosCompat_KokkosOpenMPWrapperNode; \
	typedef std::complex<double> std_complex0double0;

#endif // TPETRACORE_ETIHELPERMACROS_H
