// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Thyra_VectorStdOpsTester.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"

/** \brief Main test driver function for composite product spaces
 */
template <class Scalar>
bool run_product_space_tests(
  const int                                                     n
  ,const int                                                    numBlocks
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &tol
  ,const bool                                                   showAllTests
  ,const bool                                                   dumpAll
  ,Teuchos::FancyOStream                                        *out_arg
  )
{

  using Thyra::relErr;
  using Teuchos::OSTab;
  using Teuchos::rcp;
  using Teuchos::RCP;

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType    ScalarMag;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::fancyOStream(rcp(out_arg,false));

  if(out.get()) *out << "\n*** Entering run_product_space_tests<"<<ST::name()<<">(...) ...\n";

  bool success = true, result;

  Thyra::VectorSpaceTester<Scalar> vectorSpaceTester;
  vectorSpaceTester.warning_tol(ScalarMag(0.1)*tol);
  vectorSpaceTester.error_tol(tol);
  vectorSpaceTester.show_all_tests(showAllTests);
  vectorSpaceTester.dump_all(dumpAll);

  Thyra::VectorStdOpsTester<Scalar> vectorStdOpsTester;
  vectorStdOpsTester.warning_tol(ScalarMag(0.1)*tol);
  vectorStdOpsTester.error_tol(tol);

  Teuchos::Array<RCP<const Thyra::VectorSpaceBase<Scalar> > >
    vecSpaces(numBlocks);
  const RCP<const Teuchos::Comm<Thyra::Ordinal> >
    comm = Teuchos::DefaultComm<Thyra::Ordinal>::getComm();
  const int numProcs = size(*comm);
  RCP<const Thyra::VectorSpaceBase<Scalar> >
    spaceBlock = Thyra::defaultSpmdVectorSpace<Scalar>(comm,n,-1);
  for( int i = 0; i < numBlocks; ++i )
    vecSpaces[i] = spaceBlock;
  
  if(out.get()) *out << "\nA) Performing basic tests on product vectors with SPMD constituent vectors ...\n";

  if(out.get()) *out << "\nCreating a product space ps with numBlocks="<<numBlocks<<" and n="<<n<<"vector elements per block ...\n";

  RCP<Thyra::DefaultProductVectorSpace<Scalar> > ps =
    Thyra::productVectorSpace<Scalar>(vecSpaces());

  if(out.get()) *out << "\nps->numBlocks()=";
  result = ps->numBlocks() == numBlocks;
  if(!result) success = false;
  if(out.get()) *out
    << ps->numBlocks() << " == numBlocks=" << numBlocks
    << " : " << ( result ? "passed" : "failed" ) << std::endl;

  if(out.get()) *out << "\nTesting the product space ps ...\n";

  if(out.get()) *out << "\nps->dim()=";
  result = ps->dim() == numProcs*n*numBlocks;
  if(!result) success = false;
  if(out.get()) *out
    << ps->dim() << " == numProcs*n*numBlocks=" << numProcs*n*numBlocks
    << " : " << ( result ? "passed" : "failed" ) << std::endl;
  
  if(out.get()) *out << "\nTesting the VectorSpaceBase interface of ps ...\n";
  TEUCHOS_TEST_ASSERT(vectorSpaceTester.check(*ps, out.get()), *out, success);
  
  if(out.get()) *out << "\nTesting standard vector ops for ps ...\n";
  TEUCHOS_TEST_ASSERT(vectorStdOpsTester.checkStdOps(*ps, out.get()), *out, success);
  
  if(out.get()) *out << "\nB) Testing a nested product space of product vector spaces called pps ...\n";

  Teuchos::Array<RCP<const Thyra::VectorSpaceBase<Scalar> > >
    blockVecSpaces(numBlocks);
  for( int i = 0; i < numBlocks; ++i )
    blockVecSpaces[i] = ps;

  RCP<Thyra::DefaultProductVectorSpace<Scalar> > pps =
    Thyra::productVectorSpace<Scalar>(blockVecSpaces());
  
  if(out.get()) *out << "\nTesting the VectorSpaceBase interface of pps ...\n";
  TEUCHOS_TEST_ASSERT(vectorSpaceTester.check(*pps, out.get()), *out, success);
  
  if(out.get()) *out << "\nTesting standard vector ops for pps ...\n";
  TEUCHOS_TEST_ASSERT(vectorStdOpsTester.checkStdOps(*pps, out.get()), *out, success);
    
  if(out.get()) *out
    << "\n*** Leaving run_product_space_tests<"<<ST::name()<<">(...) ...\n";
  
  return success;

} // end run_product_space_tests() [Doxygen looks for this!]


int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;
  using Teuchos::RCP;

  bool success = true;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from command-line
    //

    int n              = 4;
    int numBlocks      = 4;
    bool showAllTests  = true;
    bool dumpAll       = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "n", &n, "Number of elements in each constituent vector." );
    clp.setOption( "num-blocks", &numBlocks, "blocks to create." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests, "Determines if all tests are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    //
    // Run the tests
    //

#ifdef HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_product_space_tests<float>(n,numBlocks,float(1e-5),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#endif // HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_product_space_tests<double>(n,numBlocks,double(1e-13),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#if defined(HAVE_THYRA_COMPLEX)
#ifdef THYRA_TEUCHOS_BLASFLOAT
    if( !run_product_space_tests<std::complex<float> >(n,numBlocks,float(1e-5),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#endif // HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_product_space_tests<std::complex<double> >(n,numBlocks,double(1e-12),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#endif // defined(HAVE_THYRA_COMPLEX)
#ifdef HAVE_TEUCHOS_GNU_MP
    //if( !run_product_space_tests<mpf_class>(n,numBlocks,mpf_class(1e-13),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
    // Above commented out code will not compile because its ScalarTraits specialization does not support eps()
#endif // HAVE_TEUCHOS_GNU_MP

  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

  if(verbose) {
    if(success)
      *out << "\nAll of the tests seem to have run successfully!\n";
    else
      *out << "\nOh no! at least one of the test failed!\n";	
  }
  
  return success ? 0 : 1;

} // end main() [Doxygen looks for this!]
