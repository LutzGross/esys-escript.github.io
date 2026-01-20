
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "SystemMatrixTestCase.h"

#include <ripley/Rectangle.h>
#include <ripley/RipleySystemMatrix.h>

#include <escript/FunctionSpaceFactory.h>

#include <cppunit/TestCaller.h>

using namespace CppUnit;
using namespace std;

// number of matrix rows (for blocksize 1)
const int rows = 20;

// diagonal offsets for full matrix
const int diag_off[] = { -6, -5, -4, -1, 0, 1, 4, 5, 6 };

// to get these reference values save the matrix via saveMM(), then in python:
// print scipy.io.mmread('mat.mtx') * range(20*blocksize)

// reference results - non-symmetric
// block size 1
const double ref_bs1[] =
    {   5800.,   6400.,  28200.,  29200.,  72100.,  73950.,
      145600., 148200., 251200., 254600., 388800., 393000.,
      558400., 563400., 608000., 612800., 360600., 363400.,
      475800., 479000.};
// block size 2
const double ref_bs2[] =
    {   24455.,   24507.,   27053.,   27105.,  118083.,  118167.,
       122281.,  122365.,  299244.,  299399.,  306990.,  307145.,
       601398.,  601614.,  612194.,  612410., 1031854., 1032134.,
      1045850., 1046130., 1590310., 1590654., 1607506., 1607850.,
      2276766., 2277174., 2297162., 2297570., 2475540., 2475931.,
      2495086., 2495477., 1468199., 1468427., 1479597., 1479825.,
      1933027., 1933287., 1946025., 1946285.};
// block size 3
const double ref_bs3[] =
    {  56174.,    56294.,    56414.,    62156.,    62276.,    62396.,
      269954.,   270146.,   270338.,   279536.,   279728.,   279920.,
      681940.,   682294.,   682648.,   699604.,   699958.,   700312.,
     1368088.,  1368580.,  1369072.,  1392652.,  1393144.,  1393636.,
     2342848.,  2343484.,  2344120.,  2374612.,  2375248.,  2375884.,
     3605608.,  3606388.,  3607168.,  3644572.,  3645352.,  3646132.,
     5156368.,  5157292.,  5158216.,  5202532.,  5203456.,  5204380.,
     5603773.,  5604658.,  5605543.,  5647987.,  5648872.,  5649757.,
     3323474.,  3323990.,  3324506.,  3349256.,  3349772.,  3350288.,
     4372454.,  4373042.,  4373630.,  4401836.,  4402424.,  4403012.};
// block size 4
const double ref_bs4[] =
    { 101334.,    101550.,    101766.,    101982.,    112062.,
      112278.,    112494.,    112710.,    484358.,    484702.,
      485046.,    485390.,    501486.,    501830.,    502174.,
      502518.,   1221084.,   1221718.,   1222352.,   1222986.,
     1252640.,   1253274.,   1253908.,   1254542.,   2446892.,
     2447772.,   2448652.,   2449532.,   2490748.,   2491628.,
     2492508.,   2493388.,   4185740.,   4186876.,   4188012.,
     4189148.,   4242396.,   4243532.,   4244668.,   4245804.,
     6436588.,   6437980.,   6439372.,   6440764.,   6506044.,
     6507436.,   6508828.,   6510220.,   9199436.,   9201084.,
     9202732.,   9204380.,   9281692.,   9283340.,   9284988.,
     9286636.,   9994708.,   9996286.,   9997864.,   9999442.,
    10073464.,  10075042.,  10076620.,  10078198.,   5927606.,
     5928526.,   5929446.,   5930366.,   5973534.,   5974454.,
     5975374.,   5976294.,   7795430.,   7796478.,   7797526.,
     7798574.,   7847758.,   7848806.,   7849854.,   7850902.};

const double* ref[] = {ref_bs1, ref_bs2, ref_bs3, ref_bs4};

// reference results - symmetric
// block size 1
const double ref_symm_bs1[] =
    {  5800.,    6400.,   28200.,   29500.,   72100.,   75600.,
     142550.,  153300.,  243850.,  263000.,  377150.,  404700.,
     542450.,  578400.,  587750.,  631100.,  336050.,  382600.,
     446950.,  501200.};
// block size 2
const double ref_symm_bs2[] =
    {  24455.,    24507.,    27202.,    27255.,   118083.,   118167.,
      123626.,   123719.,   298942.,   299096.,   314374.,   314574.,
      588036.,   588250.,   633214.,   633510.,  1001284.,  1001578.,
     1080054.,  1080446.,  1542532.,  1542906.,  1654894.,  1655382.,
     2211780.,  2212234.,  2357734.,  2358318.,  2393346.,  2393799.,
     2568842.,  2569441.,  1368797.,  1369103.,  1556820.,  1557223.,
     1816417.,  1816771.,  2035236.,  2035695.};
// block size 3
const double ref_symm_bs3[] =
    {  56174.,    56294.,    56414.,    62596.,    62722.,    62848.,
      269954.,   270146.,   270338.,   282640.,   282874.,   283108.,
      681025.,   681376.,   681727.,   716694.,   717246.,   717798.,
     1337096.,  1337600.,  1338104.,  1440210.,  1441050.,  1441890.,
     2273084.,  2273804.,  2274524.,  2451726.,  2452854.,  2453982.,
     3497072.,  3498008.,  3498944.,  3751242.,  3752658.,  3754074.,
     5009060.,  5010212.,  5011364.,  5338758.,  5340462.,  5342166.,
     5417693.,  5418878.,  5420063.,  5813769.,  5815578.,  5817387.,
     3098622.,  3099510.,  3100398.,  3522842.,  3524132.,  3525422.,
     4108830.,  4109862.,  4110894.,  4602314.,  4603784.,  4605254.};
// block size 4
const double ref_symm_bs4[] =
    { 101334.,    101550.,    101766.,    101982.,    112920.,
      113154.,    113388.,    113622.,    484358.,    484702.,
      485046.,    485390.,    507000.,    507458.,    507916.,
      508374.,   1219228.,   1219856.,   1220484.,   1221112.,
     1283176.,   1284336.,   1285496.,   1286656.,   2390844.,
     2391776.,   2392708.,   2393640.,   2575048.,   2576848.,
     2578648.,   2580448.,   4060604.,   4061984.,   4063364.,
     4064744.,   4378920.,   4381360.,   4383800.,   4386240.,
     6242364.,   6244192.,   6246020.,   6247848.,   6694792.,
     6697872.,   6700952.,   6704032.,   8936124.,   8938400.,
     8940676.,   8942952.,   9522664.,   9526384.,   9530104.,
     9533824.,   9662308.,   9664706.,   9667104.,   9669502.,
    10366660.,  10370694.,  10374728.,  10378762.,   5526118.,
     5528050.,   5529982.,   5531914.,   6280848.,   6283822.,
     6286796.,   6289770.,   7324854.,   7327106.,   7329358.,
     7331610.,   8202640.,   8206030.,   8209420.,   8212810.};

const double* ref_symm[] = {ref_symm_bs1, ref_symm_bs2, ref_symm_bs3, ref_symm_bs4};

/// helper
double lsup(const double* d0, const double* d1, int length)
{
    double result = 0.;
    for (int i=0; i<length; i++) {
        result = std::max(result, std::abs(d0[i] - d1[i]));
        //std::cerr << d0[i] << " " << d1[i] << std::endl;
    }
    return result;
}

TestSuite* SystemMatrixTestCase::suite()
{
    TestSuite *testSuite = new TestSuite("SystemMatrixTestCase");
    testSuite->addTest(new TestCaller<SystemMatrixTestCase>(
                "testSpMV_CPU_blocksize1_nonsymmetric",
                &SystemMatrixTestCase::testSpMV_CPU_blocksize1_nonsymmetric));
    testSuite->addTest(new TestCaller<SystemMatrixTestCase>(
                "testSpMV_CPU_blocksize2_nonsymmetric",
                &SystemMatrixTestCase::testSpMV_CPU_blocksize2_nonsymmetric));
    testSuite->addTest(new TestCaller<SystemMatrixTestCase>(
                "testSpMV_CPU_blocksize3_nonsymmetric",
                &SystemMatrixTestCase::testSpMV_CPU_blocksize3_nonsymmetric));
    testSuite->addTest(new TestCaller<SystemMatrixTestCase>(
                "testSpMV_CPU_blocksize4_nonsymmetric",
                &SystemMatrixTestCase::testSpMV_CPU_blocksize4_nonsymmetric));
    testSuite->addTest(new TestCaller<SystemMatrixTestCase>(
                "testSpMV_CPU_blocksize1_symmetric",
                &SystemMatrixTestCase::testSpMV_CPU_blocksize1_symmetric));
    testSuite->addTest(new TestCaller<SystemMatrixTestCase>(
                "testSpMV_CPU_blocksize2_symmetric",
                &SystemMatrixTestCase::testSpMV_CPU_blocksize2_symmetric));
    testSuite->addTest(new TestCaller<SystemMatrixTestCase>(
                "testSpMV_CPU_blocksize3_symmetric",
                &SystemMatrixTestCase::testSpMV_CPU_blocksize3_symmetric));
    testSuite->addTest(new TestCaller<SystemMatrixTestCase>(
                "testSpMV_CPU_blocksize4_symmetric",
                &SystemMatrixTestCase::testSpMV_CPU_blocksize4_symmetric));
    return testSuite;
}

void SystemMatrixTestCase::setUp()
{
    mpiInfo = escript::makeInfo(MPI_COMM_WORLD);
    domain.reset(new ripley::Rectangle(4, 3, 0., 0., 1., 1.));
}

escript::ASM_ptr SystemMatrixTestCase::createMatrix(int blocksize,
                                                    bool symmetric)
{
    escript::FunctionSpace fs(escript::solution(*domain));
    const int firstdiag = (symmetric ? 4 : 0);
    const ripley::IndexVector offsets(diag_off+firstdiag, diag_off+9);

    // create a matrix with 9 diagonals, given blocksize and symmetric flag
    escript::ASM_ptr matptr(new ripley::SystemMatrix(mpiInfo, blocksize, fs,
                             rows, offsets, symmetric));
    ripley::SystemMatrix* mat(dynamic_cast<ripley::SystemMatrix*>(matptr.get()));

    ripley::IndexVector rowIdx(4);
    std::vector<double> array(4*4*blocksize*blocksize);

    for (int i=0; i<8; i++) {
        rowIdx[0] = 2*i;
        rowIdx[1] = 2*i+1;
        rowIdx[2] = 2*i+4;
        rowIdx[3] = 2*i+5;
        for (int j=0; j<4*4; j++) {
            for (int k=0; k<blocksize; k++) {
                for (int l=0; l<blocksize; l++) {
                    // make main diagonal blocks symmetric since the current
                    // implementation actually reads full main diagonal blocks
                    // so if symmetric flag is set the matrix really has to be
                    // symmetric!
                    array[j*blocksize*blocksize + k*blocksize + l] = 
                        (j%5==0 ? 1000*i+50*j+k+l : 1000*i+50*j+blocksize*k+l);
                }
            }
        }
        mat->add(rowIdx, array);
    }
    //mat->saveMM("/tmp/test.mtx");
    return matptr;
}

escript::Data SystemMatrixTestCase::createInputVector(int blocksize)
{
    escript::FunctionSpace fs(escript::solution(*domain));
    escript::DataTypes::ShapeType shape;
    if (blocksize > 1)
        shape.push_back(blocksize);
    escript::Data x(0., shape, fs, true);
    for (int i=0; i<rows; i++) {
        double* xx= x.getSampleDataRW(i);
        for (int j=0; j<blocksize; j++)
            xx[j] = (double)i*blocksize + j;
    }
    return x;
}

void SystemMatrixTestCase::testSpMV_CPU_blocksize1_nonsymmetric()
{
    int blocksize = 1;
    bool symmetric = false;
    escript::ASM_ptr mat(createMatrix(blocksize, symmetric));
    const escript::Data x(createInputVector(blocksize));
    const escript::Data y(mat->vectorMultiply(x));
    const double* yref = ref[blocksize-1];
    const double* yy = y.getSampleDataRO(0);
    double error = lsup(yref, yy, blocksize*rows);
    CPPUNIT_ASSERT(error < 1e-12);
}

void SystemMatrixTestCase::testSpMV_CPU_blocksize2_nonsymmetric()
{
    int blocksize = 2;
    bool symmetric = false;
    escript::ASM_ptr mat(createMatrix(blocksize, symmetric));
    const escript::Data x(createInputVector(blocksize));
    const escript::Data y(mat->vectorMultiply(x));
    const double* yref = ref[blocksize-1];
    const double* yy = y.getSampleDataRO(0);
    double error = lsup(yref, yy, blocksize*rows);
    CPPUNIT_ASSERT(error < 1e-12);
}

void SystemMatrixTestCase::testSpMV_CPU_blocksize3_nonsymmetric()
{
    int blocksize = 3;
    bool symmetric = false;
    escript::ASM_ptr mat(createMatrix(blocksize, symmetric));
    const escript::Data x(createInputVector(blocksize));
    const escript::Data y(mat->vectorMultiply(x));
    const double* yref = ref[blocksize-1];
    const double* yy = y.getSampleDataRO(0);
    double error = lsup(yref, yy, blocksize*rows);
    CPPUNIT_ASSERT(error < 1e-12);
}

void SystemMatrixTestCase::testSpMV_CPU_blocksize4_nonsymmetric()
{
    int blocksize = 4;
    bool symmetric = false;
    escript::ASM_ptr mat(createMatrix(blocksize, symmetric));
    const escript::Data x(createInputVector(blocksize));
    const escript::Data y(mat->vectorMultiply(x));
    const double* yref = ref[blocksize-1];
    const double* yy = y.getSampleDataRO(0);
    double error = lsup(yref, yy, blocksize*rows);
    CPPUNIT_ASSERT(error < 1e-12);
}

void SystemMatrixTestCase::testSpMV_CPU_blocksize1_symmetric()
{
    int blocksize = 1;
    bool symmetric = true;
    escript::ASM_ptr mat(createMatrix(blocksize, symmetric));
    const escript::Data x(createInputVector(blocksize));
    const escript::Data y(mat->vectorMultiply(x));
    const double* yref = ref_symm[blocksize-1];
    const double* yy = y.getSampleDataRO(0);
    double error = lsup(yref, yy, blocksize*rows);
    CPPUNIT_ASSERT(error < 1e-12);
}

void SystemMatrixTestCase::testSpMV_CPU_blocksize2_symmetric()
{
    int blocksize = 2;
    bool symmetric = true;
    escript::ASM_ptr mat(createMatrix(blocksize, symmetric));
    const escript::Data x(createInputVector(blocksize));
    const escript::Data y(mat->vectorMultiply(x));
    const double* yref = ref_symm[blocksize-1];
    const double* yy = y.getSampleDataRO(0);
    double error = lsup(yref, yy, blocksize*rows);
    CPPUNIT_ASSERT(error < 1e-12);
}

void SystemMatrixTestCase::testSpMV_CPU_blocksize3_symmetric()
{
    int blocksize = 3;
    bool symmetric = true;
    escript::ASM_ptr mat(createMatrix(blocksize, symmetric));
    const escript::Data x(createInputVector(blocksize));
    const escript::Data y(mat->vectorMultiply(x));
    const double* yref = ref_symm[blocksize-1];
    const double* yy = y.getSampleDataRO(0);
    double error = lsup(yref, yy, blocksize*rows);
    CPPUNIT_ASSERT(error < 1e-12);
}

void SystemMatrixTestCase::testSpMV_CPU_blocksize4_symmetric()
{
    int blocksize = 4;
    bool symmetric = true;
    escript::ASM_ptr mat(createMatrix(blocksize, symmetric));
    const escript::Data x(createInputVector(blocksize));
    const escript::Data y(mat->vectorMultiply(x));
    const double* yref = ref_symm[blocksize-1];
    const double* yy = y.getSampleDataRO(0);
    double error = lsup(yref, yy, blocksize*rows);
    CPPUNIT_ASSERT(error < 1e-12);
}

