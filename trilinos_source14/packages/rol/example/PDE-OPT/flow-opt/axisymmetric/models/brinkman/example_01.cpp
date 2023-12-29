// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Brinkman porosity optimization problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Solver.hpp"
#include "ROL_SingletonVector.hpp"
#include "ROL_ConstraintFromObjective.hpp"

#include "../../../../TOOLS/meshreader.hpp"
#include "../../../../TOOLS/pdeconstraint.hpp"
#include "../../../../TOOLS/pdeobjective.hpp"
#include "../../../../TOOLS/pdevector.hpp"
#include "../../../../TOOLS/integralconstraint.hpp"

#include "pde_brinkman.hpp"
#include "obj_brinkman.hpp"

//#include "fenv.h"

typedef double RealT;

int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_OVERFLOW);
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = ROL::makePtrFromRef(std::cout);
  }
  else {
    outStream = ROL::makePtrFromRef(bhs);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    bool output = parlist->sublist("SimOpt").sublist("Solve").get("Output Iteration History",false);
    output = (iprint > 0) && (myRank==0) && output;
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History",output);

    ROL::Ptr<MeshManager<RealT>>            meshMgr;
    ROL::Ptr<PDE<RealT>>                    pde;
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> con;
    ROL::Ptr<PDE_Constraint<RealT>>         pdecon;
    ROL::Ptr<Assembler<RealT>>              assembler;
    /*** Initialize main data structure. ***/
    meshMgr = ROL::makePtr<MeshReader<RealT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    pde = ROL::makePtr<PDE_NavierStokes<RealT>>(*parlist);
    con = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    pdecon = ROL::dynamicPtrCast<PDE_Constraint<RealT>>(con);
    assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    // Create state vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<>>         u_ptr, p_ptr, z_ptr, r_ptr;
    ROL::Ptr<std::vector<RealT>>            z0_ptr;
    ROL::Ptr<ROL::StdVector<RealT>>         z0p;
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> z1p;
    ROL::Ptr<ROL::Vector<RealT>>            up, pp, zp, rp, xp;
    bool useParamVar = parlist->sublist("Problem").get("Use Optimal Constant Velocity",false);
    int dim = 2;
    u_ptr = assembler->createStateVector();    u_ptr->randomize();
    p_ptr = assembler->createStateVector();    p_ptr->randomize();
    z_ptr = assembler->createControlVector();  z_ptr->randomize();
    r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    if (useParamVar) {
      z0_ptr = ROL::makePtr<std::vector<RealT>>(dim);
      z0p = ROL::makePtr<ROL::StdVector<RealT>>(z0_ptr);
      z1p = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    }
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    if (useParamVar)
      zp = ROL::makePtr<PDE_OptVector<RealT>>(z1p,z0p,myRank);
    else
      zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    xp = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);

    // Initialize quadratic objective function.
    ROL::Ptr<QoI<RealT>>                           qoi;
    ROL::Ptr<ROL::Objective_SimOpt<RealT>>         obj;
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>> robj;
    //qoi = ROL::makePtr<QoI_Velocity_NavierStokes<RealT>>(*parlist,
    //        ROL::staticPtrCast<PDE_NavierStokes<RealT>>(pde)->getVelocityFE(),
    //        ROL::staticPtrCast<PDE_NavierStokes<RealT>>(pde)->getPressureFE(),
    //        ROL::staticPtrCast<PDE_NavierStokes<RealT>>(pde)->getVelocityBdryFE(4),
    //        ROL::staticPtrCast<PDE_NavierStokes<RealT>>(pde)->getBdryCellLocIds(4),
    //        ROL::staticPtrCast<PDE_NavierStokes<RealT>>(pde)->getStateFieldInfo());
    qoi = ROL::makePtr<QoI_VelocityTracking_NavierStokes<RealT>>(*parlist,
            ROL::staticPtrCast<PDE_NavierStokes<RealT>>(pde)->getVelocityFE(),
            ROL::staticPtrCast<PDE_NavierStokes<RealT>>(pde)->getPressureFE(),
            ROL::staticPtrCast<PDE_NavierStokes<RealT>>(pde)->getStateFieldInfo());
    obj = ROL::makePtr<PDE_Objective<RealT>>(qoi,assembler);
    robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,con,up,zp,pp,true,false);

    // Build bound constraint
    ROL::Ptr<ROL::Vector<RealT>>          lp, hp;
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd;
    if (useParamVar) {
      ROL::Ptr<ROL::StdVector<RealT>>         l0p, h0p;
      ROL::Ptr<Tpetra::MultiVector<>>         l1_ptr, h1_ptr;
      ROL::Ptr<ROL::TpetraMultiVector<RealT>> l1p, h1p;
      l0p = ROL::makePtr<ROL::StdVector<RealT>>(dim,ROL::ROL_NINF<RealT>());
      h0p = ROL::makePtr<ROL::StdVector<RealT>>(dim,ROL::ROL_INF<RealT>());
      l1_ptr = assembler->createControlVector();
      h1_ptr = assembler->createControlVector();
      l1p = ROL::makePtr<PDE_PrimalOptVector<RealT>>(l1_ptr,pde,assembler,*parlist);
      h1p = ROL::makePtr<PDE_PrimalOptVector<RealT>>(h1_ptr,pde,assembler,*parlist);
      l1p->setScalar(0.0);
      h1p->setScalar(1.0);
      lp = ROL::makePtr<PDE_OptVector<RealT>>(l1p,l0p,myRank);
      hp = ROL::makePtr<PDE_OptVector<RealT>>(h1p,h0p,myRank);
    }
    else {
      lp = zp->clone(); lp->setScalar(0.0);
      hp = zp->clone(); hp->setScalar(1.0);
    }
    bnd = ROL::makePtr<ROL::Bounds<RealT>>(lp, hp);
    // Build optimization problem
    ROL::Ptr<ROL::Problem<RealT>> optProb;
    optProb = ROL::makePtr<ROL::Problem<RealT>>(robj, zp);
    optProb->addBoundConstraint(bnd);
    optProb->finalize(false,true,*outStream);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      ROL::Ptr<ROL::Vector<RealT>> rup = up->clone(); rup->randomize(-1.0,1.0);
      ROL::Ptr<ROL::Vector<RealT>> rzp = zp->clone(); rzp->randomize( 0.0,1.0);
      ROL::Ptr<ROL::Vector<RealT>> rpp = pp->clone(); rpp->randomize(-1.0,1.0);
      ROL::Ptr<ROL::Vector<RealT>> dup = up->clone(); dup->randomize(-1.0,1.0);
      ROL::Ptr<ROL::Vector<RealT>> dzp = zp->clone(); dzp->randomize( 0.0,1.0);
      con->checkApplyJacobian_1(*rup,*rzp,*dup,*rup,true,*outStream);
      con->checkApplyJacobian_2(*rup,*rzp,*dzp,*rup,true,*outStream);
      con->checkInverseJacobian_1(*rup,*rup,*rup,*rzp,true,*outStream);
      con->checkInverseAdjointJacobian_1(*rup,*rup,*rup,*rzp,true,*outStream);
      con->checkApplyAdjointHessian_11(*rup,*rzp,*rpp,*dup,*rup,true,*outStream);
      con->checkApplyAdjointHessian_21(*rup,*rzp,*rpp,*dzp,*rup,true,*outStream);
      con->checkApplyAdjointHessian_12(*rup,*rzp,*rpp,*dup,*rzp,true,*outStream);
      con->checkApplyAdjointHessian_22(*rup,*rzp,*rpp,*dzp,*rzp,true,*outStream);
      obj->checkGradient_1(*rup,*rzp,*dup,true,*outStream);
      obj->checkGradient_2(*rup,*rzp,*dzp,true,*outStream);
      obj->checkHessVec_11(*rup,*rzp,*dup,true,*outStream);
      obj->checkHessVec_12(*rup,*rzp,*dzp,true,*outStream);
      obj->checkHessVec_21(*rup,*rzp,*dup,true,*outStream);
      obj->checkHessVec_22(*rup,*rzp,*dzp,true,*outStream);
      robj->checkGradient(*rzp,*dzp,true,*outStream);
      robj->checkHessVec(*rzp,*dzp,true,*outStream);
      //optProb->check(*outStream);
    }

    // Solve optimization problem
    zp->setScalar(0.5);
    up->zero(); pp->zero();
    bool opt = parlist->sublist("Problem").get("Solve Optimization Problem",true);
    if (opt) {
      std::ifstream infile("control.txt");
      if (infile.good())
        assembler->inputTpetraVector(z_ptr,"control.txt");
      if (useParamVar) {
        std::ifstream infile0; infile0.open("target.txt");
        if (infile0.good()) {
          for (int i = 0; i < dim; ++i) infile0 >> (*z0_ptr)[i];
          infile0.close();
        }
      }
      ROL::Solver<RealT> optSolver(optProb, *parlist);
      optSolver.solve(*outStream);
      pdecon->outputTpetraVector(z_ptr,"control.txt");
      if (useParamVar) {
        std::ofstream outfile; outfile.open("target.txt");
        for (const auto x : *z0_ptr) {
          outfile << std::scientific << std::setprecision(16);
          outfile << x << std::endl;
        }
        outfile.close();
      }
    }

    // Output
    assembler->printMeshData(*outStream);
    RealT tol(1.e-8);
    Teuchos::Array<RealT> res(1,0);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    con->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    assembler->printDataPDE(pde,u_ptr,z_ptr);
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
    //pdecon->outputTpetraData();

    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
