// $Id$
//==============================================================================
//!
//! \file main_AdvectionDiffusion.C
//!
//! \date 12 June 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Main program for the isogeometric Advection-Diffusion solver.
//!
//==============================================================================

#include "IFEM.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMExplicitRK.h"
#include "SIMExplicitRKE.h"
#include "SIMSolver.h"
#include "SIMSolverAdap.h"
#include "SIMAD.h"
#include "AppCommon.h"
#include "AdvectionDiffusionBDF.h"
#include "AdvectionDiffusionExplicit.h"
#include "AdaptiveSIM.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "Utilities.h"
#include "Profiler.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


template<template<class T> class Solver=SIMSolver, class AD>
int runSimulatorStationary(char* infile, AD& model)
{
  utl::profiler->start("Model input");
  Solver<AD> solver(model);

  // Read in model definitions
  if (!model.read(infile) || !solver.read(infile))
    return 1;

  model.opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  // Establish the FE data structures
  if (!model.preprocess())
    return 1;

  model.setQuadratureRule(model.opt.nGauss[0],true);
  model.setMode(SIM::STATIC);

  std::unique_ptr<DataExporter> exporter;
  if (model.opt.dumpHDF5(infile))
    if (model.opt.discretization < ASM::Spline && !model.opt.hdf5.empty())
      IFEM::cout <<"\n ** HDF5 output is available for spline discretization only"
        <<". Deactivating...\n"<< std::endl;
    else
      exporter.reset(SIM::handleDataOutput(model, solver, model.opt.hdf5, false, 1, 1));

  model.initSystem(model.opt.solver);
  model.init();

  return solver.solveProblem(infile, exporter.get(), "Solving Advection-Diffusion problem", false);
}


template<class Solver, class AD>
int runSimulatorTransientImpl(char* infile, TimeIntegration::Method tIt,
                              Solver& sim, AD& model)
{
  utl::profiler->start("Model input");

  SIMSolver<Solver> solver(sim);

  // Read in model definitions
  if (!model.read(infile) || !solver.read(infile))
    return 1;

  model.opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  // Establish the FE data structures
  if (!model.preprocess())
    return 1;

  model.setMode(SIM::DYNAMIC);
  model.initSystem(model.opt.solver,1,1);
  model.setQuadratureRule(model.opt.nGauss[0],true);

  std::unique_ptr<DataExporter> exporter;
  if (model.opt.dumpHDF5(infile))
    if (model.opt.discretization < ASM::Spline && !model.opt.hdf5.empty())
      IFEM::cout <<"\n ** HDF5 output is available for spline discretization only"
        <<". Deactivating...\n"<< std::endl;
    else
      exporter.reset(SIM::handleDataOutput(model, solver, model.opt.hdf5, false, 1, 1));

  model.init(solver.getTimePrm());

  if (solver.solveProblem(infile, exporter.get()))
    return 5;

  model.printFinalNorms(solver.getTimePrm());

  return 0;
}


template<class Dim>
int runSimulator(char* infile, bool adap,
                 TimeIntegration::Method tIt, int integrandType)
{
  if (tIt == TimeIntegration::NONE)  {
    AdvectionDiffusion integrand(Dim::dimension);
    SIMAD<Dim> model(integrand,true);
    if (adap)
      return runSimulatorStationary<SIMSolverAdap>(infile, model);
    else
      return runSimulatorStationary(infile, model);
  }
  else if (tIt == TimeIntegration::BE || tIt == TimeIntegration::BDF2) {
    AdvectionDiffusionBDF integrand(Dim::dimension,
                                    tIt==TimeIntegration::BE?1:2,
                                    integrandType);
    SIMAD<Dim,AdvectionDiffusionBDF> model(integrand, true);
    return runSimulatorTransientImpl(infile, tIt, model, model);
  }
  else {
    AdvectionDiffusionExplicit integrand(Dim::dimension, integrandType);
    typedef SIMAD<Dim,AdvectionDiffusionExplicit> ADSIM;
    ADSIM model(integrand, true);
    if (tIt >= TimeIntegration::HEUNEULER) {
      TimeIntegration::SIMExplicitRKE<ADSIM> sim(model, tIt);
      return runSimulatorTransientImpl(infile, tIt, sim, model);
    }
    else {
      TimeIntegration::SIMExplicitRK<ADSIM> sim(model, tIt);
      return runSimulatorTransientImpl(infile, tIt, sim, model);
    }
  }
}


/*!
  \brief Main program for the isogeometric Advection-Diffusion solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -2D : Use two-parametric simulation driver
  \arg -adap : Use adaptive simulation driver with LR-splines discretization
  \arg -BE : Time dependent (BE) simulation
  \arg -BDF2 : Time dependent (BDF2) simulation
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  int integrandType = Integrand::STANDARD;
  TimeIntegration::Method tInt = TimeIntegration::NONE;
  bool twoD = false;
  bool adap = false;
  char* infile = nullptr;

  IFEM::Init(argc,argv,"Advection-Diffusion solver");

  for (int i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-adap"))
      adap = true;
    else if (!strcmp(argv[i],"-2D"))
      twoD = true;
    else if (!strcmp(argv[i],"-be"))
      tInt = TimeIntegration::BE;
    else if (!strcmp(argv[i],"-bdf2"))
      tInt = TimeIntegration::BDF2;
    else if (!strcmp(argv[i],"-euler"))
      tInt = TimeIntegration::EULER;
    else if (!strcmp(argv[i],"-heun"))
      tInt = TimeIntegration::HEUN;
    else if (!strcmp(argv[i],"-heuneuler"))
      tInt = TimeIntegration::HEUNEULER;
    else if (!strcmp(argv[i],"-bs"))
      tInt = TimeIntegration::BOGACKISHAMPINE;
    else if (!strcmp(argv[i],"-fehlberg"))
      tInt = TimeIntegration::FEHLBERG;
    else if (!strcmp(argv[i],"-rk3"))
      tInt = TimeIntegration::RK3;
    else if (!strcmp(argv[i],"-rk4"))
      tInt = TimeIntegration::RK4;
    else if (!strcmp(argv[i],"-supg"))
      integrandType = Integrand::SECOND_DERIVATIVES | Integrand::G_MATRIX;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-2D] [-nGauss <n>] [-adap]\n"
	      <<"       [-hdf5] [-vtf <format> [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]]\n";
    return 0;
  }

  if (adap)
    IFEM::getOptions().discretization = ASM::LRSpline;

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout) << std::endl;
  utl::profiler->stop("Initialization");

  if (twoD)
    return runSimulator<SIM2D>(infile,adap,tInt,integrandType);
  else
    return runSimulator<SIM3D>(infile,adap,tInt,integrandType);
}
