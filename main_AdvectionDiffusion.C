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

#include "AdvectionDiffusion.h"
#include "AdvectionDiffusionArgs.h"
#include "AdvectionDiffusionBDF.h"
#include "AdvectionDiffusionExplicit.h"
#include "AdvectionDiffusionImplicit.h"
#include "SIMAD.h"

#include "ASMenums.h"
#include "IFEM.h"
#include "LogStream.h"
#include "Profiler.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMconfigure.h"
#include "SIMExplicitRK.h"
#include "SIMExplicitRKE.h"
#include "SIMExplicitLMM.h"
#include "SIMImplicitLMM.h"
#include "SIMoptions.h"
#include "SIMSolver.h"
#include "SIMSolverAdap.h"
#include "TimeStep.h"
#include "TimeIntUtils.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>


/*!
  \brief Runs a stationary advection-diffusion problem.
*/

template<template<class T> class Solver=SIMSolverStat, class AD>
int runSimulatorStationary (char* infile, AD& model, bool printNorms)
{
  utl::profiler->start("Model input");

  Solver<AD> solver(model);
  if constexpr (std::is_same<Solver<AD>, SIMSolverAdap<AD>>::value)
    solver.setExporterName("T");

  int res = ConfigureSIM(model, infile, typename AD::SetupProps());
  if (res)
    return res;

  // Read in model definitions
  if (!solver.read(infile))
    return 1;

  model.opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  if (model.opt.dumpHDF5(infile))
    solver.handleDataOutput(model.opt.hdf5,model.getProcessAdm());

  res = solver.solveProblem(infile,"Solving Advection-Diffusion problem");
  if (!res && printNorms) model.printFinalNorms(TimeStep());

  return res;
}


/*!
  \brief Runs a transient advection-diffusion problem.
*/

template<class Solver, class AD>
int runSimulatorTransientImpl (char* infile, Solver& sim, AD& model)
{
  utl::profiler->start("Model input");

  SIMSolver<Solver> solver(sim);

  int res = ConfigureSIM(model, infile, typename AD::SetupProps());
  if (res)
    return res;

  // Read in model definitions
  if (!solver.read(infile))
    return 1;

  model.opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  if (solver.restart(model.opt.restartFile,model.opt.restartStep) < 0)
    return 2;

  if (model.opt.dumpHDF5(infile))
    solver.handleDataOutput(model.opt.hdf5,model.getProcessAdm(),
                            model.opt.saveInc,model.opt.restartInc);

  res = solver.solveProblem(infile,"Solving Advection-Diffusion problem");
  if (!res) model.printFinalNorms(solver.getTimePrm());

  return res;
}


/*!
  \brief Creates and runs the advection-diffusion problem.
*/

template<class Dim>
int runSimulator(char* infile, const AdvectionDiffusionArgs& args)
{
  if (args.timeMethod == TimeIntegration::NONE)  {
    AdvectionDiffusion integrand(Dim::dimension);
    SIMAD<Dim,AdvectionDiffusion> model(integrand,true);
    if (args.adap)
      return runSimulatorStationary<SIMSolverAdap>(infile, model, false);
    else
      return runSimulatorStationary(infile, model, true);
  }
  else if (args.timeMethod == TimeIntegration::BE ||
           args.timeMethod == TimeIntegration::BDF2 ||
           args.timeMethod == TimeIntegration::THETA) {
    AdvectionDiffusionBDF integrand(Dim::dimension,
                                    args.timeMethod,
                                    args.integrandType);
    using ADSIM = SIMAD<Dim,AdvectionDiffusionBDF>;
    ADSIM model(integrand, true);
    return runSimulatorTransientImpl(infile, model, model);
  } else if (args.timeMethod >= TimeIntegration::AM1 &&
             args.timeMethod <= TimeIntegration::AM4) {
    AdvectionDiffusionImplicit integrand(Dim::dimension,
                                         args.timeMethod,
                                         args.integrandType);
    using ADSIM = SIMAD<Dim,AdvectionDiffusionImplicit>;
    ADSIM model(integrand, true);
    TimeIntegration::SIMImplicitLMM<ADSIM> sim(model, args.timeMethod, true, "temperature");
    return runSimulatorTransientImpl(infile, sim, model);
  } else {
    AdvectionDiffusionExplicit integrand(Dim::dimension,
                                         args.timeMethod,
                                         args.integrandType);
    using ADSIM = SIMAD<Dim,AdvectionDiffusionExplicit>;
    ADSIM model(integrand, true);
    if (args.timeMethod >= TimeIntegration::AB1 &&
        args.timeMethod <= TimeIntegration::AB5) {
      TimeIntegration::SIMExplicitLMM<ADSIM> sim(model, args.timeMethod, true, "temperature");
      sim.setLinear(true);
      return runSimulatorTransientImpl(infile, sim, model);
    } else if (args.timeMethod >= TimeIntegration::HEUNEULER) {
      TimeIntegration::SIMExplicitRKE<ADSIM> sim(model, args.timeMethod, args.errTol);
      return runSimulatorTransientImpl(infile, sim, model);
    } else {
      TimeIntegration::SIMExplicitRK<ADSIM> sim(model, args.timeMethod);
      sim.setLinear(true);
      return runSimulatorTransientImpl(infile, sim, model);
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
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  char* infile = nullptr;
  AdvectionDiffusionArgs args;

  IFEM::Init(argc,argv,"Advection-Diffusion solver");
  for (int i = 1; i < argc; i++)
    if (argv[i] == infile || args.parseArg(argv[i]))
      ; // ignore the input file on the second pass
    else if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!infile) {
      infile = argv[i];
      if (!args.readXML(infile,false))
        return 1;
      i = 0;
    }
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

  if (args.adap)
    IFEM::getOptions().discretization = ASM::LRSpline;

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout) << std::endl;
  utl::profiler->stop("Initialization");

  if (args.dim == 2)
    return runSimulator<SIM2D>(infile,args);
  else
    return runSimulator<SIM3D>(infile,args);
}
