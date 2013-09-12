// $Id$
//==============================================================================
//!
//! \file main_AdvectionDiffusion.C
//!
//! \date 12 June 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Main program for the isogeometric solver for the Advection-Diffusion
//! equation.
//!
//==============================================================================

#include "IFEM.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMAD.h"
#include "SIMExplicitRK.h"
#include "SIMExplicitRKE.h"
#include "SIMSolver.h"
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


/*!
  \brief Main program for the isogeometric Advection-Diffusion equation solver.

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

  template<class Dim>
int runSimulatorStationary(bool adap, char* infile)
{
  // Create the simulation model
  SIMAD<Dim>* model = new SIMAD<Dim>(new AdvectionDiffusion(Dim::dimension), true);

  SIMinput* theSim = model;
  AdaptiveSIM* aSim = 0;
  if (adap) {
    theSim = aSim = new AdaptiveSIM(model);
    IFEM::getOptions().discretization = ASM::LRSpline;
  }

  utl::profiler->start("Model input");
  // Read in model definitions
  if (!theSim->read(infile))
    return 1;

  utl::profiler->stop("Model input");

  // Establish the FE data structures
  if (!model->preprocess())
    return 1;

  SIMoptions::ProjectionMap& pOpt = model->opt.project;
  SIMoptions::ProjectionMap::const_iterator pit;

  // Set default projection method (tensor splines only)
  if (model->opt.discretization < ASM::Spline)
    pOpt.clear(); // No projection if Lagrange/Spectral
  else if (model->opt.discretization == ASM::Spline && pOpt.empty())
    pOpt[SIMoptions::GLOBAL] = "Greville point projection";

  model->setQuadratureRule(model->opt.nGauss[0],true);
  model->setMode(SIM::STATIC);

  Matrix eNorm, ssol;
  Vector sol, load;
  Vectors projs, gNorm;

  DataExporter* exporter = NULL;
  if (model->opt.dumpHDF5(infile))
  {
    exporter = new DataExporter(true);
    exporter->registerField("u","velocity",DataExporter::SIM,
                            DataExporter::PRIMARY);
    exporter->setFieldValue("u",model, &sol);
    exporter->registerWriter(new HDF5Writer(model->opt.hdf5));
    exporter->registerWriter(new XMLWriter(model->opt.hdf5));
  }

  model->initSystem(model->opt.solver,1,1);
  model->setAssociatedRHS(0,0);
  model->init();

  if (adap) {
    aSim->initAdaptor(0,2);
    aSim->setupProjections();
    for (int iStep = 1; aSim->adaptMesh(iStep); iStep++)
    {
      if (!aSim->solveStep(infile,iStep))
        return 5;

      // New added for print norm
      model->solutionNorms(Vectors(1,sol),projs,eNorm,gNorm);
      // print norm of solution
      NormBase* norm = model->getNormIntegrand();
      std::cout << norm->getName(1,1) << ": " << gNorm[0](1) << std::endl;
      std::cout << norm->getName(1,2) << ": " << gNorm[0](2) << std::endl;
      if (gNorm[0].size() > 3) {
        std::cout << norm->getName(1,3) << ": " << gNorm[0](3) << std::endl;
        std::cout << norm->getName(1,4) << ": " << gNorm[0](4) << " :"<< std::endl;
        std::cout << norm->getName(1,5) << ": " << gNorm[0](4)/gNorm[0](3) << std::endl;
      }
      delete norm;
      //end here

      if (!aSim->writeGlv(infile,iStep,2))
        return 6;
    }
  }
  else {
    if (!model->assembleSystem())
      return 2;

    // Solve the linear system of equations
    if (!model->solveSystem(sol,1))
      return 3;

    // Project onto the splines basis
    for (pit = pOpt.begin(); pit != pOpt.end(); pit++)
      if (!model->project(ssol,sol,pit->first))
        return 4;
      else
        projs.push_back(ssol);

    model->setMode(SIM::RECOVERY);
    model->setQuadratureRule(model->opt.nGauss[1]);
    model->solutionNorms(Vectors(1,sol),projs,eNorm,gNorm);

    // print norm of solution
    NormBase* norm = model->getNormIntegrand();
    std::cout << norm->getName(1,1) << ": " << gNorm[0](1) << std::endl;
    std::cout << norm->getName(1,2) << ": " << gNorm[0](2) << std::endl;
    if (gNorm[0].size() > 3) {
      std::cout << norm->getName(1,3) << ": " << gNorm[0](3) << std::endl;
      std::cout << norm->getName(1,4) << ": " << gNorm[0](4) << std::endl;
      std::cout << norm->getName(1,5) << ": " << gNorm[0](4)/gNorm[0](3) << std::endl;
    }
    delete norm;

    const char* prefix[pOpt.size()+1];
    prefix[pOpt.size()] = 0;

    size_t j=1;
    for (pit = pOpt.begin(); pit != pOpt.end(); pit++)
      prefix[j++] = pit->second.c_str();

    // Print the norms 
    for ( j=1, pit = pOpt.begin(); pit != pOpt.end() && j<=gNorm[0].size(); pit++, j++)
    {
      std::cout <<"\n\n>>> Error estimates based on "<< pit->second<<"" <<" <<<";

      std::cout <<"\n |e|_H^1, e=u^r-u^h : " <<gNorm[j](2)<<std::endl;
      if (model->haveAnaSol() && j <= gNorm.size())
      {
        std::cout <<"\n |e|_H^1, e=u-u^r : " <<gNorm[j](3)<<std::endl;
        std::cout <<"\n |e|_H^1, e=u-u^h : " <<gNorm[0](3)<<std::endl;
        std::cout <<"\nEffectivity index (recovery) : "<< gNorm[j](2)/gNorm[0](3)<<std::endl;
      }
    }

    if (model->opt.format >= 0)
    {
      int geoBlk = 0, nBlock = 0;

      // Write VTF-file with model geometry
      if (!model->writeGlvG(geoBlk,infile))
        return 7;

      // Write Dirichlet boundary conditions
      if (!model->writeGlvBC(nBlock))
        return 8;

      // Write load vector to VTF-file
      if (!model->writeGlvV(load,"Source vector",1,nBlock))
        return 9;

      // Write solution fields to VTF-file
      if (!model->writeGlvS(sol,1,nBlock))
        return 10;

      // Write projected solution fields to VTF-file
      size_t i = 0;
      int iBlk = 100;
      for (pit = pOpt.begin(); pit != pOpt.end(); pit++, i++, iBlk += 10)
        if (!model->writeGlvP(projs[i],1,nBlock,iBlk,pit->second.c_str()))
          return 11;

      // Write element norms
      if (!model->writeGlvN(eNorm,1,nBlock,prefix))
        return 12;

      model->writeGlvStep(1,0.0,1);
    }
    model->closeGlv();
    if (exporter)
      exporter->dumpTimeLevel();
  }

  delete exporter;
  delete theSim;
  return 0;
}


  template<class Solver, class AD>
int runSimulatorTransientImpl(char* infile, TimeIntegration::Method tIt,
                              Solver& sim, AD& model)
{
  SIMSolver<Solver> solver(sim);
  // Read in model definitions
  utl::profiler->start("Model input");
  if (!model.read(infile) || !solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");

  // Establish the FE data structures
  if (!model.preprocess())
    return 1;

  model.setMode(SIM::DYNAMIC);
  model.initSystem(model.opt.solver,1,1);
  model.setQuadratureRule(model.opt.nGauss[0],true);

  DataExporter* exporter=NULL;
  if (model.opt.dumpHDF5(infile))
  {
    exporter = new DataExporter(true, model.getDumpInterval(),
                                TimeIntegration::Steps(tIt));
    exporter->registerField("theta","temperature",DataExporter::SIM,
                            DataExporter::PRIMARY);
    exporter->setFieldValue("theta", &model, &model.getSolution());
    exporter->registerWriter(new HDF5Writer(model.opt.hdf5));
    exporter->registerWriter(new XMLWriter(model.opt.hdf5));
  }

  model.init(solver.getTimePrm());

  if (!solver.solveProblem(infile, exporter))
    return 5;

  delete exporter;
  return 0;
}


  template<class Dim>
int runSimulatorTransient(char* infile, TimeIntegration::Method tIt,
                          int integrandType)
{
  if (tIt == TimeIntegration::BE || tIt == TimeIntegration::BDF2) {
    SIMAD<Dim> model(new AdvectionDiffusionBDF(Dim::dimension,
                                               tIt==TimeIntegration::BE?1:2,
                                               integrandType), true);
    return runSimulatorTransientImpl<SIMAD<Dim>, SIMAD<Dim> >(infile, tIt,
                                                              model, model);
  } else if (tIt >= TimeIntegration::HEUNEULER) {
    SIMAD<Dim> model(new AdvectionDiffusionExplicit(Dim::dimension,
                                                    integrandType), true);
    TimeIntegration::SIMExplicitRKE<SIMAD<Dim> > sim(model, tIt);
    return runSimulatorTransientImpl<TimeIntegration::SIMExplicitRKE<SIMAD<Dim> >,
                                     SIMAD<Dim> >(infile, tIt, sim, model);
  } else {
    SIMAD<Dim> model(new AdvectionDiffusionExplicit(Dim::dimension,
                                                    integrandType), true);
    TimeIntegration::SIMExplicitRK<SIMAD<Dim> > sim(model, tIt);
    return runSimulatorTransientImpl<TimeIntegration::SIMExplicitRK<SIMAD<Dim> >,
                                     SIMAD<Dim> >(infile, tIt, sim, model);
  }

  return 1; // should not be here
}


int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  SIMoptions dummy;
  char* infile = 0;
  bool adap = false;
  int i, ndim = 3;
  TimeIntegration::Method tInt = TimeIntegration::NONE;
  int integrandType = Integrand::STANDARD;

  int myPid = IFEM::Init(argc, argv);

  for (i = 1; i < argc; i++)
    if (dummy.parseOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-adap"))
      adap = true;
    else if (!strcmp(argv[i],"-2D"))
      ndim = 2;
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
	      <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
	      <<" [-lag|-spec|-LR] [-2D] [-adap] [-nGauss <n>]"
	      <<"\n       [-vtf <format> [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]\n";
    return 0;
  }

  if (adap)
    IFEM::getOptions().discretization = ASM::LRSpline;

  if (myPid == 0)
  {
    const SIMoptions& opts = IFEM::getOptions();
    std::cout <<"\n >>> IFEM Advection-Diffusion equation solver <<<"
	      <<"\n ====================================\n"
	      <<"\n Executing command:\n";
    for (i = 0; i < argc; i++) std::cout <<" "<< argv[i];
    std::cout <<"\n\nInput file: "<< infile
	      <<"\nEquation solver: "<< opts.solver
	      <<"\nNumber of Gauss points: "<< opts.nGauss[0];
    if (opts.format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (opts.format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "<< opts.nViz[0];
      if (ndim > 1) std::cout <<" "<< opts.nViz[1];
      if (ndim > 2) std::cout <<" "<< opts.nViz[2];
    }

    if (opts.discretization == ASM::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (opts.discretization == ASM::Spectral)
      std::cout <<"\nSpectral basis functions are used";
    else if (opts.discretization == ASM::LRSpline)
      std::cout <<"\nLR-spline basis functions are used";
    std::cout << std::endl;
  }
  utl::profiler->stop("Initialization");

  if (tInt == TimeIntegration::NONE) {
    if (ndim == 2)
      return runSimulatorStationary<SIM2D>(adap, infile);
    if (ndim == 3)
      return runSimulatorStationary<SIM3D>(adap, infile);
  } else {
    if (ndim == 2)
      return runSimulatorTransient<SIM2D>(infile, tInt, integrandType);
    if (ndim == 3)
      return runSimulatorTransient<SIM3D>(infile, tInt, integrandType);
  }

  return 1; // Should not be here
}
