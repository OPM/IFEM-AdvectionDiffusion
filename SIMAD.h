// $Id$
//==============================================================================
//!
//! \file SIMAD.h
//!
//! \date Jun 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for Advection-Diffusion problems.
//!
//==============================================================================

#ifndef _SIM_AD_H
#define _SIM_AD_H

#include "AdvectionDiffusion.h"
#include "AnaSol.h"
#include "ASMstruct.h"
#include "Functions.h"
#include "Property.h"
#include "SIMoutput.h"
#include "SIMSolver.h"
#include "IFEM.h"
#include "TimeStep.h"
#include "Profiler.h"
#include "Utilities.h"
#include "DataExporter.h"
#include "tinyxml.h"


/*!
  \brief Driver class for an Advection-diffusion simulator.
  \details The class incapsulates data and methods for solving a
  Advection-Diffusion problem using NURBS-based finite elements.
*/

template<class Dim, class Integrand=AdvectionDiffusion> class SIMAD : public Dim
{
public:
  //! \brief Setup properties.
  struct SetupProps
  {
    bool shareGrid = false; //!< True to share grid with some other simulator
    Integrand* integrand = nullptr; //!< Integrand to use
    SIMoutput* share = nullptr; //!< Simulator to share grid with
  };

  //! \brief Default constructor.
  //! \param[in] ad Integrand for advection-diffusion problem
  //! \param[in] alone Integrand is used stand-alone (controls time stepping)
  SIMAD(Integrand& ad, bool alone = false) :
    Dim(1), AD(ad), weakDirBC(Dim::dimension, 4.0, 1.0), inputContext("advectiondiffusion")
  {
    standalone = alone;
    Dim::myProblem = &AD;
    this->myHeading = "Advection-Diffusion solver";
  }

  //! \brief Construct from properties
  //! \param props The properties
  SIMAD(SetupProps& props) :
    Dim(1), AD(*props.integrand),
    weakDirBC(Dim::dimension, 4.0, 1.0), inputContext("advectiondiffusion")
  {
    standalone = false;
    Dim::myProblem = &AD;
    this->myHeading = "Advection-Diffusion solver";
  }

  //! \brief The destructor clears the VTF-file pointer, unless stand-alone.
  //! \details This is needed when the VTF-file is assumed to be owned by
  //! another SIM-object and is deleted by the SIMbase destructor of that one.
  virtual ~SIMAD()
  {
    if (!standalone)
      this->setVTF(nullptr);
    Dim::myProblem = nullptr;
    Dim::myInts.clear();
  }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),inputContext.c_str()))
      return this->Dim::parse(elem);

    const char* value = 0;
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())

      if (strcasecmp(child->Value(),"stabilization") == 0) {
        std::string type;
        utl::getAttribute(child,"type",type,true);
        if (type == "supg") {
          IFEM::cout <<"SUPG stabilization is enabled."<< std::endl;
          AD.setStabilization(AdvectionDiffusion::SUPG);
        }
        else if (type == "gls") {
          IFEM::cout <<"GLS stabilization is enabled."<< std::endl;
          AD.setStabilization(AdvectionDiffusion::GLS);
        }
        else if (type == "ms") {
          IFEM::cout <<"MS stabilization is enabled."<< std::endl;
          AD.setStabilization(AdvectionDiffusion::MS);
        }
        double Cinv;
        if (utl::getAttribute(child,"Cinv",Cinv))
          AD.setCinv(Cinv);
      }
      else if (!strcasecmp(child->Value(),"fluidproperties")) {
        AD.getFluidProperties().parse(child);
        weakDirBC.getFluidProperties().parse(child);
        AD.getFluidProperties().printLog();
      }
      else if ((value = utl::getValue(child,"advectionfield"))) {
        std::string variables;
        utl::getAttribute(child,"variables",variables);
        AD.setAdvectionField(new VecFuncExpr(value,variables));
        weakDirBC.setAdvectionField(new VecFuncExpr(value,variables));
        IFEM::cout <<"Advection field: "<< value;
        if (!variables.empty())
          IFEM::cout <<" (variables: "<< variables <<")";
        IFEM::cout << std::endl;
      }
      else if ((value = utl::getValue(child,"reactionfield"))) {
        AD.setReactionField(new EvalFunction(value));
        IFEM::cout <<"Reaction field: "<< value << std::endl;
      }
      else if ((value = utl::getValue(child,"source"))) {
        AD.setSource(new EvalFunction(value));
        IFEM::cout <<"Source field: "<< value << std::endl;
      }
      else if (strcasecmp(child->Value(),"anasol") == 0) {
        IFEM::cout <<"\tAnalytical solution: Expression"<< std::endl;
        if (!Dim::mySol)
          Dim::mySol = new AnaSol(child);

        // Define the analytical boundary traction field
        int code = 0;
        if (utl::getAttribute(child,"code",code)) {
          if (code > 0 && Dim::mySol->getScalarSecSol())
          {
            this->setPropertyType(code,Property::NEUMANN);
            Dim::myVectors[code] = Dim::mySol->getScalarSecSol();
          }
        }
      }
      else
        this->Dim::parse(child);

    return true;
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "AdvectionDiffusion"; }

  //! \brief Initialize simulator.
  bool init(const TimeStep& tp=TimeStep())
  {
    int p1, p2, p3;
    Dim::myModel.front()->getOrder(p1,p2,p3);

    AD.setOrder(p1); // assumes equal ordered basis
    AD.setElements(this->getNoElms());

    // Initialize temperature solution vectors
    size_t n, nSols = this->getNoSolutions();
    temperature.resize(3);
    temperature[0].resize(this->getNoDOFs(),true);
    std::string str = "temperature1";
    for (n = 0; n < nSols; n++, str[11]++) {
      temperature[n].resize(this->getNoDOFs(),true);
      this->registerField(str,temperature[n]);
    }
    return true;
  }

  //! \brief Preprocessing performed before the FEM model generation.
  virtual void preprocessA()
  {
    Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

    // Couple the weak Dirichlet integrand to the generic Neumann property codes
    PropertyVec::iterator p;
    for (p = Dim::myProps.begin(); p != Dim::myProps.end(); p++)
      if (p->pcode == Property::NEUMANN_GENERIC)
        if (Dim::myInts.find(p->pindx) == Dim::myInts.end())
          Dim::myInts.insert(std::make_pair(p->pindx,&weakDirBC));
  }

  //! \brief Defines the global number of elements.
  virtual bool preprocessB() { AD.setElements(this->getNoElms());return true; }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  virtual bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0) return true;

    nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp)
  {
    // Update temperature vectors between time steps
    const int nNusols = temperature.size();
    for (int n = nNusols-1; n > 0; n--)
      temperature[n] = temperature[n-1];

    AD.advanceStep();

    return true;
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    PROFILE1("SIMAD::solveStep");

    if (Dim::msgLevel >= 0 && standalone)
      IFEM::cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    Vector dummy;
    this->updateDirichlet(tp.time.t, &dummy);

    if (!this->assembleSystem(tp.time, temperature))
      return false;

    if (!this->solveSystem(temperature.front(), Dim::msgLevel-1,"temperature "))
      return false;

    if (Dim::msgLevel == 1)
    {
      size_t iMax[1];
      double dMax[1];
      double normL2 = this->solutionNorms(temperature.front(),dMax,iMax,1);
      IFEM::cout <<"\n  Temperature summary:  L2-norm        : "<< normL2
                 <<"\n                       Max temperature : "<< dMax[0]
                 << std::endl;
    }

    return true;
  }

  //! \brief No solution postprocessing.
  bool postSolve(const TimeStep& tp,bool) { return true; }

  //! \brief Evaluates and prints out solution norms.
  void printFinalNorms(const TimeStep& tp)
  {
    Vectors gNorm;
    this->setMode(SIM::RECOVERY);
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (!this->solutionNorms(tp.time,temperature,gNorm))
      return;
    else if (gNorm.empty())
      return;

    IFEM::cout <<"L2 norm |t^h| = a(t^h,t^h)^0.5      : "<< gNorm[0](1);
    IFEM::cout <<"\nH1 norm |t^h| = a(t^h,t^h)^0.5      : "<< gNorm[0](2);
    if (this->haveAnaSol() && gNorm[0].size() >= 6)
      IFEM::cout <<"\nL2 norm |t|   = (t,t)^0.5           : "<< gNorm[0](3)
                 <<"\nH1 norm |t|   = a(t,t)^0.5          : "<< gNorm[0](5)
                 <<"\nL2 norm |e|   = (e,e)^0,5, e=t-t^h  : "<< gNorm[0](4)
                 <<"\nH1 norm |e|   = a(e,e)^0.5, e=t-t^h : "<< gNorm[0](6)
                 <<"\nExact relative error (%)            : "
                 << gNorm[0](6)/gNorm[0](5)*100.0;
    IFEM::cout << std::endl;
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    PROFILE1("SIMAD::saveStep");

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;
    if (!this->writeGlvS1(this->getSolution(0),iDump,nBlock,
                          tp.time.t,"temperature",89))
      return false;
    else if (!standalone)
      return true;

    return this->writeGlvStep(iDump,tp.time.t);
  }

  //! \brief Obtain a reference to a solution vector.
  Vector& getSolution(int n=0) { return temperature[n]; }

  //! \brief Obtain a reference to a solution vector.
  //! \details If set, returns external vector for index 0.
  const Vector& getSolution(int n=0) const
  {
    if (n == 0 && solution)
      return *solution;

    return temperature[n];
  }

  //! \brief Register fields for simulation result export.
  void registerFields(DataExporter& exporter, const std::string& prefix="")
  {
    exporter.registerField("theta","temperature",DataExporter::SIM,
                           DataExporter::PRIMARY|DataExporter::RESTART,
                           prefix);
    exporter.setFieldValue("theta", this, &this->getSolution(0));
  }

  double externalEnergy(const Vectors&) const { return 0.0; }

  //! \brief Set context to read from input file
  void setContext(int ctx)
  {
    std::stringstream str;
    str <<"advectiondiffusion-"<< ctx;
    inputContext = str.str();
  }

#ifdef HAS_PETSC
  //! \brief Set MPI communicator for the linear equation solvers
  //! \param comm The communicator to use
  void setCommunicator(const MPI_Comm* comm)
  {
    this->adm.setCommunicator(comm);
  }
#endif

  //! \brief Set externally provided solution vector (adaptive).
  void setSol(const Vector* sol) { solution = sol; }

protected:
  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::SclFuncMap::const_iterator tit = Dim::myScalars.find(propInd);
    if (tit == Dim::myScalars.end()) return false;

    weakDirBC.setFlux(tit->second);
    AD.setFlux(tit->second);
    return true;
  }

private:
  Integrand& AD; //!< Problem integrand definition
  AdvectionDiffusion::WeakDirichlet weakDirBC; //!< Weak Dirichlet integrand

  const Vector* solution = nullptr; //!< Externally provided solution vector (adaptive simulators).
  Vectors temperature; //!< Temperature solutioni vectors
  bool    standalone; //!< True if simulator runs standalone (i.e. we own the VTF object).
  std::string inputContext; //!< Input context
};


//! \brief Partial specialization for configurator
template<class Dim, class Integrand>
struct SolverConfigurator< SIMAD<Dim,Integrand> > {
  int setup(SIMAD<Dim,Integrand>& ad,
            const typename SIMAD<Dim,Integrand>::SetupProps& props, char* infile)
  {
    utl::profiler->start("Model input");

    if (props.shareGrid)
      // Let the turbulence solver use the same grid as the velocity solver
      ad.clonePatches(props.share->getFEModel(), props.share->getGlob2LocMap());

    // Reset the global element and node numbers
    ASMstruct::resetNumbering();
    if (!ad.read(infile))
      return 2;

    utl::profiler->stop("Model input");

    // Preprocess the model and establish data structures for the algebraic system
    if (!ad.preprocess())
      return 3;

    // Initialize the linear solvers
    ad.setMode(SIM::DYNAMIC);
    ad.initSystem(ad.opt.solver,1,1,false);
    ad.setQuadratureRule(ad.opt.nGauss[0]);

    // Time-step loop
    ad.init(TimeStep());
    if (props.share)
      ad.setVTF(props.share->getVTF());
    ad.setInitialConditions();

    return 0;
  }
};

#endif
