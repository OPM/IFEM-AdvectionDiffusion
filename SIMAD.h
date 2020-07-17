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

#include "SIMMultiPatchModelGen.h"
#include "SIMsolution.h"
#include "SIMoutput.h"
#include "SIMconfigure.h"
#include "Property.h"
#include "ASMstruct.h"
#include "AdvectionDiffusion.h"
#include "AdvectionDiffusionSource.h"
#include "AnaSol.h"
#include "Functions.h"
#include "ExprFunctions.h"
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

template<class Dim, class Integrand=AdvectionDiffusion>
class SIMAD : public SIMMultiPatchModelGen<Dim>, public SIMsolution
{
public:
  //! \brief Setup properties.
  struct SetupProps
  {
    bool shareGrid = false; //!< True to share grid with some other simulator
    Integrand* integrand = nullptr; //!< Integrand to use
    SIMoutput* share = nullptr; //!< Simulator to share grid with
    bool standalone = false; //!< Simulator runs standalone
  };

  //! \brief Norm to use to measure subiteration convergence.
  enum SubItNorm {
    RESIDUAL_L2 = 0, //!< L2-norm of residual
    RESIDUAL_LINFTY = 1, //!< L-infty norm of residual
    RESIDUAL_L2_SCALED = 2 //!< L2-norm scaled with sqrt(n)
  };

  //! \brief Default constructor.
  //! \param[in] ad Integrand for advection-diffusion problem
  //! \param[in] alone Integrand is used stand-alone (controls time stepping)
  explicit SIMAD(Integrand& ad, bool alone = false) :
    SIMMultiPatchModelGen<Dim>(1),
    AD(ad),
    robinBC(Dim::dimension, ad),
    inputContext("advectiondiffusion")
  {
    standalone = alone;
    Dim::myProblem = &AD;
    Dim::myHeading = "Advection-Diffusion solver";
  }

  //! \brief Constructs from given properties.
  explicit SIMAD(const SetupProps& props) :
    SIMMultiPatchModelGen<Dim>(1),
    AD(*props.integrand),
    robinBC(Dim::dimension, AD),
    inputContext("advectiondiffusion")
  {
    standalone = props.standalone;
    Dim::myProblem = &AD;
    Dim::myHeading = "Advection-Diffusion solver";
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

    // To prevent the SIMbase destructor try to delete already deleted functions
    if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
    if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
  }

  //! \brief Initializes the property containers of the model.
  void clearProperties() override
  {
    // To prevent SIMbase::clearProperties deleting the analytical solution
    if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
    if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
    aCode[0] = aCode[1] = 0;

    ++this->isRefined;

    this->Dim::clearProperties();
  }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem) override
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
        AD.getFluidProperties().printLog();
      }
      else if ((value = utl::getValue(child,"advectionfield"))) {
        std::string variables;
        utl::getAttribute(child,"variables",variables);
        AD.setAdvectionField(new VecFuncExpr(value,variables));
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
        std::string type;
        utl::getAttribute(child, "type", type);
        if (type == "expression") {
          AD.setSource(new EvalFunction(value));
          IFEM::cout <<"Source field: "<< value << std::endl;
        } else if (type == "components") {
           IFEM::cout << "\tSource function:";
           AD.setSource(new AD::AdvectionDiffusionSource(child, AD.getFluidProperties()));
        }
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
      else if (strcasecmp(child->Value(),"subiterations") == 0) {
        utl::getAttribute(child,"max",maxSubIt);
        utl::getAttribute(child,"tol",subItTol);
        utl::getAttribute(child,"relax",subItRelax);
        utl::getAttribute(child,"aitken",subItAitken);
        std::string func;
        if (utl::getAttribute(child,"maxFunc",func))
          maxSubItFunc.reset(utl::parseRealFunc(func,"expression",false));
        utl::getAttribute(child,"continue_on_failure",continue_on_failure);
        std::string norm;
        utl::getAttribute(child,"norm",norm);
        std::string normType = "L2-norm of residual";
        if (norm == "L2")
          subItNorm = RESIDUAL_L2;
        else if (norm == "Linfty") {
          normType = "Linfty-norm of residual";
          subItNorm = RESIDUAL_LINFTY;
        }
        else if (norm == "L2n") {
          normType = "Scaled L2-norm of residual";
          subItNorm = RESIDUAL_L2_SCALED;
        }
        IFEM::cout << "\tUsing subiterations";
        if (func.empty())
          IFEM::cout <<"\n\t\tmax = " << maxSubIt;
        else
          IFEM::cout << "\n\t\tmax = " << func;

        IFEM::cout <<"\n\t\ttol = " << subItTol
                  <<"\n\t\tnorm = " << normType;
        if (subItRelax != 1.0) {
          IFEM::cout <<"\n\t\trelaxation = " << subItRelax;
          if (subItAitken)
            IFEM::cout << " (aitken)";
          IFEM::cout << "\n";
        }
      }
      else if (strcasecmp(child->Value(),"advection") == 0) {
        const char* value = child->FirstChild()->Value();
        if (!strcasecmp(value, "convective")) {
          IFEM::cout << "\n\tUsing convective advection operator\n";
          AD.setAdvectionForm(WeakOperators::CONVECTIVE);
        } else if (!strcasecmp(value, "conservative")) {
          IFEM::cout << "\n\tUsing conservative advection operator\n";
          AD.setAdvectionForm(WeakOperators::CONSERVATIVE);
        }
      }
      else
        this->Dim::parse(child);

    return true;
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "AdvectionDiffusion"; }

  //! \brief Initialize simulator.
  bool init(const TimeStep&)
  {
    int p1, p2, p3;
    Dim::myModel.front()->getOrder(p1,p2,p3);

    AD.setOrder(p1); // assumes equal ordered basis
    AD.setElements(this->getNoElms());

    // Initialize temperature solution vectors
    size_t n, nSols = this->getNoSolutions();
    this->initSolution(this->getNoDOFs(),nSols);
    std::string str = "temperature1";
    for (n = 0; n < nSols; n++, str[11]++)
      this->registerField(str,solution[n]);

    if (this->opt.project.find(SIMoptions::NONE) != this->opt.project.end())
      AD.setResidualNorm(true);

    return true;
  }

  //! \brief Preprocessing performed before the FEM model generation.
  //! \details This method is reimplemented to couple the weak Dirichlet
  //! integrand to the Robin property codes.
  void preprocessA() override
  {
    Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

    for (Property& p : Dim::myProps)
      if (p.pcode == Property::ROBIN) {
        if (Dim::myInts.find(p.pindx) == Dim::myInts.end())
          Dim::myInts.insert(std::make_pair(p.pindx,&robinBC));
      } else if (p.pcode == Property::DIRICHLET_ANASOL) {
        if (!Dim::mySol->getScalarSol())
          p.pcode = Property::UNDEFINED;
        else if (aCode[0] == abs(p.pindx))
          p.pcode = Property::DIRICHLET_INHOM;
        else if (aCode[0] == 0)
        {
          aCode[0] = abs(p.pindx);
          Dim::myScalars[aCode[0]] = Dim::mySol->getScalarSol();
          p.pcode = Property::DIRICHLET_INHOM;
        }
        else
          p.pcode = Property::UNDEFINED;
      } else if (p.pcode == Property::NEUMANN_ANASOL) {
        if (!Dim::mySol->getScalarSecSol())
          p.pcode = Property::UNDEFINED;
        else if (aCode[1] == p.pindx)
          p.pcode = Property::NEUMANN;
        else if (aCode[1] == 0)
        {
          aCode[1] = p.pindx;
          Dim::myVectors[aCode[1]] = Dim::mySol->getScalarSecSol();
          p.pcode = Property::NEUMANN;
        }
        else
          p.pcode = Property::UNDEFINED;
      }
  }

  //! \brief Defines the global number of elements.
  bool preprocessB() override { AD.setElements(this->getNoElms());return true; }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0) return true;

    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep&)
  {
    this->pushSolution(); // Update solution vectors between time steps
    AD.advanceStep();
    return true;
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp, bool = false)
  {
    PROFILE1("SIMAD::solveStep");

    this->setMode(tp.multiSteps() ? SIM::DYNAMIC : SIM::STATIC);
    if (Dim::msgLevel >= 0 && standalone && tp.multiSteps())
      IFEM::cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    Vector dummy;
    this->updateDirichlet(tp.time.t,&dummy);

    if (!this->assembleSystem(tp.time,solution))
      return false;

    if (!this->solveSystem(solution.front(),Dim::msgLevel-1,"temperature "))
      return false;

    if (Dim::msgLevel == 1)
    {
      size_t iMax[1];
      double dMax[1];
      double normL2 = this->solutionNorms(solution.front(),dMax,iMax,1);
      IFEM::cout <<"\n  Temperature summary:  L2-norm        : "<< normL2
                 <<"\n                       Max temperature : "<< dMax[0]
                 << std::endl;
    }

    return true;
  }

  //! \brief No solution postprocessing.
  bool postSolve(const TimeStep&) { return true; }

  //! \brief Evaluates and prints out solution norms.
  void printFinalNorms(const TimeStep& tp)
  {
    Vectors gNorm;
    this->setMode(SIM::NORMS);
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (!this->solutionNorms(tp.time,solution,gNorm))
      return;
    else if (gNorm.empty())
      return;

    this->printExactNorms(gNorm.front());
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    PROFILE1("SIMAD::saveStep");

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    TimeIntegration::Method method = AD.getTimeMethod();
    int iDump = tp.step/Dim::opt.saveInc +
                (method == TimeIntegration::NONE ? 0 : 1);
    if (!this->writeGlvS(this->getSolution(0),iDump,nBlock,
                          tp.time.t,"temperature",70))
      return false;
    else if (!standalone)
      return true;

    double param2 = method==TimeIntegration::NONE ? iDump : tp.time.t;
    return this->writeGlvStep(iDump, param2,
                              method==TimeIntegration::NONE ? 1 : 0);
  }

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(SerializeMap& data) const override
  {
    return this->saveSolution(data,this->getName());
  }

  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const SerializeMap& data) override
  {
    if (!this->restoreSolution(data,this->getName()))
      return false;

    AD.advanceStep();
    return true;
  }

  //! \brief Returns a reference to current solution vector.
  Vector& getSolution() { return solution.front(); }
  //! \brief Returns a const reference to current solution vector.
  //! \details If set, returns external vector for \a idx=0.
  const Vector& getSolution(int idx) const override
  {
    if (idx == 0 && extsol)
      return *extsol;

    return solution[idx];
  }

  //! \brief Register fields for simulation result export.
  void registerFields(DataExporter& exporter, const std::string& prefix="")
  {
    int results = DataExporter::PRIMARY;

    if (!Dim::opt.pSolOnly)
      results |= DataExporter::SECONDARY;

    if (Dim::opt.saveNorms)
      results |= DataExporter::NORMS;

    exporter.registerField("u","temperature",DataExporter::SIM,
                           results, prefix);
    exporter.setFieldValue("u", this, &this->getSolution(0),&projs,&eNorm);
  }

  //! \brief Set context to read from input file
  void setContext(int ctx)
  {
    std::stringstream str;
    str <<"advectiondiffusion-"<< ctx;
    inputContext = str.str();
  }

  //! \brief Returns the maximum number of iterations.
  int getMaxit(int iStep = 0) const
  {
    if (maxSubItFunc) {
      Vec4 X;
      X.t = iStep;
      return static_cast<int>((*maxSubItFunc)(X));
    }

    return maxSubIt;
  }

  //! \brief True to continue even if subiterations failed to converge.
  bool getSubItContinue() const { return continue_on_failure; }

  //! \brief Norm to use to measure subiteration convergence.
  SubItNorm getSubItNorm() const { return subItNorm; }

  //! \brief Returns the sub-iteration tolerance.
  double getSubItTol() const { return subItTol; }

  //! \brief Returns the sub-iteration relaxation factor.
  double getSubItRelax(double aScale)
  {
    if (subItAitken)
      subItRelax *= -aScale;
    return subItRelax;
  }

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  //!
  //! \details Since this solver is linear, this is just a normal solve.
  SIM::ConvStatus solveIteration(TimeStep& tp)
  {
    if (Dim::msgLevel == 1 && tp.multiSteps())
      IFEM::cout <<"\n  step="<< tp.step <<"  time="<< tp.time.t << std::endl;
    return this->solveStep(tp,false) ? SIM::CONVERGED : SIM::FAILURE;
  }

#ifdef HAS_PETSC
  //! \brief Sets MPI communicator for the linear equation solvers.
  void setCommunicator(const MPI_Comm* comm) { Dim::adm.setCommunicator(comm); }
#endif

  //! \brief Sets the externally provided solution vector (adaptive simulation).
  void setSol(const Vector* sol) { extsol = sol; }

  //! \brief Prints a summary of the calculated solution to std::cout.
  //! \param[in] solution The solution vector
  //! \param[in] printSol Print solution only if size is less than this value
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  void printSolutionSummary(const Vector& solution, int printSol, const char*,
                            std::streamsize outPrec = 0) override
  {
    this->SIMbase::printSolutionSummary(solution, printSol,
                                        "temperature ", outPrec);
  }

  //! \brief Set time scaling factor from ODE solver.
  //! \param scale Time scaling factor
  void setTimeScale(double scale) { AD.setTimeScale(scale); }

  //! \brief Prints integrated solution norms to the log stream.
  //! \param[in] norms The norm values
  //! \param[in] w Total number of characters in the norm labels
  void printNorms (const Vectors& norms, size_t w) const override
  {
    if (norms.empty()) return;

    NormBase* norm = this->getNormIntegrand();
    const Vector& n = norms.front();

    this->printExactNorms(norms.front());

    size_t j = 0;
    for (const auto& prj : Dim::opt.project)
      if (++j < norms.size())
        this->printNormGroup(norms[j],n,prj.second);

    IFEM::cout << std::endl;
    delete norm;
  }

  //! \brief Returns a reference to the element norm matrix.
  Matrix& getElmNorms() { return eNorm; }

  //! \brief Returns a reference to the projection vectors.
  Vectors& getProjections() { return projs; }

protected:
  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  bool initNeumann(size_t propInd) override
  {
    auto tit = Dim::myScalars.find(propInd);
    auto vit = Dim::myVectors.find(propInd);
    if (tit != Dim::myScalars.end()) {
      AD.setFlux(tit->second);
      robinBC.setFlux(tit->second);
    } else if (vit != Dim::myVectors.end())
      robinBC.setAlpha(vit->second);
    else
      return false;

    return true;
  }

  //! \brief Prints solution norms to the log stream.
  //! \param[in] gNorm The norm values
  void printExactNorms(const Vector& gNorm) const
  {
    IFEM::cout <<">>> Norm summary for AdvectionDiffusion <<<"
               <<"\n  L2 norm |T^h| = (T^h,T^h)^0.5       : "<< gNorm(1)
               <<"\n  H1 norm |T^h| = a(T^h,T^h)^0.5      : "<< gNorm(2);
    if (this->haveAnaSol() && gNorm.size() >= 6)
      IFEM::cout <<"\n  L2 norm |T|   = (T,T)^0.5           : "<< gNorm(5)
                 <<"\n  H1 norm |T|   = a(T,T)^0.5          : "<< gNorm(3)
                 <<"\n  L2 norm |e|   = (e,e)^0,5, e=T-T^h  : "<< gNorm(6)
                 <<"\n  H1 norm |e|   = a(e,e)^0.5, e=T-T^h : "<< gNorm(4)
                 <<"\n  Exact relative error (%)            : "
                 << gNorm(4)/gNorm(3)*100.0;
    IFEM::cout << std::endl;

  }

  //! \brief Prints a norm group to the log stream.
  void printNormGroup(const Vector& rNorm, const Vector& fNorm,
                      const std::string& name) const override
  {
    IFEM::cout << "\nError estimates based on >>> " << name << " <<<";
    size_t w = 36;
    if (name == "Pure residuals")
      IFEM::cout << "\n  Residual norm"
                 << utl::adjustRight(w-15,"") << rNorm[1]
                 << "\n  Effectivity index eta^res"
                 << utl::adjustRight(w-27,"")
                 << this->getEffectivityIndex(Vectors{fNorm, rNorm},1,2);
    else {
      IFEM::cout << "\n  H1 norm |T^r-T^h|"
                 << utl::adjustRight(w-19,"") << rNorm[1];
      if (rNorm.size() > 3) {
        IFEM::cout << "\n  H1 norm |T^r-T|"
                   << utl::adjustRight(w-17,"") << rNorm[2]
                   << "\n  Effectivity index eta^rec"
                   << utl::adjustRight(w-27,"")
                   << this->getEffectivityIndex(Vectors{fNorm, rNorm},1,2);
      }
    }
  }

private:
  Integrand& AD; //!< Problem integrand definition
  typename Integrand::Robin robinBC; //!< Robin integrand definition
  Vectors projs; //!< Projection vectors
  Matrix eNorm; //!< Element norms

  const Vector* extsol = nullptr; //!< Solution vector for adaptive simulators
  bool standalone = false; //!< If \e true, this simulator owns the VTF object
  std::string inputContext; //!< Input context
  double subItTol = 1e-4; //!< Sub-iteration tolerance
  int maxSubIt = 50; //!< Maximum number of sub-iterations
  SubItNorm subItNorm = RESIDUAL_L2; //!< Norm to use for measure subiteration convergence
  double subItRelax = 1.0; //!< Relaxation factor in subiterations
  std::unique_ptr<RealFunc> maxSubItFunc; //!< Maximum number of sub-iterations as a function
  bool subItAitken = false; //!< Use Aitken relaxation
  bool continue_on_failure = false; //!< Continue simulation if subiterations fails.
  int aCode[2] = {0}; //!< Analytical BC code (used by destructor)
};


//! \brief Partial specialization for configurator
template<class Dim, class Integrand>
struct SolverConfigurator< SIMAD<Dim,Integrand> > {
  //! \brief Configure a SIMAD instance.
  //! \param ad The SIMAD instance to configure
  //! \param[in] props Configuration properties
  //! \param[in] infile The input file to read
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
    ad.initSystem(ad.opt.solver);
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
