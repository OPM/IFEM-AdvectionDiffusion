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

  //! \brief Default constructor.
  //! \param[in] ad Integrand for advection-diffusion problem
  //! \param[in] alone Integrand is used stand-alone (controls time stepping)
  explicit SIMAD(Integrand& ad, bool alone = false) :
    SIMMultiPatchModelGen<Dim>(1), AD(ad),
    weakDirBC(Dim::dimension, 4.0, 1.0), inputContext("advectiondiffusion")
  {
    standalone = alone;
    Dim::myProblem = &AD;
    Dim::myHeading = "Advection-Diffusion solver";
  }

  //! \brief Constructs from given properties.
  explicit SIMAD(const SetupProps& props) :
    SIMMultiPatchModelGen<Dim>(1), AD(*props.integrand),
    weakDirBC(Dim::dimension, 4.0, 1.0), inputContext("advectiondiffusion")
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
      else if (strcasecmp(child->Value(),"subiterations") == 0) {
       utl::getAttribute(child,"max",maxSubIt);
       utl::getAttribute(child,"tol",subItTol);
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
    this->initSolution(this->getNoDOFs(),3);
    size_t n, nSols = this->getNoSolutions();
    std::string str = "temperature1";
    for (n = 0; n < nSols && n < 2; n++, str[11]++)
      this->registerField(str,solution[n]);
    return true;
  }

  //! \brief Preprocessing performed before the FEM model generation.
  void preprocessA() override
  {
    Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

    // Couple the weak Dirichlet integrand to the generic Neumann property codes
    for (const Property& p : Dim::myProps)
      if (p.pcode == Property::NEUMANN_GENERIC)
        if (Dim::myInts.find(p.pindx) == Dim::myInts.end())
          Dim::myInts.insert(std::make_pair(p.pindx,&weakDirBC));
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

    if (!tp.multiSteps())
      printFinalNorms(tp);

    return true;
  }

  //! \brief No solution postprocessing.
  bool postSolve(const TimeStep&) { return true; }

  //! \brief Evaluates and prints out solution norms.
  void printFinalNorms(const TimeStep& tp)
  {
    Vectors gNorm;
    this->setMode(SIM::RECOVERY);
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (!this->solutionNorms(tp.time,solution,gNorm))
      return;
    else if (gNorm.empty())
      return;

    IFEM::cout << ">>> Norm summary for AdvectionDiffusion <<<" << std::endl;
    IFEM::cout <<"  L2 norm |T^h| = (T^h,T^h)^0.5       : "<< gNorm[0](1);
    IFEM::cout <<"\n  H1 norm |T^h| = a(T^h,T^h)^0.5      : "<< gNorm[0](2);
    if (this->haveAnaSol() && gNorm[0].size() >= 6)
      IFEM::cout <<"\n  L2 norm |T|   = (T,T)^0.5           : "<< gNorm[0](5)
                 <<"\n  H1 norm |T|   = a(T,T)^0.5          : "<< gNorm[0](3)
                 <<"\n  L2 norm |e|   = (e,e)^0,5, e=T-T^h  : "<< gNorm[0](6)
                 <<"\n  H1 norm |e|   = a(e,e)^0.5, e=T-T^h : "<< gNorm[0](4)
                 <<"\n  Exact relative error (%)            : "
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

    exporter.registerField("u","temperature",DataExporter::SIM,results, prefix);
    exporter.setFieldValue("u", this, &this->getSolution(0));
  }

  //! \brief Set context to read from input file
  void setContext(int ctx)
  {
    std::stringstream str;
    str <<"advectiondiffusion-"<< ctx;
    inputContext = str.str();
  }

  //! \brief Returns the maximum number of iterations.
  int getMaxit() const { return maxSubIt; }
  //! \brief Returns the sub-iteration tolerance.
  double getSubItTol() const { return subItTol; }

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  //!
  //! \details Since this solver is linear, this is just a normal solve.
  SIM::ConvStatus solveIteration(TimeStep& tp)
  {
    if (Dim::msgLevel == 1)
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

protected:
  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  bool initNeumann(size_t propInd) override
  {
    typename Dim::SclFuncMap::const_iterator tit = Dim::myScalars.find(propInd);
    if (tit == Dim::myScalars.end()) return false;

    weakDirBC.setFlux(tit->second);
    AD.setFlux(tit->second);
    return true;
  }

  //! \brief Prints integrated solution norms to the log stream.
  //! \param[in] norms The norm values
  //! \param[in] w Total number of characters in the norm labels
  void printNorms (const Vectors& norms, size_t w) const override
  {
    if (norms.empty()) return;

    NormBase* norm = this->getNormIntegrand();
    const Vector& n = norms.front();

    IFEM::cout <<"Energy norm"
      << utl::adjustRight(w-8,norm->getName(1,2)) << n(2);

    if (this->haveAnaSol() && n.size() >= 4)
      IFEM::cout <<"\nExact norm"
        << utl::adjustRight(w-7,norm->getName(1,3)) << n(3)
        <<"\nExact error"
        << utl::adjustRight(w-8,norm->getName(1,4)) << n(4)
        <<"\nExact relative error (%)                : "<< 100.0*n(4)/n(3);

    size_t j = 0;
    for (const auto& prj : Dim::opt.project)
      if (++j < norms.size())
        this->printNormGroup(norms[j],n,prj.second);

    IFEM::cout << std::endl;
    delete norm;
  }

private:
  Integrand& AD; //!< Problem integrand definition
  AdvectionDiffusion::WeakDirichlet weakDirBC; //!< Weak Dirichlet integrand

  const Vector* extsol = nullptr; //!< Solution vector for adaptive simulators
  bool standalone = false; //!< If \e true, this simulator owns the VTF object
  std::string inputContext; //!< Input context
  double subItTol = 1e-4; //!< Sub-iteration tolerance
  int maxSubIt = 50; //!< Maximum number of sub-iterations
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
