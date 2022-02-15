// $Id$
//==============================================================================
//!
//! \file SIMAD.C
//!
//! \date Jun 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for Advection-Diffusion problems.
//!
//==============================================================================

#include "SIMAD.h"

#include "AdvectionDiffusion.h"
#include "AdvectionDiffusionBDF.h"
#include "AdvectionDiffusionExplicit.h"
#include "AdvectionDiffusionImplicit.h"
#include "AdvectionDiffusionSource.h"

#include "AnaSol.h"
#include "ASMstruct.h"
#include "DataExporter.h"
#include "EqualOrderOperators.h"
#include "ExprFunctions.h"
#include "IFEM.h"
#include "IntegrandBase.h"
#include "LogStream.h"
#include "Profiler.h"
#include "Property.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMoptions.h"
#include "TimeDomain.h"
#include "TimeIntUtils.h"
#include "TimeStep.h"
#include "Utilities.h"

#include <cmath>
#include <cstdlib>
#include <map>
#include <ostream>
#include "strings.h"
#include <vector>
#include "tinyxml.h"
#include <utility>


template<class Dim, class Integrand>
SIMAD<Dim,Integrand>::SIMAD (Integrand& ad, bool alone) :
  SIMMultiPatchModelGen<Dim>(1),
  AD(ad),
  robinBC(Dim::dimension, ad),
  inputContext("advectiondiffusion")
{
  standalone = alone;
  Dim::myProblem = &AD;
  Dim::myHeading = "Advection-Diffusion solver";
}


template<class Dim, class Integrand>
SIMAD<Dim,Integrand>::~SIMAD ()
{
  if (!standalone)
    this->setVTF(nullptr);
  Dim::myProblem = nullptr;
  Dim::myInts.clear();

  // To prevent the SIMbase destructor try to delete already deleted functions
  if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
  if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
}


template<class Dim, class Integrand>
void SIMAD<Dim,Integrand>::clearProperties ()
{
  // To prevent SIMbase::clearProperties deleting the analytical solution
  if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
  if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
  aCode[0] = aCode[1] = 0;

  ++this->isRefined;

  this->Dim::clearProperties();
}


template<class Dim, class Integrand>
bool SIMAD<Dim,Integrand>::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),inputContext.c_str()))
    return this->Dim::parse(elem);

  const char* value;
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
    else if ((value = utl::getValue(child, "advection"))) {
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


template<class Dim, class Integrand>
bool SIMAD<Dim,Integrand>::init (const TimeStep&)
{
  int p1 = 1, p2, p3;
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


template<class Dim, class Integrand>
void SIMAD<Dim,Integrand>::preprocessA ()
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


template<class Dim, class Integrand>
bool SIMAD<Dim,Integrand>::preprocessB ()
{
  AD.setElements(this->getNoElms());
  return true;
}


template<class Dim, class Integrand>
bool SIMAD<Dim,Integrand>::saveModel (char* fileName, int& geoBlk, int&)
{
  if (Dim::opt.format < 0) return true;

  return this->writeGlvG(geoBlk,fileName);
}


template<class Dim, class Integrand>
bool SIMAD<Dim,Integrand>::advanceStep (TimeStep&)
{
  this->pushSolution(); // Update solution vectors between time steps
  AD.advanceStep();
  return true;

}


template<class Dim, class Integrand>
bool SIMAD<Dim,Integrand>::solveStep (const TimeStep& tp, bool)
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

  if (Dim::msgLevel >= 1)
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


template<class Dim, class Integrand>
void SIMAD<Dim,Integrand>::printFinalNorms (const TimeStep& tp)
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


template<class Dim, class Integrand>
bool SIMAD<Dim,Integrand>::saveStep (const TimeStep& tp, int& nBlock)
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


template<class Dim, class Integrand>
bool SIMAD<Dim,Integrand>::serialize (SerializeMap& data) const
{
  return this->saveSolution(data,this->getName());
}


template<class Dim, class Integrand>
bool SIMAD<Dim,Integrand>::deSerialize (const SerializeMap& data)
{
  if (!this->restoreSolution(data,this->getName()))
    return false;

  AD.advanceStep();
  return true;
}


template<class Dim, class Integrand>
void SIMAD<Dim,Integrand>::registerFields (DataExporter& exporter,
                                           const std::string& prefix)
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


template<class Dim, class Integrand>
void SIMAD<Dim,Integrand>::setContext (int ctx)
{
  std::stringstream str;
  str <<"advectiondiffusion-"<< ctx;
  inputContext = str.str();
}


template<class Dim, class Integrand>
SIM::ConvStatus SIMAD<Dim,Integrand>::solveIteration (TimeStep& tp)
{
  if (Dim::msgLevel == 1 && tp.multiSteps())
    IFEM::cout <<"\n  step="<< tp.step <<"  time="<< tp.time.t << std::endl;
  return this->solveStep(tp,false) ? SIM::CONVERGED : SIM::FAILURE;
}


template<class Dim, class Integrand>
void SIMAD<Dim,Integrand>::printNorms (const Vectors& norms, size_t) const
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


template<class Dim, class Integrand>
bool SIMAD<Dim,Integrand>::initNeumann (size_t propInd)
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


template<class Dim, class Integrand>
void SIMAD<Dim,Integrand>::printExactNorms (const Vector& gNorm) const
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


template<class Dim, class Integrand>
void SIMAD<Dim,Integrand>::printNormGroup (const Vector& rNorm,
                                           const Vector& fNorm,
                                           const std::string& name) const
{
  IFEM::cout << "\nError estimates based on >>> " << name << " <<<";
  size_t w = 36;
  if (name == "Pure residuals") {
    IFEM::cout << "\n  Residual norm"
                 << utl::adjustRight(w-15,"") << rNorm[1];
    if (fNorm.size() > 2)
      IFEM::cout << "\n  Effectivity index eta^res"
                 << utl::adjustRight(w-27,"")
                 << this->getEffectivityIndex(Vectors{fNorm, rNorm},1,2);
  } else {
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


template<class Dim, class Integrand>
int SolverConfigurator<SIMAD<Dim,Integrand>>::
setup (SIMAD<Dim,Integrand>& ad,
       const typename SIMAD<Dim,Integrand>::SetupProps& props,
       char* infile)
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

//! \brief Helper macro for explicit template instancing.
#define INSTANCE(T) \
  template class SIMAD<SIM1D,T>; \
  template class SIMAD<SIM2D,T>; \
  template class SIMAD<SIM3D,T>; \
  template struct SolverConfigurator<SIMAD<SIM1D,T>>; \
  template struct SolverConfigurator<SIMAD<SIM2D,T>>; \
  template struct SolverConfigurator<SIMAD<SIM3D,T>>;

INSTANCE(AdvectionDiffusion)
INSTANCE(AdvectionDiffusionBDF)
INSTANCE(AdvectionDiffusionExplicit)
INSTANCE(AdvectionDiffusionImplicit)
