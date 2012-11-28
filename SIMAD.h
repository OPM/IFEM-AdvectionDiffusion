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

#include "AlgEqSystem.h"
#include "AnaSol.h"
#include "Functions.h"
#include "Profiler.h"
#include "Property.h"
#include "TimeStep.h"
#include "Utilities.h"

#include "tinyxml.h"

/*!
  \brief Driver class for an Advection-diffusion simulator.
  \details The class incapsulates data and methods for solving a
  Advection-Diffusion problem using NURBS-based finite elements.
*/


  template<class Dim>
class SIMAD : public Dim
{
public:
  //! \brief Default constructor.
  //! \param[in] ad Integrand for advection-diffusion problem
  //! \param[in] standalone Integrand is used standalone (controls time stepping)
  SIMAD(AdvectionDiffusion* ad, bool standalone_=false) : 
    Dim(1), AD(ad), standalone(standalone_)
  {
    Dim::myProblem = AD;
  }

  //! \brief The destructor frees the dynamically allocated data.
  virtual ~SIMAD()
  {
    if (!standalone)
      Dim::setVTF(0);
  }

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is)
  {
    return false;
  }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"advectiondiffusion"))
      return Dim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    while (child) {
      if (strcasecmp(child->Value(),"stabilization") == 0) {
        std::string type;
        utl::getAttribute(child,"type",type,true);
        if (type == "supg") {
          std::cout << "SUPG stabilization is enabled." << std::endl;
          AD->setStabilization(AdvectionDiffusion::SUPG);
        }
        if (type == "gls") {
          std::cout << "GLS stabilization is enabled." << std::endl;
          AD->setStabilization(AdvectionDiffusion::GLS);
        }
        if (type == "ms") {
          std::cout << "MS stabilization is enabled." << std::endl;
          AD->setStabilization(AdvectionDiffusion::MS);
        }
        double Cinv;
        if (utl::getAttribute(child,"Cinv",Cinv))
          AD->setCinv(Cinv);
      } else if (strcasecmp(child->Value(),"kappa") == 0) {
        double kappa = atof(child->FirstChild()->Value());
        std::cout << "Kappa: " << kappa << std::endl;
        AD->setKappa(kappa);
      } else if (strcasecmp(child->Value(),"prandtl") == 0) {
        double Pr = atof(child->FirstChild()->Value());
        std::cout << "Prandtl number: " << Pr << std::endl;
        AD->setPrandtlNumber(Pr);
      } else if (strcasecmp(child->Value(),"advectionfield") == 0) {
        VecFunc* func = new VecFuncExpr(utl::getValue(child,"advectionfield"),"");
        AD->setAdvectionField(func);
        std::cout << "Advection field: " << utl::getValue(child,"advectionfield") << std::endl;
      } else if (strcasecmp(child->Value(),"reactionfield") == 0) {
        RealFunc* func = new EvalFunction(utl::getValue(child,"reactionfield"));
        AD->setReactionField(func);
        std::cout << "Reaction field: " << utl::getValue(child,"reactionfield") << std::endl;
      } else if (strcasecmp(child->Value(),"source") == 0) {
        RealFunc* src = new EvalFunction(utl::getValue(child,"source"));
        AD->setSource(src);
        std::cout << "Source field: " << utl::getValue(child,"source") << std::endl;
      } else if (strcasecmp(child->Value(),"anasol") == 0) {
        int code = 0;
        std::string type;
        utl::getAttribute(child,"type",type,true);
        std::cout <<"\tAnalytical solution: Expression"<< std::endl;
        if (!Dim::mySol)
          Dim::mySol = new AnaSol(child);

        // Define the analytical boundary traction field
        if (utl::getAttribute(child,"code",code)) {
          if (code > 0 && Dim::mySol && Dim::mySol->getScalarSecSol())
          {
            Dim::setPropertyType(code,Property::NEUMANN);
            Dim::myVectors[code] = Dim::mySol->getScalarSecSol();
          }
        }
      } else
        Dim::parse(child);

      child = child->NextSiblingElement();
    }

    return true;
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "AdvectionDiffusion"; }

  void init()
  {
    TimeStep dummy;
    init(dummy);
  }

  void init(const TimeStep& tp)
  {
    AD->setElements(Dim::getNoElms());

    // Initialize temperature solution vectors
    size_t n, nSols = Dim::getNoSolutions();
    temperature.resize(nSols);
    std::string str = "temperature1";
    for (n = 0; n < nSols; n++, str[11]++) {
      temperature[n].resize(Dim::getNoDOFs(),true);
      Dim::registerField(str,temperature[n]);
    }
  }

  virtual bool preprocess(const std::vector<int>& ignored = std::vector<int>(),
			  bool fixDup = false)
  {
    bool result = Dim::preprocess(ignored, fixDup);
    AD->setElements(Dim::getNoElms());
    return result;
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  virtual bool saveModel(char* fileName, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;
    return Dim::writeGlvG(nBlock, fileName);
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp)
  {
    // Update temperature vectors between time steps
    const int nNusols = temperature.size();
    for (int n = nNusols-1; n > 0; n--)
      temperature[n] = temperature[n-1];

    return true;
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    PROFILE1("SIMAD::solveStep");

    if (Dim::myPid == 0 && Dim::msgLevel >= 0 && standalone)
      std::cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

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
      double normL2 = this->solutionNorms(temperature.front(),dMax,iMax);
      if (Dim::myPid == 0)
        std::cout <<"Temperature summary: L2-norm        : "<< normL2
                 <<"\n                   Max temperature : "<< dMax[0]
                  << std::endl;
    }

    return true;
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

    // Write solution fields
    bool result = this->writeGlvS(temperature.front(), iDump, nBlock,
                                  tp.time.t, true, "temperature", 89);

    return (!standalone || this->writeGlvStep(iDump, tp.time.t)) && result;
  }

  Vector& getSolution(int n=0)
  {
    return temperature[n];
  }

  const Vector& getSolution(int n=0) const
  {
    return temperature[n];
  }

protected:
  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::SclFuncMap::const_iterator tit = Dim::myScalars.find(propInd);
    if (tit != Dim::myScalars.end())
      AD->setFlux(tit->second);

    return tit != Dim::myScalars.end();
  }
private:
  AdvectionDiffusion* AD;
  Vectors temperature;
  bool standalone;
};

#endif
