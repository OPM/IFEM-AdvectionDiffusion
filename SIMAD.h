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
#include "ASMbase.h"
#include "DataExporter.h"
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

template<class Dim> class SIMAD : public Dim
{
public:
  //! \brief Default constructor.
  //! \param[in] ad Integrand for advection-diffusion problem
  //! \param[in] alone Integrand is used stand-alone (controls time stepping)
  SIMAD(AdvectionDiffusion* ad, bool alone = false) :
    Dim(1), AD(ad), weakDirBC(Dim::dimension, 4.0, 1.0)
  {
    standalone = alone;
    Dim::myProblem = AD;
  }

  //! \brief The destructor clears the VTF-file pointer, unless stand-alone.
  //! \details This is needed when the VTF-file is assumed to be owned by
  //! another SIM-object, is deleted by the SIMbase destructor of that object.
  virtual ~SIMAD()
  {
    if (!standalone)
      this->setVTF(NULL);
  }

  //! \brief Parses a data section from an input stream (depreciated).
  virtual bool parse(char*, std::istream&) { return false; }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"advectiondiffusion"))
      return this->Dim::parse(elem);

    const char* value = 0;
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())

      if (strcasecmp(child->Value(),"stabilization") == 0) {
        std::string type;
        utl::getAttribute(child,"type",type,true);
        if (type == "supg") {
          std::cout << "SUPG stabilization is enabled." << std::endl;
          AD->setStabilization(AdvectionDiffusion::SUPG);
        }
        else if (type == "gls") {
          std::cout << "GLS stabilization is enabled." << std::endl;
          AD->setStabilization(AdvectionDiffusion::GLS);
        }
        else if (type == "ms") {
          std::cout << "MS stabilization is enabled." << std::endl;
          AD->setStabilization(AdvectionDiffusion::MS);
        }
        double Cinv;
        if (utl::getAttribute(child,"Cinv",Cinv))
          AD->setCinv(Cinv);
      }
      else if ((value = utl::getValue(child,"kappa"))) {
        std::cout <<"Kappa: "<< value << std::endl;
        AD->setKappa(atof(value));
        weakDirBC.setKappa(atof(value));
      }
      else if ((value = utl::getValue(child,"prandtl"))) {
        std::cout <<"Prandtl number: "<< value << std::endl;
        AD->setPrandtlNumber(atof(value));
      }
      else if ((value = utl::getValue(child,"advectionfield"))) {
        AD->setAdvectionField(new VecFuncExpr(value,""));
        weakDirBC.setAdvectionField(new VecFuncExpr(value,""));
        std::cout <<"Advection field: "<< value << std::endl;
      }
      else if ((value = utl::getValue(child,"reactionfield"))) {
        AD->setReactionField(new EvalFunction(value));
        std::cout <<"Reaction field: "<< value << std::endl;
      }
      else if ((value = utl::getValue(child,"source"))) {
        AD->setSource(new EvalFunction(value));
        std::cout <<"Source field: "<< value << std::endl;
      }
      else if (strcasecmp(child->Value(),"anasol") == 0) {
        std::cout <<"\tAnalytical solution: Expression"<< std::endl;
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

  void init() { TimeStep dummy; this->init(dummy); }

  void init(const TimeStep& tp)
  {
    int p1, p2, p3;
    this->getPatch(0)->getOrder(p1,p2,p3);

    AD->setOrder(p1); // assumes equal ordered basis
    AD->setElements(this->getNoElms());

    // Initialize temperature solution vectors
    size_t n, nSols = this->getNoSolutions();
    temperature.resize(nSols);
    std::string str = "temperature1";
    for (n = 0; n < nSols; n++, str[11]++) {
      temperature[n].resize(this->getNoDOFs(),true);
      this->registerField(str,temperature[n]);
    }
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
  virtual bool preprocessB() { AD->setElements(this->getNoElms());return true; }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param nBlock Running VTF block counter
  virtual bool saveModel(char* fileName, int& nBlock)
  {
    return Dim::opt.format < 0 ? true : this->writeGlvG(nBlock,fileName);
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
      double normL2 = this->solutionNorms(temperature.front(),dMax,iMax,1);
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

  Vector& getSolution(int n=0) { return temperature[n]; }

  const Vector& getSolution(int n=0) const { return temperature[n]; }

  void registerFields(DataExporter& exporter)
  {
    exporter.registerField("theta","temperature",DataExporter::SIM);
    exporter.setFieldValue("theta", this, &temperature.front());
  }

  double externalEnergy(const Vectors&) const { return 0.0; }

protected:
  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::SclFuncMap::const_iterator tit = Dim::myScalars.find(propInd);
    if (tit != Dim::myScalars.end()) {
      weakDirBC.setFlux(tit->second);
      AD->setFlux(tit->second);
    }

    return tit != Dim::myScalars.end();
  }

private:
  AdvectionDiffusion* AD;
  AdvectionDiffusion::WeakDirichlet weakDirBC;

  Vectors temperature;
  bool    standalone;
};

#endif
