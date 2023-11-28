// $Id$
//==============================================================================
//!
//! \file ADFluidProperties.C
//!
//! \date Aug 11 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Advection-Diffusion fluid properties.
//!
//==============================================================================


#include "ADFluidProperties.h"

#include "IFEM.h"
#include "Functions.h"
#include "LogStream.h"
#include "Function.h"
#include "Utilities.h"
#include "Vec3.h"

#include <cassert>
#include <map>
#include <ostream>
#include "tinyxml2.h"
#include <utility>


namespace AD {

FluidProperties::FluidProperties() :
  rho(1.0), kappa(1.0), alpha(1.0), C(1.0), Ra(1.0), Pr(1.0), scaling(PHYSICAL)
{
}


FluidProperties::~FluidProperties()
{
}


void FluidProperties::parse (const tinyxml2::XMLElement* elem)
{
  // Type of scaling to use
  std::string type;
  utl::getAttribute(elem,"type",type);
  if (type == "physical")
    scaling = PHYSICAL;
  if (type == "prandtlrayleigh")
    scaling = PR_RA;

  utl::getAttribute(elem,"rho",rho);
  utl::getAttribute(elem,"C",C);
  utl::getAttribute(elem,"alpha",alpha);
  utl::getAttribute(elem,"kappa",kappa);
  utl::getAttribute(elem,"Ra",Ra);
  utl::getAttribute(elem,"Pr",Pr);
  const tinyxml2::XMLElement* raf = elem->FirstChildElement("rayleigh");
  if (raf) {
    type.clear();
    utl::getAttribute(raf,"type",type);
    RaFdef = utl::getValue(raf,"rayleigh");
    RaF.reset(utl::parseRealFunc(RaFdef, type, false));
  }
  const tinyxml2::XMLElement* kappaf = elem->FirstChildElement("kappa");
  if (kappaf) {
    type.clear();
    utl::getAttribute(kappaf,"type",type);
    kappaFdef = utl::getValue(kappaf,"kappa");
    kappaF.reset(utl::parseRealFunc(kappaFdef,type,false));
  }
}


void FluidProperties::printLog() const
{
  static const std::map<Scaling,std::string> name_map =
        {{PHYSICAL, "Physical dimensions"},
         {PR_RA,    "Prandtl/Rayleigh numbers"}};

  IFEM::cout << "\tFluid properties:"
             << "\n\t\tScaling = " << name_map.find(scaling)->second;
  if (scaling == PHYSICAL) {
    IFEM::cout << "\n\t\t\t            Density,   rho = " << rho
               << "\n\t\t\tThermal diffusivity, kappa = ";
    if (kappaF)
      IFEM::cout << kappaFdef;
    else
      IFEM::cout << kappa;
  } else if (scaling == PR_RA) {
    IFEM::cout << "\n\t\t\t Prandtl number, Pr = " << Pr;
    IFEM::cout << "\n\t\t\tRayleigh number, Ra = ";
    if (RaFdef.empty())
      IFEM::cout << Ra;
    else
      IFEM::cout << RaFdef;
  }
  IFEM::cout << std::endl;
}


double FluidProperties::getMassAdvectionConstant() const
{
  if (scaling == PR_RA)
    return 1.0;

  return getMassDensity()*getHeatCapacity();
}


double FluidProperties::getDiffusionConstant(const Vec3& X) const
{
  if (scaling == PR_RA)
    return 1.0;

  if (kappaF)
    return (*kappaF)(X);

  return getDiffusivity();
}


double FluidProperties::getReactionConstant(const Vec3& X) const
{
  return 1.0;
}


double FluidProperties::getThermalExpansion(double T) const
{
  if (scaling == PR_RA)
    return -T;

  assert(scaling == PHYSICAL);
  return rho*alpha*T;
}

}
