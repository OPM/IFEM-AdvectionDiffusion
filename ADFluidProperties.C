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
#include "tinyxml.h"
#include "Utilities.h"
#include "Vec3.h"

namespace AD {

FluidProperties::FluidProperties() :
  rho(1.0), kappa(1.0), alpha(1.0), C(1.0), Ra(1.0), Pr(1.0), scaling(PHYSICAL)
{
}


void FluidProperties::parse(const TiXmlElement* elem)
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
  const TiXmlElement* raf = elem->FirstChildElement("rayleigh");
  if (raf) {
    std::string type;
    utl::getAttribute(raf,"type",type);
    RaFdef = utl::getValue(raf,"rayleigh");
    RaF.reset(utl::parseRealFunc(RaFdef, type, false));
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
               << "\n\t\t\tThermal diffusivity, kappa = " << kappa;
  } else if (scaling == PR_RA) {
    IFEM::cout << "\n\t\t\t Prandtl number, Pr = " << Pr;
    IFEM::cout << "\n\t\t\tRayleigh number, Ra = ";
    if (RaFdef.empty())
      IFEM::cout << Ra;
    else
      IFEM::cout << RaFdef;
    IFEM::cout << std::endl;
  }
}


double FluidProperties::getMassAdvectionConstant() const
{
  if (scaling == PR_RA)
    return 1.0;

  return getMassDensity()*getHeatCapacity();
}


double FluidProperties::getDiffusionConstant(const Vec3& X) const
{
  if (scaling == PR_RA) {
    double eRa = RaF ? (*RaF)(X) : Ra;
    return 1.0/sqrt(eRa*Pr);
  }

  return getDiffusivity();
}


double FluidProperties::getReactionConstant(const Vec3& X) const
{
  if (scaling == PR_RA) {
    double eRa = RaF ? (*RaF)(X) : Ra;
    return 1.0/(sqrt(eRa*Pr)*getDiffusivity());
  }

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
