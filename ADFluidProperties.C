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
  if (type == "pr/ra")
    scaling = PR_RA;

  utl::getAttribute(elem,"rho",rho);
  utl::getAttribute(elem,"C",C);
  utl::getAttribute(elem,"alpha",alpha);
  utl::getAttribute(elem,"kappa",kappa);
  utl::getAttribute(elem,"Ra",Ra);
  utl::getAttribute(elem,"Pr",Pr);
}


void FluidProperties::printLog() const
{
  static const std::map<Scaling,std::string> name_map =
        {{PHYSICAL, "Physical dimensions"},
         {PR_RA,    "Prandtl/Rayleigh"}};

  IFEM::cout << "\tFluid properties:"
             << "\n\t\tScaling = " << name_map.find(scaling)->second;
  if (scaling == PHYSICAL) {
    IFEM::cout << "\n\t\t\t            Density,   rho = " << rho
               << "\n\t\t\tThermal diffusivity, kappa = " << kappa;
  } else if (scaling == PR_RA) {
    IFEM::cout << "\n\t\t\tRayleigh number, Ra = " << Ra;
    IFEM::cout << "\n\t\t\t Prandtl number, Pr = " << Pr;
  }
  IFEM::cout << std::endl;
}


double FluidProperties::getMassAdvectionConstant() const
{
  if (scaling == PR_RA)
    return 1.0;

  return getMassDensity()*getHeatCapacity();
}


double FluidProperties::getDiffusionConstant() const
{
  if (scaling == PR_RA)
    return 1.0/sqrt(Ra*Pr);

  return getDiffusivity();
}


double FluidProperties::getReactionConstant() const
{
  if (scaling == PR_RA)
    return 1.0/(sqrt(Ra*Pr)*getDiffusivity());

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
