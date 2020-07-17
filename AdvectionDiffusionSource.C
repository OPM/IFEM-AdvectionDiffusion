// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionSource.C
//!
//! \date Jul 17 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Stokes source function.
//!
//==============================================================================

#include "AdvectionDiffusionSource.h"
#include "ExprFunctions.h"
#include "IFEM.h"
#include "ADFluidProperties.h"
#include "Functions.h"
#include "Vec3Oper.h"
#include "Utilities.h"

#include <tinyxml.h>


namespace AD {

AdvectionDiffusionSource::AdvectionDiffusionSource(const TiXmlElement* elem,
                                                   const FluidProperties& prop) :
  props(prop)
{
  const TiXmlElement* tg = elem->FirstChildElement("temperature_grad");
  if (tg) {
    std::string type;
    utl::getAttribute(tg, "type", type);
    std::string func = utl::getValue(tg, "temperature_grad");
    IFEM::cout << "\n\t\tTemperature gradient (" << type << ") = " << func;
    gradT.reset(utl::parseVecFunc(func, type));
  }

  const TiXmlElement* lap = elem->FirstChildElement("laplacian");
  if (lap) {
    std::string type;
    utl::getAttribute(lap, "type", type);
    std::string func = utl::getValue(lap, "laplacian");
    IFEM::cout << "\n\t\tLaplacian (" << type << ") = " << func;
    lapT.reset(utl::parseVecFunc(func, type));
  }

  const TiXmlElement* grad = elem->FirstChildElement("velocity");
  if (grad) {
    std::string type;
    utl::getAttribute(grad, "type", type);
    std::string func = utl::getValue(grad, "velocity");
    IFEM::cout << "\n\t\tVelocity (" << type <<  ") = " << func;
    U.reset(utl::parseVecFunc(func, type));
  }

  IFEM::cout << std::endl;

  ncmp = gradT->dim();
}


double AdvectionDiffusionSource::evaluate(const Vec3& X) const
{
  Vec3 lap = (*lapT)(X);
  Vec3 grad = (*gradT)(X);
  Vec3 u = (*U)(X);

  return -props.getDiffusivity() * lap.sum() + props.getMassAdvectionConstant()*(u*grad);
}

}
