// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionSource.C
//!
//! \date Jul 17 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Advection-Diffusion source function.
//!
//==============================================================================

#include "AdvectionDiffusionSource.h"
#include "ADFluidProperties.h"

#include "AnaSol.h"
#include "Functions.h"
#include "IFEM.h"
#include "LogStream.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"

#include <ostream>
#include <string>
#include <tinyxml2.h>


namespace AD {

AdvectionDiffusionAnaSolSource::
AdvectionDiffusionAnaSolSource (const AnaSol& aSol,
                                const VecFunc& U,
                                const FluidProperties& prop,
                                const RealFunc* rField,
                                bool stationary) :
  anaSol(aSol), adVel(U), props(prop), reaction(rField), stat(stationary)
{
}


Real AdvectionDiffusionAnaSolSource::evaluate (const Vec3& X) const
{
  const RealFunc& T = *anaSol.getScalarSol();
  const Vec3 dT = T.gradient(X);
  const SymmTensor d2T = T.hessian(X);

  Real result = -props.getDiffusionConstant(X) * d2T.trace();
  if (props.kappaFunc())
    result -= dT*props.kappaFunc()->gradient(X);

  const Vec3 u = adVel(X);
  result += props.getMassAdvectionConstant()*(u*dT);

  if (reaction)
    result += props.getReactionConstant(X)*(*reaction)(X) * T(X);

  if (!stat)
    result += T.timeDerivative(X);

  return result;
}


AdvectionDiffusionSource::AdvectionDiffusionSource (const tinyxml2::XMLElement* elem,
                                                    const FluidProperties& prop,
                                                    const RealFunc* react) :
  reaction(react), props(prop)
{
  const tinyxml2::XMLElement* tt = elem->FirstChildElement("temperature");
  if (tt) {
    std::string type;
    utl::getAttribute(tt, "type", type);
    std::string func = utl::getValue(tt, "temperature");
    IFEM::cout << "\n\t\t         Temperature (" << type << ")";
    T.reset(utl::parseRealFunc(func, type));
  }

  const tinyxml2::XMLElement* tg = elem->FirstChildElement("temperature_grad");
  if (tg) {
    std::string type;
    utl::getAttribute(tg, "type", type);
    std::string func = utl::getValue(tg, "temperature_grad");
    IFEM::cout << "\n\t\tTemperature gradient (" << type << ")";
    gradT.reset(utl::parseVecFunc(func, type));
  }

  const tinyxml2::XMLElement* lap = elem->FirstChildElement("laplacian");
  if (lap) {
    std::string type;
    utl::getAttribute(lap, "type", type);
    std::string func = utl::getValue(lap, "laplacian");
    IFEM::cout << "\n\t\t           Laplacian (" << type << ")";
    lapT.reset(utl::parseVecFunc(func, type));
  }

  const tinyxml2::XMLElement* grad = elem->FirstChildElement("velocity");
  if (grad) {
    std::string type;
    utl::getAttribute(grad, "type", type);
    std::string func = utl::getValue(grad, "velocity");
    IFEM::cout << "\n\t\t            Velocity (" << type <<  ")";
    U.reset(utl::parseVecFunc(func, type));
  }

  IFEM::cout << std::endl;

  ncmp = gradT->dim();
}


Real AdvectionDiffusionSource::evaluate (const Vec3& X) const
{
  Vec3 lap = (*lapT)(X);
  Vec3 grad = (*gradT)(X);
  Vec3 u = (*U)(X);

  Real result = -props.getDiffusionConstant(X)*lap.sum() +
                 props.getMassAdvectionConstant()*(u*grad);
  if (T && reaction)
    result += props.getReactionConstant(X)*(*reaction)(X)*(*T)(X);

  return result;
}

}
