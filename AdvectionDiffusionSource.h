// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionSource.h
//!
//! \date Jul 16 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Advection-Diffusion source function.
//!
//==============================================================================
#ifndef AD_SOURCE_H_
#define AD_SOURCE_H_

#include "Function.h"

#include <memory>

class TiXmlElement;


namespace AD {

class FluidProperties;

/*!
  \brief Class that derives the Advection-Diffusion source function from solution components.
 */

class AdvectionDiffusionSource : public RealFunc
{
public:
  //! \brief Constructor.
  //! \param elem XML element to parse
  //! \param props Fluid properties
  AdvectionDiffusionSource(const TiXmlElement* elem,
                           const FluidProperties& props);

protected:
  //! \brief Evaluates the function.
  double evaluate(const Vec3& X) const override;

  std::unique_ptr<VecFunc> lapT; //!< Laplacian of T
  std::unique_ptr<VecFunc> gradT; //!< Gradient of T
  std::unique_ptr<VecFunc> U; //!< Advection velocity
  const FluidProperties& props; //!< Fluid properties
};

}

#endif
