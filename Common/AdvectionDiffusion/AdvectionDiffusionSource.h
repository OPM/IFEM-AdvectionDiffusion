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


class AnaSol;
class Vec3;
namespace tinyxml2 { class XMLElement; }


namespace AD {

class FluidProperties;

/*!
  \brief Class that derives the Advection-Diffusion source function from the analytic solution.
 */

class AdvectionDiffusionAnaSolSource : public RealFunc
{
public:
  //! \brief Constructor.
  //! \param aSol Analytic solution to use
  //! \param U Advecting field
  //! \param props Fluid properties
  //! \param rField Reaction field
  //! \param stationary True for a stationary solution
  AdvectionDiffusionAnaSolSource(const AnaSol& aSol,
                                 const VecFunc& U,
                                 const FluidProperties& props,
                                 const RealFunc* rField,
                                 bool stationary);

protected:
  //! \brief Evaluates the function.
  Real evaluate(const Vec3& X) const override;

  const AnaSol& anaSol; //!< Reference to analytic solution
  const VecFunc& adVel; //!< Advecting velocity field
  const FluidProperties& props; //!< Fluid properties
  const RealFunc* reaction; //!< Reaction field
  bool stat; //!< True for a stationary solution
};


/*!
  \brief Class that derives the Advection-Diffusion source function from solution components.
 */

class AdvectionDiffusionSource : public RealFunc
{
public:
  //! \brief Constructor.
  //! \param elem XML element to parse
  //! \param props Fluid properties
  //! \param react Reaction field
  AdvectionDiffusionSource(const tinyxml2::XMLElement* elem,
                           const FluidProperties& props,
                           const RealFunc* react = nullptr);

protected:
  //! \brief Evaluates the function.
  Real evaluate(const Vec3& X) const override;

  std::unique_ptr<RealFunc> T; //!< Temperature T
  const RealFunc* reaction; //!< Reaction field
  std::unique_ptr<VecFunc> lapT; //!< Laplacian of T
  std::unique_ptr<VecFunc> gradT; //!< Gradient of T
  std::unique_ptr<VecFunc> U; //!< Advection velocity
  const FluidProperties& props; //!< Fluid properties
};

}

#endif
