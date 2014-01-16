// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionBDF.h
//!
//! \date Jun 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for time-dependent Advection-Diffusion
//!         problems.
//!
//==============================================================================

#ifndef _ADVECTIONDIFFUSION_EXPLICIT_H
#define _ADVECTIONDIFFUSION_EXPLICIT_H

#include "AdvectionDiffusion.h"
#include "ElmMats.h"
#include "Vec3.h"


/*!
  \brief Class representing the integrand of a time-dependent 
         Advection-Diffusion problem.
  \details Time stepping is done using an explicit (RK type) method
*/

class AdvectionDiffusionExplicit : public AdvectionDiffusion
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] order Temporal order (1,2)
  //! \param[in] stab Integrand formulation
  AdvectionDiffusionExplicit(unsigned short int n = 3,
                             int itg_type_=Integrand::STANDARD,
                             int form = 0);

  //! \brief Empty destructor.
  virtual ~AdvectionDiffusionExplicit();

  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents.
  //! It can also be used to implement multiple integration point loops within
  //! the same element, provided the necessary integration point values are
  //! stored internally in the object during the first integration loop.
  virtual bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t = 0);

  using AdvectionDiffusion::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const TimeDomain& time, const Vec3& X) const;

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  //! \brief Returns the integrand type
  virtual int getIntegrandType() const
  { 
    if (stab == NONE)
      return Integrand::STANDARD;

    return Integrand::SECOND_DERIVATIVES | Integrand::G_MATRIX;
  }
};

#endif
