// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionImplicit.h
//!
//! \date Oct 29 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Advection-Diffusion problems with
//!        implicit time stepping.
//!
//==============================================================================

#ifndef _ADVECTION_DIFFUSION_IMPLICIT_H
#define _ADVECTION_DIFFUSION_IMPLICIT_H

#include "AdvectionDiffusion.h"
#include "TimeIntUtils.h"


/*!
  \brief Class representing the integrand of Advection-Diffusion problem
         with implicit time stepping.
*/

class AdvectionDiffusionImplicit : public AdvectionDiffusion
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] method The time integration method to use
  //! \param[in] itg_type The integrand type to use
  //! \param[in] form Integrand formulation
  AdvectionDiffusionImplicit(unsigned short int n = 3,
                             TimeIntegration::Method method = TimeIntegration::AM1,
                             int itg_type = STANDARD, int form = 0);

  //! \brief Empty destructor.
  virtual ~AdvectionDiffusionImplicit() {}

  using AdvectionDiffusion::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt,
               const FiniteElement& fe,
               const TimeDomain& time,
               const Vec3& X) const override;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol = 0) const override;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  void setMode(SIM::SolutionMode mode) override;

protected:
  TimeIntegration::Method timeMethod; //!< Time stepping method used
};

#endif
