// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionExplicit.h
//!
//! \date Oct 29 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for time-dependent Advection-Diffusion.
//!
//==============================================================================

#ifndef _ADVECTION_DIFFUSION_EXPLICIT_H
#define _ADVECTION_DIFFUSION_EXPLICIT_H

#include "AdvectionDiffusion.h"


/*!
  \brief Class representing the integrand of a time-dependent
         Advection-Diffusion problem.
  \details Time stepping is done using an explicit (RK type) method.
*/

class AdvectionDiffusionExplicit : public AdvectionDiffusion
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] method The time integration method to use
  //! \param[in] itg_type The integrand type to use
  //! \param[in] form Integrand formulation
  explicit AdvectionDiffusionExplicit(unsigned short int n,
                                      TimeIntegration::Method method,
                                      int itg_type = STANDARD, int form = 0);

  using AdvectionDiffusion::finalizeElement;
  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents.
  bool finalizeElement(LocalIntegral&) override;

  using AdvectionDiffusion::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  using AdvectionDiffusion::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                  bool neumann) const override;

  //! \brief Returns the integrand type.
  int getIntegrandType() const override
  {
    return stab == NONE ? STANDARD : SECOND_DERIVATIVES | G_MATRIX;
  }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol = 0) const override;

protected:
  TimeIntegration::Method timeMethod; //!< Time stepping method used
};

#endif
