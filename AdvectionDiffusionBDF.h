// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionBDF.h
//!
//! \date Jun 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for time-dependent Advection-Diffusion.
//!
//==============================================================================

#ifndef _ADVECTION_DIFFUSION_BDF_H
#define _ADVECTION_DIFFUSION_BDF_H

#include "AdvectionDiffusion.h"
#include "BDF.h"
#include "TimeIntUtils.h"
#include "Fields.h"
#include "TimeIntUtils.h"
#include <array>
#include <memory>


/*!
  \brief Class representing the integrand of a time-dependent
         Advection-Diffusion problem.

  \details Time stepping is done using BDF1/BDF2
*/

class AdvectionDiffusionBDF : public AdvectionDiffusion
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] method The time integration method to use
  //! \param[in] itg_type The integrand type to use
  //! \param[in] useALE If \e true, use ALE formulation
  AdvectionDiffusionBDF(unsigned short int n = 3,
                        TimeIntegration::Method method = TimeIntegration::BE,
                        int itg_type = STANDARD, bool useALE = false);

  //! \brief Empty destructor.
  virtual ~AdvectionDiffusionBDF() {}

  using AdvectionDiffusion::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param A Local integral for element
  bool initElement(const std::vector<int>& MNPC, LocalIntegral& A) override;

  using AdvectionDiffusion::finalizeElement;
  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents.
  //! It is used here to add stabilization terms to the element matrices.
  bool finalizeElement(LocalIntegral&) override;

  using AdvectionDiffusion::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;

  //! \brief Returns the integrand type.
  int getIntegrandType() const override
  {
    return stab == NONE ? STANDARD : SECOND_DERIVATIVES | G_MATRIX;
  }

  //! \brief Advances the time stepping scheme.
  void advanceStep() override { bdf.advanceStep(); }

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
  bool      ALEformulation; //!< Formulation switch (STANDARD or ALE)
  TimeIntegration::BDF bdf; //!< BDF helper class
  TimeIntegration::Method timeMethod; //!< Time integration method

  Vector  ux;       //!< Grid velocity (ALE)
};

#endif
