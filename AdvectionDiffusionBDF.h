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

#ifndef _ADVECTIONDIFFUSIONBDF_H
#define _ADVECTIONDIFFUSIONBDF_H

#include "AdvectionDiffusion.h"
#include "BDF.h"
#include "ElmMats.h"
#include "Vec3.h"


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
  //! \param[in] order Temporal order (1,2)
  //! \param[in] stab Integrand formulation
  AdvectionDiffusionBDF(unsigned short int n = 3, int order = 1,
                        int itg_type_=Integrand::STANDARD,
                        int form = 0);

  //! \brief Empty destructor.
  virtual ~AdvectionDiffusionBDF();

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& A);

  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents.
  //! It can also be used to implement multiple integration point loops within
  //! the same element, provided the necessary integration point values are
  //! stored internally in the object during the first integration loop.
  virtual bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t = 0);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const TimeDomain& time, const Vec3& X) const;

  //! \brief Returns the integrand type
  virtual int getIntegrandType() const
  { 
    if (stab == NONE)
      return Integrand::STANDARD;

    return Integrand::SECOND_DERIVATIVES | Integrand::G_MATRIX;
  }

  virtual void advanceStep()
  {
    bdf.advanceStep();
  }
protected:
  int formulation;  //!< Formulation (STANDARD, RANS, ALE)
  BDF bdf;          //!< BDF helper class

  Vectors velocity; //!< The advecting velocity field
  Vector  nut;      //!< The turbulent viscosity field
  Vector  ux;       //!< Grid velocity (ALE)
};

#endif
