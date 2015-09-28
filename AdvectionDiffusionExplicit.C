// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionExplicit.C
//!
//! \date Oct 29 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for time-dependent Advection-Diffusion
//!        problems.
//!
//==============================================================================

#include "AdvectionDiffusionExplicit.h"
#include "FiniteElement.h"
#include "WeakOperators.h"


AdvectionDiffusionExplicit::AdvectionDiffusionExplicit (unsigned short int n,
                                                        int itg_type, int form) :
  AdvectionDiffusion(n, itg_type == Integrand::STANDARD?NONE:SUPG)
{
  primsol.resize(1);
}


LocalIntegral* AdvectionDiffusionExplicit::getLocalIntegral (size_t nen, size_t,
                                                             bool neumann) const
{
  ElmMats* result = new ElmMats(!neumann);
  result->resize(neumann ? 0 : 2, 1);
  result->redim(nen);

  return result;
}


bool AdvectionDiffusionExplicit::evalInt (LocalIntegral& elmInt,
                                          const FiniteElement& fe,
                                          const TimeDomain& time,
                                          const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  Vec3 U;
  if (Uad) U = (*Uad)(X);
  double nut = props.getDiffusivity();
  double react = reaction ? (*reaction)(X) : 0.0;

  // Integrate source, if defined
  if (source)
    WeakOperators::Source(elMat.b.front(), fe, (*source)(X));

  WeakOperators::Laplacian(elMat.A[1], fe, -nut);
  WeakOperators::Mass(elMat.A[0], fe, 1.0);
  WeakOperators::Mass(elMat.A[1], fe, -react);
  WeakOperators::Advection(elMat.A[1], fe, U, -1.0);

  return true;
}


bool AdvectionDiffusionExplicit::finalizeElement (LocalIntegral& A)
{
  ElmMats& elMat = static_cast<ElmMats&>(A);
  return elMat.A[1].multiply(elMat.vec[0], elMat.b[0], false, true);
}
