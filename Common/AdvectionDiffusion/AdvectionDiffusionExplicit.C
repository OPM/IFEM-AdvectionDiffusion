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


AdvectionDiffusionExplicit::AdvectionDiffusionExplicit (unsigned short int n,
                                                        TimeIntegration::Method method,
                                                        int itg_type, int form) :
  AdvectionDiffusion(n, itg_type == Integrand::STANDARD?NONE:SUPG),
  timeMethod(method)
{
  primsol.resize(1 + TimeIntegration::Steps(method));
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
                                          const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Integrate source, if defined
  if (source)
    WeakOps::Source(elMat.b.front(), fe, (*source)(X));

  WeakOps::Laplacian(elMat.A[1], fe, -props.getDiffusivity());
  WeakOps::Mass(elMat.A[0], fe, 1.0);
  if (reaction)
    WeakOps::Mass(elMat.A[1], fe, -(*reaction)(X));

  // Integrate advection, if defined
  if (Uad)
    WeakOps::Advection(elMat.A[1], fe, (*Uad)(X), -1.0);

  return true;
}


bool AdvectionDiffusionExplicit::finalizeElement (LocalIntegral& A)
{
  ElmMats& elMat = static_cast<ElmMats&>(A);
  return elMat.A[1].multiply(elMat.vec[0], elMat.b[0], false, true);
}


NormBase* AdvectionDiffusionExplicit::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new AdvectionDiffusionNorm(*const_cast<AdvectionDiffusionExplicit*>(this),asol);
  else
    return new AdvectionDiffusionNorm(*const_cast<AdvectionDiffusionExplicit*>(this));
}
