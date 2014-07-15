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
  double nut = kappa;
  double react = reaction ? (*reaction)(X) : 0.0;

  // Integrate source, if defined
  if (source)
    elMat.b.front().add(fe.N,(*source)(X)*fe.detJxW);

  // loop over test functions (i) and basis functions (j)
  for (size_t i = 1; i <= fe.N.size(); ++i) {
    for (size_t j = 1; j <= fe.N.size(); ++j) {
      double laplace = 0.0, advect = 0.0;
      for (size_t k = 1;k <= nsd; ++k) {
        laplace += fe.dNdX(i,k)*fe.dNdX(j,k);
        advect += U[k-1]*fe.dNdX(j,k);
      }
      elMat.A[1](i,j) -= (nut*laplace+fe.N(i)*(advect+react*fe.N(j)))*fe.detJxW;
      elMat.A[0](i,j) += fe.N(i)*fe.N(j)*fe.detJxW;
    }
  }

  return true;
}


bool AdvectionDiffusionExplicit::finalizeElement (LocalIntegral& A)
{
  ElmMats& elMat = static_cast<ElmMats&>(A);
  return elMat.A[1].multiply(elMat.vec[0], elMat.b[0], false, true);
}
