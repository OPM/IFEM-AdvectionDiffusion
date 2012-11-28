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
#include "CFDenums.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "StabilizationUtils.h"
#include "TimeDomain.h"
#include "Vec3Oper.h"


AdvectionDiffusionExplicit::AdvectionDiffusionExplicit (unsigned short int n,
                                                        int itg_type, int form) :
  AdvectionDiffusion(n, itg_type == Integrand::STANDARD?NONE:SUPG)
{
  primsol.resize(1);
}


AdvectionDiffusionExplicit::~AdvectionDiffusionExplicit ()
{
}


LocalIntegral* AdvectionDiffusionExplicit::getLocalIntegral (size_t nen, 
                                                             size_t,
                                                             bool neumann) const
{
  ElementInfo* result = new ElementInfo;

  if (neumann) {
    result->withLHS = false;
    result->resize(0,1);
    result->b[0].resize(nen,true);
  } else {
    result->withLHS = true;
    result->resize(2,1);
    result->A[0].resize(nen,nen,true);
    result->A[1].resize(nen,nen,true);
    result->b[0].resize(nen,true);
  }

  return result;
}


bool AdvectionDiffusionExplicit::evalInt (LocalIntegral& elmInt,
                                          const FiniteElement& fe,
                                          const TimeDomain& time,
                                          const Vec3& X) const
{
  ElementInfo& elMat = static_cast<ElementInfo&>(elmInt);

  Vector U;
  U.resize(nsd);
  double nut = kappa;
  Vec3 Ua;
  if (Uad)
    Ua = (*Uad)(X);
  double react=0.0;
  if (reaction)
    react = (*reaction)(X);

  // Integrate source, if defined
  double theta=0.0;
  if (source)
    theta += (*source)(X);

  // loop over test functions (i) and basis functions (j)
  for (size_t i = 1; i <= fe.N.size(); ++i) {
    for (size_t j = 1; j <= fe.N.size(); ++j) {
      double laplace = 0.0, advect = 0.0;
      for (size_t k = 1;k <= nsd; ++k) {
        laplace += fe.dNdX(i,k)*fe.dNdX(j,k);
        advect += U[k-1]*fe.dNdX(j,k);
      }
      advect *= fe.N(i);

      elMat.A[1](i,j) -= (nut*laplace+advect+
                         react*fe.N(i)*fe.N(j))*fe.detJxW;
      elMat.A[0](i,j) += fe.N(i)*fe.N(j)*fe.detJxW;
    }
  }

  elMat.b.front().add(fe.N, theta*fe.detJxW);

  return true;
}


bool AdvectionDiffusionExplicit::finalizeElement(LocalIntegral& A, 
                                                 const TimeDomain&, size_t)
{
  ElementInfo& elMat = static_cast<ElementInfo&>(A);
  elMat.A[1].multiply(elMat.vec[0], elMat.b[0], false, true);
  return true;
}
