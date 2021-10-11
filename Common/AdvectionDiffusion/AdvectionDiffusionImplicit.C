// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionImplicit.C
//!
//! \date Aug 21 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementation for Advection-Diffusion problems with
//!        implicit time stepping.
//!
//==============================================================================

#include "AdvectionDiffusionImplicit.h"
#include "ADFluidProperties.h"

#include "ElmMats.h"
#include "FiniteElement.h"
#include "Function.h"
#include "MatVec.h"
#include "matrix.h"
#include "Vec3.h"
#include "Vec3Oper.h"

#include <cstddef>
#include <memory>
#include <vector>


class AnaSol;
class LocalIntegral;
class NormBase;


AdvectionDiffusionImplicit::AdvectionDiffusionImplicit (unsigned short int n,
                                                        TimeIntegration::Method method,
                                                        int itg_type, int form) :
  AdvectionDiffusion(n, itg_type == Integrand::STANDARD?NONE:SUPG),
  timeMethod(method)
{
  primsol.resize(2);
}


bool AdvectionDiffusionImplicit::evalInt (LocalIntegral& elmInt,
                                          const FiniteElement& fe,
                                          const TimeDomain& time,
                                          const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  double T = elMat.vec.empty()     ? 0.0 : elMat.vec[0].dot(fe.N);
  double Tp = elMat.vec.size() < 2 ? 0.0 : elMat.vec[1].dot(fe.N);

  Vec3 dTdX;
  for (size_t k = 0; k < nsd && !elMat.vec.empty(); ++k)
    dTdX[k] = elMat.vec[0].dot(fe.grad(1).getColumn(k+1));

  Vec3 vel;
  if (Uad)
    vel = (*Uad)(X);

  // Integrate source, if defined
  if (source)
    WeakOps::Source(elMat.b[0], fe, (*source)(X) * timeScale);

  if (!elMat.rhsOnly) {
    WeakOps::Mass(elMat.A[0], fe, 1.0); // Diagonal term
    WeakOps::Laplacian(elMat.A[0], fe, timeScale * props.getDiffusionConstant(X)); // Diffusion

    // Integrate advection, if defined
    if (Uad)
      WeakOps::Advection(elMat.A[0], fe, vel, timeScale * props.getMassAdvectionConstant()); // Advection
  }

  // Residual
  if (!elMat.rhsOnly) {
    WeakOps::Source(elMat.b[0], fe, -T); // Diagonal term
    WeakOps::Source(elMat.b[0], fe, Tp); // Diagonal term
  }

  ResidualOps::Laplacian(elMat.b[0], fe, dTdX, -timeScale*props.getDiffusionConstant(X)); // Diffusion
  WeakOps::Source(elMat.b[0], fe, -(vel*dTdX)*timeScale*props.getMassAdvectionConstant()); //Advection

  return true;
}


NormBase* AdvectionDiffusionImplicit::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new AdvectionDiffusionNorm(*const_cast<AdvectionDiffusionImplicit*>(this),asol);
  else
    return new AdvectionDiffusionNorm(*const_cast<AdvectionDiffusionImplicit*>(this));
}


void AdvectionDiffusionImplicit::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  if (mode == SIM::STATIC)
    primsol.clear();
  else
    primsol.resize(2);
}
