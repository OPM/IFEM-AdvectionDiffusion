// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionBDF.C
//!
//! \date Oct 29 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for time-dependent Advection-Diffusion
//!        problems.
//!
//==============================================================================

#include "AdvectionDiffusionBDF.h"

#include "Fields.h"
#include "FiniteElement.h"
#include "StabilizationUtils.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3Oper.h"

#include <ext/alloc_traits.h>
#include <iostream>


AdvectionDiffusionBDF::AdvectionDiffusionBDF (unsigned short int n,
                                              TimeIntegration::Method method,
                                              int itg_type, bool useALE)
  : AdvectionDiffusion(n, itg_type == STANDARD ? NONE : SUPG),
    bdf(TimeIntegration::Steps(method)),
    timeMethod(method)
{
  ALEformulation = useALE;
  primsol.resize(1+TimeIntegration::Steps(method));
  registerVector("grid velocity", &ux);
}


bool AdvectionDiffusionBDF::initElement (const std::vector<int>& MNPC,
                                         LocalIntegral& A)
{
  size_t nfield = primsol.size() + velocity.size();
  if (ALEformulation)
    nfield++;

  A.vec.resize(nfield);
  int ierr = 0;
  for (size_t i = 0; i < primsol.size() && ierr == 0; i++) {
    ierr = utl::gather(MNPC,1,primsol[i],A.vec[i]);
    if (!Uad && !uFields[0] && i < velocity.size() && ierr == 0)
      ierr = utl::gather(MNPC,nsd,velocity[i],A.vec[primsol.size()+i]);
  }

  if (ALEformulation && ierr == 0)
    ierr = utl::gather(MNPC,nsd,ux,A.vec[primsol.size()+velocity.size()]);

  if (ierr == 0)
    return true;

  std::cerr <<" *** AdvectionDiffusionBDF::initElement: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


bool AdvectionDiffusionBDF::evalInt (LocalIntegral& elmInt,
                                     const FiniteElement& fe,
                                     const TimeDomain& time,
                                     const Vec3& X) const
{
  ElementInfo& elMat = static_cast<ElementInfo&>(elmInt);

  Vec4 Xt(static_cast<const Vec4&>(X));
  if (timeMethod == TimeIntegration::THETA)
    Xt.t -= time.dt;

  std::array<Vec3,2> U;
  if (Uad)
    U = { (*Uad)(X), (*Uad)(Xt) };
  else if (uFields[0])
    for (int i = 0; i < 2 && uFields[i]; ++i) {
      Vector u;
      uFields[i]->valueFE(fe, u);
      U[i] = Vec3(u);
    }
  else {
    std::vector<Vec3> tmp(bdf.getOrder());
    for (int j=0; j<bdf.getOrder();++j)
      for (size_t i=0;i<nsd;++i)
        tmp[j][i] = elMat.vec[primsol.size()+j].dot(fe.N,i,nsd);
    U[0] = bdf.extrapolate(tmp);
  }

  if (ALEformulation)
    for (size_t i=0;i<nsd;++i)
      U[0][i] -= elMat.vec[primsol.size()+velocity.size()].dot(fe.N,i,nsd);

  double theta = 0.0;
  for (int t=0;t<bdf.getOrder();++t) {
    double val = fe.N.dot(elMat.vec[t+1]);
    theta += -props.getMassAdvectionConstant()*bdf[1+t]/time.dt*val;
  }

  double tau=0;
  if (stab == SUPG)
    tau = StabilizationUtils::getTauPt(time.dt, props.getDiffusivity(),
                                       Vector(U[0].ptr(),nsd), fe.G);

  // Integrate source, if defined
  if (source) {
    if (timeMethod == TimeIntegration::THETA)
      theta += 0.5*((*source)(X) + (*source)(Xt));
    else
      theta += (*source)(X);
  }

  double timeCoef = timeMethod == TimeIntegration::THETA ? 0.5 : 1.0;
  double mu = props.getDiffusionConstant(X)*timeCoef;
  double s = props.getMassAdvectionConstant()*bdf[0]/time.dt;
  double r = reaction ? props.getReactionConstant(X)*(*reaction)(X) : 0.0;

  WeakOps::Laplacian(elMat.A[0], fe, mu);
  WeakOps::Mass(elMat.A[0], fe, s + r*timeCoef);

  if (timeMethod == TimeIntegration::THETA) {
    RealArray g;
    fe.dNdX.multiply(elMat.vec[1], g, true);
    Vec3 grad(g);
    ResidualOps::Laplacian(elMat.b[0], fe, grad, mu);
    theta -= 0.5*props.getMassAdvectionConstant()*U[1]*grad;
    if (reaction)
      theta -= 0.5*props.getReactionConstant(X)*(*reaction)(Xt);
  }

  WeakOps::Source(elMat.b.front(), fe, theta);

  if (stab == NONE) {
    WeakOps::Advection(elMat.A[0], fe, U[0], timeCoef);
    return true;
  }

  // loop over test functions (i) and basis functions (j)
  for (size_t i = 1; i <= fe.N.size(); ++i) {
    double convI = 0.0;
    if (stab == SUPG) {
      // Convection for test functions
      for (size_t k = 1;k <= nsd;k++)
        convI += U[0][k-1]*fe.dNdX(i,k);
      convI *= tau*fe.detJxW;
    }
    for (size_t j = 1; j <= fe.N.size(); ++j) {
      double advect = 0.0;
      for (size_t k = 1;k <= nsd; ++k)
        advect += U[0][k-1]*fe.dNdX(j,k);
      advect *= fe.N(i)*props.getMassAdvectionConstant();

      elMat.A[0](i,j) += advect*fe.detJxW;

      if (stab == SUPG) {
        // Diffusion
        double laplace = -fe.d2NdX2.trace(j);
        elMat.eMs(i,j) += (props.getDiffusionConstant(X)*laplace+advect+
                           (props.getMassAdvectionConstant()*bdf[0]/time.dt+r)*
                            fe.N(j))*convI;
      }
    }
    if (stab == SUPG)
      elMat.eSs(i) += theta*convI;
  }

  return true;
}


bool AdvectionDiffusionBDF::finalizeElement (LocalIntegral& A)
{
  if (stab != NONE) {
    ElementInfo& E = static_cast<ElementInfo&>(A);

    // Add stabilization terms
    E.A[0] += E.eMs;
    E.b[0] += E.eSs;
  }

  return true;
}


NormBase* AdvectionDiffusionBDF::getNormIntegrand (AnaSol* asol) const
{
  return new AdvectionDiffusionNorm(*const_cast<AdvectionDiffusionBDF*>(this),asol);
}
