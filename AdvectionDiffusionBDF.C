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
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "StabilizationUtils.h"
#include "CFDenums.h"
#include "WeakOperators.h"


AdvectionDiffusionBDF::AdvectionDiffusionBDF (unsigned short int n, int order,
                                              int itg_type, int form) :
  AdvectionDiffusion(n, itg_type == STANDARD ? NONE : SUPG),
  formulation(form), bdf(order)
{
  primsol.resize(order);
  velocity.resize(order);
  registerVector("velocity1", &velocity[0]);
  if (order == 2)
    registerVector("velocity2", &velocity[1]);
  registerVector("grid velocity", &ux);
}


bool AdvectionDiffusionBDF::initElement(const std::vector<int>& MNPC,
                                        LocalIntegral& A)
{
  size_t nvec   = primsol.size() + velocity.size();
  size_t nfield = nvec;
  if (formulation & CFD::ALE)
    nfield++;

  A.vec.resize(nfield);
  int ierr = 0;
  for (size_t i = 0; i < primsol.size() && ierr == 0; i++) {
    ierr = utl::gather(MNPC,1,primsol[i],A.vec[i]);
    if (!Uad)
      ierr = utl::gather(MNPC,nsd,velocity[i],A.vec[i+primsol.size()]);
  }

  if (formulation & CFD::ALE)
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

  Vec3 U;
  if (Uad)
    U = (*Uad)(X);
  else {
    for (size_t i=1;i<=nsd;++i) {
      double tmp[2] = {0};
      for (int j=0; j<bdf.getOrder();++j)
        tmp[j] = elMat.vec[primsol.size()+j].dot(fe.N,i-1,nsd);
      U[i-1] = bdf.extrapolate(tmp);
    }
    if (formulation & CFD::ALE) {
      for (size_t i=1;i<=nsd;++i)
        U[i-1] -= elMat.vec[primsol.size()+velocity.size()].dot(fe.N,i-1,nsd);
    }
  }
  double react = 0;
  if (reaction)
    react = (*reaction)(X);

  double rhoC = props.getMassDensity()*props.getHeatCapacity();
  double theta=0;
  for (int t=0;t<bdf.getOrder();++t) {
    double val = fe.N.dot(elMat.vec[t]);
    theta += -rhoC*bdf[1+t]/time.dt*val;
  }

  double tau=0;
  if (stab == SUPG)
    tau = StabilizationUtils::getTauPt(time.dt, props.getDiffusivity(), Vector(U.ptr(),nsd), fe.G);

  // Integrate source, if defined
  if (source)
    theta += (*source)(X);

  WeakOperators::Laplacian(elMat.A[0], fe, props.getDiffusivity());
  WeakOperators::Mass(elMat.A[0], fe, bdf[0]/time.dt+react);

  // loop over test functions (i) and basis functions (j)
  for (size_t i = 1; i <= fe.N.size(); ++i) {
    double convI = 0.0;
    if (stab == SUPG) {
      // Convection for test functions
      for (size_t k = 1;k <= nsd;k++)
        convI += U[k-1]*fe.dNdX(i,k);
      convI *= tau*fe.detJxW;
    }
    for (size_t j = 1; j <= fe.N.size(); ++j) {
      double advect = 0.0;
      for (size_t k = 1;k <= nsd; ++k)
        advect += U[k-1]*fe.dNdX(j,k);
      advect *= fe.N(i)*rhoC;

      elMat.A[0](i,j) += advect*fe.detJxW;

      if (stab == SUPG) {
        double laplace = 0.0;
        for (size_t k = 1;k <= nsd;k++) {
          // Diffusion
          laplace -= fe.d2NdX2(j,k,k);
        }
        elMat.eMs(i,j) += (props.getDiffusivity()*laplace+advect+
                           (rhoC*bdf[0]/time.dt+react)*fe.N(j))*convI;
      }
    }
    if (stab == SUPG)
      elMat.eSs(i) += theta*convI;
  }

  WeakOperators::Source(elMat.b.front(), fe, theta);

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
