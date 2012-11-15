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
#include "CFDenums.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "StabilizationUtils.h"
#include "TimeDomain.h"
#include "Vec3Oper.h"


AdvectionDiffusionBDF::AdvectionDiffusionBDF (unsigned short int n, int order,
                                              int itg_type, int form) :
  AdvectionDiffusion(n, itg_type == Integrand::STANDARD?NONE:SUPG),
  formulation(form), bdf(order)
{
  primsol.resize(order);
  velocity.resize(order);
  registerVector("velocity1", &velocity[0]);
  if (order == 2)
    registerVector("velocity2", &velocity[1]);

  registerVector("nut", &nut);
  registerVector("grid velocity", &ux);
}


AdvectionDiffusionBDF::~AdvectionDiffusionBDF()
{
}


bool AdvectionDiffusionBDF::initElement(const std::vector<int>& MNPC,
                                        LocalIntegral& A)
{
  size_t nvec   = primsol.size() + velocity.size();
  if (formulation & CFD::RANS)
    nvec++;
//  if (formulation & SIM::ALE)
//    nvec++;

  A.vec.resize(nvec);
  int ierr = 0;
  for (size_t i = 0; i < primsol.size() && ierr == 0; i++) {
    ierr = utl::gather(MNPC,1,primsol[i],A.vec[i]);
    if (!Uad)
      ierr = utl::gather(MNPC,nsd,velocity[i],A.vec[i+primsol.size()]);
  }

  if (formulation & CFD::RANS)
    ierr = utl::gather(MNPC,1,nut,A.vec[nvec-1]);

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

  Vector U;
  U.resize(nsd);
  double nut = kappa;
  if (Uad) {
    Vec3 Ua = (*Uad)(X);
    for (int i=0;i<nsd;++i)
      U[i] = Ua[i];
  } else {
    for (size_t i=1;i<=nsd;++i) {
      double tmp[2] = {0};
      for (int j=0; j<bdf.getOrder();++j)
        tmp[j] = elMat.vec[primsol.size()+j].dot(fe.N,i-1,nsd);
      U[i-1] = bdf.extrapolate(tmp);
    }
    if (formulation & CFD::RANS)
      nut += fe.N.dot(elMat.vec[2*primsol.size()])/Pr;
  }
  double react = 0;
  if (reaction)
    react = (*reaction)(X);

  double tau=0;
  if (stab == SUPG)
    tau = StabilizationUtils::getTauPt(time.dt, nut, U, fe.G);

  double theta=0;
  for (int t=0;t<bdf.getOrder();++t)
    theta += -bdf[1+t]/time.dt*fe.N.dot(elMat.vec[t]);

  // Integrate source, if defined
  if (source)
    theta += (*source)(X);

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
      double laplace = 0.0, advect = 0.0;
      for (size_t k = 1;k <= nsd; ++k) {
        laplace += fe.dNdX(i,k)*fe.dNdX(j,k);
        advect += U[k-1]*fe.dNdX(j,k);
      }
      advect *= fe.N(i);

      elMat.A[0](i,j) += (nut*laplace+advect+
                          (bdf[0]/time.dt+react)*fe.N(i)*fe.N(j))*fe.detJxW;

      if (stab == SUPG) {
        laplace = 0.0;
        for (size_t k = 1;k <= nsd;k++) {
          // Diffusion
          laplace -= fe.d2NdX2(j,k,k);
        }
        elMat.eMs(i,j) += (nut*laplace+advect+
                           (bdf[0]/time.dt+react)*fe.N(j))*convI;
      }
    }
    if (stab == SUPG)
      elMat.eSs(i) += theta*convI;
  }

  elMat.b.front().add(fe.N, theta*fe.detJxW);

  return true;
}


bool AdvectionDiffusionBDF::finalizeElement(LocalIntegral& A, 
                                            const TimeDomain&, size_t)
{
  if (stab != NONE) {
    ElementInfo& E = static_cast<ElementInfo&>(A);

    // Add stabilization terms
    E.A[0] += E.eMs;
    E.b[0] += E.eSs;
  }

  return true;
}
