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
#include "WeakOperators.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "Vec3Oper.h"


AdvectionDiffusionBDF::AdvectionDiffusionBDF (unsigned short int n, int order,
                                              int itg_type, bool useALE)
  : AdvectionDiffusion(n, itg_type == STANDARD ? NONE : SUPG), bdf(order)
{
  ALEformulation = useALE;
  primsol.resize(1+order);
  velocity.resize(order);
  registerVector("velocity1", &velocity[0]);
  if (order == 2)
    registerVector("velocity2", &velocity[1]);
  registerVector("grid velocity", &ux);
}


bool AdvectionDiffusionBDF::initElement (const std::vector<int>& MNPC,
                                         LocalIntegral& A)
{
  size_t nvec   = primsol.size() + velocity.size();
  size_t nfield = nvec;
  if (ALEformulation)
    nfield++;

  A.vec.resize(nfield);
  int ierr = 0;
  for (size_t i = 0; i < primsol.size() && ierr == 0; i++) {
    ierr = utl::gather(MNPC,1,primsol[i],A.vec[i]);
    if (!Uad && i < primsol.size()-1)
      ierr = utl::gather(MNPC,nsd,velocity[i],A.vec[i+primsol.size()]);
  }

  if (ALEformulation)
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
    if (ALEformulation)
      for (size_t i=1;i<=nsd;++i)
        U[i-1] -= elMat.vec[primsol.size()+velocity.size()].dot(fe.N,i-1,nsd);
  }
  double react = 0;
  if (reaction)
    react = (*reaction)(X);

  double theta=0;
  for (int t=0;t<bdf.getOrder();++t) {
    double val = fe.N.dot(elMat.vec[t+1]);
    theta += -props.getMassAdvectionConstant()*bdf[1+t]/time.dt*val;
  }

  double tau=0;
  if (stab == SUPG)
    tau = StabilizationUtils::getTauPt(time.dt, props.getDiffusivity(),
                                       Vector(U.ptr(),nsd), fe.G);

  // Integrate source, if defined
  if (source)
    theta += (*source)(X);

  WeakOperators::Laplacian(elMat.A[0], fe, props.getDiffusionConstant());
  WeakOperators::Mass(elMat.A[0], fe, props.getMassAdvectionConstant()*bdf[0]/time.dt +
                                      react*props.getReactionConstant());

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
      advect *= fe.N(i)*props.getMassAdvectionConstant();

      elMat.A[0](i,j) += advect*fe.detJxW;

      if (stab == SUPG) {
        double laplace = 0.0;
        for (size_t k = 1;k <= nsd;k++) {
          // Diffusion
          laplace -= fe.d2NdX2(j,k,k);
        }
        elMat.eMs(i,j) += (props.getDiffusionConstant()*laplace+advect+
                           (props.getMassAdvectionConstant()*bdf[0]/time.dt+
                            props.getReactionConstant()*react)*fe.N(j))*convI;
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

ADNorm::ADNorm(AdvectionDiffusion& p, AnaSol* a) : NormBase(p)
{
  nrcmp = myProblem.getNoFields(2);
  anasol = a;
}


bool ADNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                      const Vec3& X) const
{
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  AdvectionDiffusion& hep = static_cast<AdvectionDiffusion&>(myProblem);

  // Evaluate the FE temperature and thermal conductivity at current point
  double Uh = fe.N.dot(elmInt.vec.front());
  double kappa = hep.getFluidProperties().getDiffusionConstant();

  // Evaluate the FE heat flux vector, gradU = dNdX^T * eV
  Vector gradUh;
  if (!fe.dNdX.multiply(elmInt.vec.front(),gradUh,true))
    return false;

  size_t ip = 0;
  // Integrate the L2 norm, (U^h, U^h)
  pnorm[ip++] += Uh*Uh*fe.detJxW;

  // Integrate the energy norm, a(U^h,U^h)
  pnorm[ip++] = 0.5*kappa*gradUh.dot(gradUh)*fe.detJxW;

  if (anasol && anasol->getScalarSol()) {
    double T = (*anasol->getScalarSol())(X);
    pnorm[ip++] += T*T*fe.detJxW; // L2 norm of analytical solution
    pnorm[ip++] += (T-Uh)*(T-Uh)*fe.detJxW; // L2 norm of error
  }

  if (anasol && anasol->getScalarSecSol()) {
    Vec3 dT = (*anasol->getScalarSecSol())(X);
    pnorm[ip++] += 0.5*kappa*dT*dT*fe.detJxW;
    pnorm[ip++] += 0.5*kappa*(dT-gradUh)*(dT-gradUh)*fe.detJxW;
  }

  return true;
}


size_t ADNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else
    return anasol ? 6 : 2;
}


std::string ADNorm::getName (size_t i, size_t j,
                             const char* prefix) const
{
  if (i == 0 || j == 0 || j > 4)
    return this->NormBase::getName(i,j,prefix);

  static const char* s[] = {
    "(theta^h,theta^h)^0.5",
    "a(theta^h,theta^h)^0.5",
    "(q,theta^h)^0.5",
    "(theta, theta)^0.5",
    "a(theta,theta)^0.5",
    "a(e,e)^0.5, e=theta-theta^h",
    "a(theta^r,theta^r)^0.5",
    "a(e,e)^0.5, e=theta^r-theta^h",
    "a(e,e)^0.5, e=theta-theta^r",
    "effectivity index"
  };

  size_t k = i > 1 ? j+3 : j-1;

  if (!prefix)
    return s[k];

  return prefix + std::string(" ") + s[k];
}


bool ADNorm::hasElementContributions (size_t i, size_t j) const
{
  return i > 1 || j != 2;
}


NormBase* AdvectionDiffusionBDF::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new ADNorm(*const_cast<AdvectionDiffusionBDF*>(this),asol);
  else
    return new ADNorm(*const_cast<AdvectionDiffusionBDF*>(this));
}
