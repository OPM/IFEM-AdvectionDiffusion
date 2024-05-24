// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusion.C
//!
//! \date Jun 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Advection-Diffusion problems.
//!
//==============================================================================

#include "AdvectionDiffusion.h"

#include "AnaSol.h"
#include "ElmNorm.h"
#include "FiniteElement.h"
#include "Function.h"
#include "Fields.h"
#include "Integrand.h"
#include "LocalIntegral.h"
#include "Vec3.h"
#include "Vec3Oper.h"

#include <algorithm>
#include <cmath>
#include <ext/alloc_traits.h>
#include <iostream>
#include <vector>


AdvectionDiffusion::AdvectionDiffusion (unsigned short int n,
                                        AdvectionDiffusion::Stabilization s)
  : IntegrandBase(n), order(1), stab(s), Cinv(5.0), residualNorm(false)
{
  primsol.resize(1);
  flux = nullptr;

  velocity.resize(2);
  registerVector("velocity1", &velocity[0]);
  registerVector("velocity2", &velocity[1]);
}


AdvectionDiffusion::~AdvectionDiffusion()
{
}


int AdvectionDiffusion::getIntegrandType () const
{
  return stab == NONE ? STANDARD : SECOND_DERIVATIVES | ELEMENT_CORNERS;
}


LocalIntegral* AdvectionDiffusion::getLocalIntegral (size_t nen, size_t,
                                                     bool neumann) const
{
  ElementInfo* result = new ElementInfo(!neumann);
  result->rhsOnly = m_mode >= SIM::RHS_ONLY;
  result->resize(neumann ? 0 : 1, 1);
  result->redim(nen);

  if (stab != NONE) {
    result->eMs.resize(nen,nen);
    result->eSs.resize(nen);
    result->Cv.resize(nsd+1);
  }

  return result;
}


bool AdvectionDiffusion::evalInt (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const Vec3& X) const
{
  ElementInfo& elMat = static_cast<ElementInfo&>(elmInt);
  elMat.iEl = fe.iel;
  elMat.hk = fe.h;

  double f = source ? (*source)(X) : 0.0;

  if (!elMat.A.empty()) {
    // evaluate reaction field
    double react = reaction ? (*reaction)(X) : 0.0;

    // evaluate advection field
    Vec3 U = this->getAdvectionVelocity(fe, X);

    WeakOps::Laplacian(elMat.A[0], fe, props.getDiffusionConstant(X));
    WeakOps::Mass(elMat.A[0], fe, react);
    WeakOps::Advection(elMat.A[0], fe, U, props.getMassAdvectionConstant(), advForm);

    // loop over test functions (i) and basis functions (j)
    for (size_t i = 1; i <= fe.N.size(); i++)
      for (size_t j = 1; j <= fe.N.size(); j++)
        if (stab == SUPG) {
          double convV=0, Lu=0;
          for (size_t k = 1;k <= nsd;k++) {
            convV += U[k-1]*fe.dNdX(i,k);
            Lu += U[k-1]*fe.dNdX(j,k)-props.getDiffusionConstant(X)*fe.d2NdX2(j,k,k);
            elMat.Cv(k) += U[k-1]*fe.detJxW;
          }
          Lu += react*fe.N(j);
          elMat.Cv(nsd+1) += fe.detJxW;
          elMat.eMs(i,j) += convV*Lu*fe.detJxW;
          if (source && j == 1)
            elMat.eSs(i) += convV*f*fe.detJxW;
        }

        else if (stab == GLS) {
          double Lv=0, Lu=0;
          for (size_t k = 1;k <= nsd;k++) {
            Lv += U[k-1]*fe.dNdX(i,k)-props.getDiffusionConstant(X)*fe.d2NdX2(i,k,k);
            Lu += U[k-1]*fe.dNdX(j,k)-props.getDiffusionConstant(X)*fe.d2NdX2(j,k,k);
            elMat.Cv(k) += U[k-1]*fe.detJxW;
          }
          Lv += react*fe.N(i);
          Lu += react*fe.N(j);
          elMat.Cv(nsd+1) += fe.detJxW;
          elMat.eMs(i,j) += Lv*Lu*fe.detJxW;

          if (source && j == 1)
            elMat.eSs(i) += Lv*f*fe.detJxW;
        }

        else if (stab == MS) {
          double Lav=0, Lu=0;
          for (size_t k = 1;k <= nsd;k++) {
            Lav += -U[k-1]*fe.dNdX(i,k)-props.getDiffusionConstant(X)*fe.d2NdX2(i,k,k);
            Lu += U[k-1]*fe.dNdX(j,k)-props.getDiffusionConstant(X)*fe.d2NdX2(j,k,k);
            elMat.Cv(k) += U[k-1]*fe.detJxW;
          }
          Lav += react*fe.N(i);
          Lu +=  react*fe.N(j);
          elMat.Cv(nsd+1) += fe.detJxW;
          elMat.eMs(i,j) += -Lav*Lu*fe.detJxW;
          if (source && j == 1)
            elMat.eSs(i) += -Lav*f*fe.detJxW;
        }
  }

  // Integrate source, if defined
  if (source)
    WeakOps::Source(elMat.b.front(), fe, f);

  return true;
}


bool AdvectionDiffusion::evalBou (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const Vec3& X, const Vec3& normal) const
{
  if (!flux) {
    std::cerr <<" *** AdvectionDiffusion::evalBou: No fluxes."<< std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  if (elMat.b.empty()) {
    std::cerr <<" *** AdvectionDiffusion::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Evaluate the Neumann value
  double T = (*flux)(X);

  // Integrate the Neumann value
  WeakOps::Source(elMat.b.front(), fe, T);

  return true;
}


bool AdvectionDiffusion::evalSol2 (Vector& s, const Vectors& eV,
                                   const FiniteElement& fe, const Vec3&) const
{
  return eV.empty() ? false : fe.dNdX.multiply(eV.front(),s,true);
}


std::string AdvectionDiffusion::getField1Name (size_t, const char* prefix) const
{
  if (!prefix)
    return "T";

  return prefix + std::string(" T");
}


std::string AdvectionDiffusion::getField2Name (size_t i,
                                               const char* prefix) const
{
  if (i > 2)
    return "";

  static const char* s[3] = { "T,x", "T,y", "T,z" };
  if (!prefix)
    return s[i];

  return prefix + std::string(" ") + s[i];
}


void AdvectionDiffusion::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  primsol.resize(1);
}


bool AdvectionDiffusion::finalizeElement (LocalIntegral& A)
{
  if (stab == NONE)
    return true;

  ElementInfo& E = static_cast<ElementInfo&>(A);

  // Compute stabilization parameter
  double tau = E.getTau(props.getDiffusionConstant(Vec3()), Cinv, order);

  E.eMs *= tau;
  E.eSs *= tau;

  tauE(E.iEl) = tau;

  // Add stabilization terms
  E.A[0] += E.eMs;
  E.b[0] += E.eSs;

  return true;
}


void AdvectionDiffusion::setNamedFields (const std::string& name, Fields* field)
{
  if (name == "velocity1")
    uFields[0].reset(field);
  else
    uFields[1].reset(field);
}


Vec3 AdvectionDiffusion::getAdvectionVelocity (const FiniteElement& fe,
                                               const Vec3& X) const
{
  Vec3 U;
  if (uFields[0]) {
    Vector u;
    uFields[0]->valueFE(fe, u);
    U = u;
  } else if (Uad)
    U = (*Uad)(X);

 return U;
}


double AdvectionDiffusion::ElementInfo::getTau(double kappa,
                                               double Cinv, int p) const
{
  // eqns (211)-(214)
  double mk = std::min(1.0/(3.0*p*p), 2.0/Cinv);

  // find mean element advection velocity
  double vel = 0.0;
  for (size_t k = 1; k < Cv.size(); k++)
    vel += pow(Cv(k)/Cv.back(),2.0);
  vel = sqrt(vel);

  double Xi = std::min(mk*vel*hk/(2.0*kappa),1.0);

  return hk/(2.0*vel)*Xi;
}


NormBase* AdvectionDiffusion::getNormIntegrand (AnaSol* asol) const
{
  return new AdvectionDiffusionNorm(*const_cast<AdvectionDiffusion*>(this), asol);
}


AdvectionDiffusionNorm::AdvectionDiffusionNorm(AdvectionDiffusion& p, AnaSol* a) :
  NormBase(p)
{
  nrcmp = myProblem.getNoFields(2);
  anasol = a;
}


bool AdvectionDiffusionNorm::evalInt (LocalIntegral& elmInt,
                                      const FiniteElement& fe,
                                      const Vec3& X) const
{
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  const AdvectionDiffusion& hep = static_cast<const AdvectionDiffusion&>(myProblem);

  // Evaluate the FE temperature and thermal conductivity at current point
  double Th = fe.N.dot(elmInt.vec.front());
  double kappa = hep.getFluidProperties().getDiffusionConstant(X);

  // Evaluate the FE heat flux vector, gradU = dNdX^T * eV
  Vector dTh;
  if (!fe.dNdX.multiply(elmInt.vec.front(),dTh,true))
    return false;

  size_t ip = 0;
  // Integrate the L2 norm, (T^h, T^h)
  pnorm[ip++] += Th*Th*fe.detJxW;

  // Integrate the energy norm, a(T^h,T^h)
  pnorm[ip++] += kappa*dTh.dot(dTh)*fe.detJxW;

  Vec3 dT;
  if (anasol && anasol->getScalarSecSol()) {
    dT = (*anasol->getScalarSecSol())(X);
    pnorm[ip++] += kappa*dT*dT*fe.detJxW;
    pnorm[ip++] += kappa*(dT-dTh)*(dT-dTh)*fe.detJxW;
  }

  if (anasol && anasol->getScalarSol()) {
    double T = (*anasol->getScalarSol())(X);
    pnorm[ip++] += T*T*fe.detJxW; // L2 norm of analytical solution
    pnorm[ip++] += (T-Th)*(T-Th)*fe.detJxW; // L2 norm of error
  }

  for (const Vector& psol : pnorm.psol)
    if (psol.empty()) { // residual
      Vec3 hess;
      if (this->getIntegrandType() & Integrand::SECOND_DERIVATIVES) {
        for (size_t k = 1; k <= hep.getNoSpaceDim(); k++)
          for (size_t j = 1; j <= fe.N.size(); j++) {
            if (this->getIntegrandType() & Integrand::SECOND_DERIVATIVES)
              hess[k-1] += fe.d2NdX2(j,k,k)*pnorm.vec.front()(j);
        }
        double f = hep.source ? (*hep.source)(X) : 0.0;
        Vec3 U;
        if (hep.Uad)
          U = (*hep.Uad)(X);
        double react = hep.reaction ? (*hep.reaction)(X) : 0.0;
        double res = f + kappa*hess.sum() - U*dTh - react*Th;
        double kk;
        if (hep.useModifiedElmSize()) {
          if (hep.getCbar() > 0.0)
            kk = std::min(fe.h/sqrt(kappa), 1.0/sqrt(hep.getCbar()));
          else
            kk = fe.h/sqrt(kappa);
        } else
          kk = fe.h;
        ip++; // unused
        pnorm[ip++] += kk*kk*res*res*fe.detJxW;
        if (anasol && anasol->getScalarSecSol())
          ip += 2;
      }
    } else {
      Vec3 dTr;
      for (size_t k = 0; k < hep.nsd; k++)
        dTr[k] += psol.dot(fe.N,k,hep.nsd);

      pnorm[ip++] += kappa*dTr*dTr*fe.detJxW;
      // Recovery based estimate
      pnorm[ip++] += kappa*(dTr-dTh)*(dTr-dTh)*fe.detJxW;

      if (anasol && anasol->getScalarSecSol()) {
        pnorm[ip++] += kappa*(dT-dTr)*(dT-dTr)*fe.detJxW;
        ip++; // effectivity index
      }
    }

  return true;
}


bool AdvectionDiffusionNorm::finalizeElement (LocalIntegral& elmInt)
{
  if (!anasol) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as a(e^r,e^r)/a(e,e)
  // with e^r = u^r - u^h  and  e = u - u^h
  size_t g = 2;
  for (size_t ip = this->getNoFields(1); ip < pnorm.size(); ip += 4)
  {
    size_t gsize = this->getNoFields(g++);
    if (gsize != 4)
      continue;
    pnorm[ip+3] = pnorm[ip+1] / pnorm[3];
  }

  return true;
}


size_t AdvectionDiffusionNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else if (group == 1)
    return anasol ? 6 : 2;
  else
    return anasol ? 4 : 2;
}


std::string AdvectionDiffusionNorm::getName (size_t i, size_t j,
                                             const char* prefix) const
{
  static const char* s[] = {
    "(T^h,T^h)^0.5",
    "a(T^h,T^h)^0.5",
    "a(T,T)^0.5",
    "a(e,e)^0.5, e=T-T^h",
    "(T,T)^0.5",
    "(e,e)^0.5, e=T-T^h",
  };

  static const char* rec[] = {
    "a(T^r,T^r)^0.5",
    "a(e,e)^0.5, e=T^r-T^h",
    "a(e,e)^0.5, e=T-T^r",
    "effectivity index"
  };

  static const char* res[] = {
    "",
    "|T^h|_res",
    "",
    "effectivity index"
  };

  const AdvectionDiffusion& hep = static_cast<const AdvectionDiffusion&>(myProblem);
  const char** n = s;
  if (i > 1) {
    if (hep.doResidualNorm())
      n = res;
    else
      n = rec;
  }

  if (!prefix)
    return n[j-1];

  return prefix + std::string(" ") + n[j-1];
}


int AdvectionDiffusionNorm::getIntegrandType () const
{
  const AdvectionDiffusion& hep = static_cast<const AdvectionDiffusion&>(myProblem);
  return hep.doResidualNorm() ? SECOND_DERIVATIVES | ELEMENT_CORNERS : STANDARD;
}


AdvectionDiffusion::Robin::Robin(unsigned short int n, const AdvectionDiffusion& itg) :
  IntegrandBase(n),
  integrand(itg)
{
}


bool AdvectionDiffusion::Robin::evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                                        const Vec3& X, const Vec3& normal) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  if (integrand.getAdvForm() != WeakOperators::CONSERVATIVE) {
    std::cerr << "Robin boundary conditions only implemented for conservative form." << std::endl;
    return false;
  }

  Vec3 ax = alpha ? (*alpha)(X) : Vec3(1.0, 1.0, 1.0);
  if (g)
    ax[1] = (*g)(X);

  elMat.A[0].outer_product(fe.N, fe.N, true, fe.detJxW * ax[0]); // mass
  elMat.b[0].add(fe.N, fe.detJxW * ax[1]); // source

  return true;
}
