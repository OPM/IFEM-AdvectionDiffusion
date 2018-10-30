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
#include "FiniteElement.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "Utilities.h"


AdvectionDiffusion::AdvectionDiffusion (unsigned short int n,
                                        AdvectionDiffusion::Stabilization s)
  : IntegrandBase(n), order(1), stab(s), Cinv(5.0)
{
  primsol.resize(1);

  Uad = nullptr;
  reaction = source = flux = nullptr;
}


AdvectionDiffusion::~AdvectionDiffusion()
{
  delete Uad;
  delete reaction;
  delete source;
}


int AdvectionDiffusion::getIntegrandType () const
{
  return stab == NONE ? STANDARD : SECOND_DERIVATIVES | ELEMENT_CORNERS;
}


LocalIntegral* AdvectionDiffusion::getLocalIntegral (size_t nen, size_t,
                                                     bool neumann) const
{
  ElementInfo* result = new ElementInfo(!neumann);
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
    Vec3 U;
    if (Uad)
      U = (*Uad)(X);

    WeakOps::Laplacian(elMat.A[0], fe, props.getDiffusivity());
    WeakOps::Mass(elMat.A[0], fe, react);
    WeakOps::Advection(elMat.A[0], fe, U, 1.0);

    // loop over test functions (i) and basis functions (j)
    for (size_t i = 1; i <= fe.N.size(); i++)
      for (size_t j = 1; j <= fe.N.size(); j++)
        if (stab == SUPG) {
          double convV=0, Lu=0;
          for (size_t k = 1;k <= nsd;k++) {
            convV += U[k-1]*fe.dNdX(i,k);
            Lu += U[k-1]*fe.dNdX(j,k)-props.getDiffusivity()*fe.d2NdX2(j,k,k);
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
            Lv += U[k-1]*fe.dNdX(i,k)-props.getDiffusivity()*fe.d2NdX2(i,k,k);
            Lu += U[k-1]*fe.dNdX(j,k)-props.getDiffusivity()*fe.d2NdX2(j,k,k);
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
            Lav += -U[k-1]*fe.dNdX(i,k)-props.getDiffusivity()*fe.d2NdX2(i,k,k);
            Lu += U[k-1]*fe.dNdX(j,k)-props.getDiffusivity()*fe.d2NdX2(j,k,k);
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


bool AdvectionDiffusion::evalSol (Vector& s, const FiniteElement& fe,
                                  const Vec3&,
                                  const std::vector<int>& MNPC) const
{
  Vector ePhi;
  if (utl::gather(MNPC,1,primsol.front(),ePhi) > 0)
    return false;

  return fe.dNdX.multiply(ePhi,s,true);
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

  if (mode >= SIM::RECOVERY)
    primsol.resize(1);
  else
    primsol.clear();
}


bool AdvectionDiffusion::finalizeElement (LocalIntegral& A)
{
  if (stab == NONE)
    return true;

  ElementInfo& E = static_cast<ElementInfo&>(A);

  // Compute stabilization parameter
  double tau = E.getTau(props.getDiffusivity(), Cinv, order);

  E.eMs *= tau;
  E.eSs *= tau;

  tauE(E.iEl) = tau;

  // Add stabilization terms
  E.A[0] += E.eMs;
  E.b[0] += E.eSs;

  return true;
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


AdvectionDiffusion::WeakDirichlet::WeakDirichlet (unsigned short int n,
                                                  double CBI_, double gamma_)
  : IntegrandBase(n), CBI(CBI_), gamma(gamma_), Uad(nullptr), flux(nullptr)
{
  // Need current solution only
  primsol.resize(1);
}


AdvectionDiffusion::WeakDirichlet::~WeakDirichlet ()
{
  delete Uad;
}


LocalIntegral* AdvectionDiffusion::WeakDirichlet::getLocalIntegral (size_t nen,
								    size_t,
								    bool) const
{
  ElmMats* result = new ElmMats;
  result->resize(1,1);
  result->redim(nen);

  return result;
}


bool AdvectionDiffusion::WeakDirichlet::initElementBou (const std::vector<int>& MNPC,
                                                        LocalIntegral& A)
{
  int ierr = 0;
  size_t nvec = primsol.size();
  for (size_t i = 0; i < nvec && !primsol[i].empty() && ierr == 0; i++)
    ierr = utl::gather(MNPC,1,primsol[i],A.vec[i]);

  if (ierr == 0) return true;

  std::cerr <<" *** AdvectionDiffusion::WeakDirichlet::initElementBou: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


bool AdvectionDiffusion::WeakDirichlet::evalBou (LocalIntegral& elmInt,
					         const FiniteElement& fe,
					         const Vec3& X,
					         const Vec3& normal) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // evaluate advection field
  Vec3 U;
  if (Uad)
    U = (*Uad)(X);

  double g = flux ? (*flux)(X) : 0.0;
  double An = U*normal;
  double C = CBI*fabs(props.getDiffusivity())/fe.h;
  double kap = props.getDiffusivity();

  // loop over test functions (i) and basis functions (j)
  for (size_t i = 1; i <= fe.N.size(); ++i) {
    double addI = 0.0;
    for (int k=1; k <= nsd; ++k)
      addI += fe.dNdX(i,k)*normal[k-1];
    addI *= gamma*kap;
    for (size_t j = 1; j <= fe.N.size(); ++j) {
      // adjoint
      double addJ = 0.0;
      for (int k=1; k <= nsd; ++k)
        addJ += fe.dNdX(j,k)*normal[k-1];
      elMat.A[0](i,j) += (-kap*addJ+An*fe.N(j))*fe.N(i)*fe.detJxW;

      // inflow
      if (An < 0)
        elMat.A[0](i,j) += fe.N(j)*(-addI-An*fe.N(i))*fe.detJxW;

      // outflow
      if (An > 0)
        elMat.A[0](i,j) += -addI*fe.N(j)*fe.detJxW;

      elMat.A[0](i,j) += C*fe.N(i)*fe.N(j)*fe.detJxW;
    }
    if (flux) {
      // inflow
      if (An < 0)
        elMat.b[0](i) += -g*(addI+An*fe.N(i))*fe.detJxW;

      // outflow
      if (An > 0)
        elMat.b[0](i) += -g*addI*fe.detJxW;

      elMat.b[0](i) += C*g*fe.N(i)*fe.detJxW;
    }
  }

  return true;
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
  pnorm[ip++] += kappa*gradUh.dot(gradUh)*fe.detJxW;

  Vec3 dT;
  if (anasol && anasol->getScalarSecSol()) {
    dT = (*anasol->getScalarSecSol())(X);
    pnorm[ip++] += kappa*dT*dT*fe.detJxW;
    pnorm[ip++] += kappa*(dT-gradUh)*(dT-gradUh)*fe.detJxW;
  }

  if (anasol && anasol->getScalarSol()) {
    double T = (*anasol->getScalarSol())(X);
    pnorm[ip++] += T*T*fe.detJxW; // L2 norm of analytical solution
    pnorm[ip++] += (T-Uh)*(T-Uh)*fe.detJxW; // L2 norm of error
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
        double res = -kappa*hess.sum() + U*gradUh + react*Uh - f;
        double kk = fe.h*std::min(1.0/(sqrt(2.0)*sqrt(13.0)),fe.h/(3.0*sqrt(10.0)*hep.props.getDiffusivity()));
        ip++; // unused
        pnorm[ip++] += kk*kk*res*res*fe.detJxW;
        if (anasol && anasol->getScalarSecSol())
          ip += 2;
      }
    } else {
      Vec3 rGrad;
      for (size_t k = 0; k < hep.nsd; k++)
        rGrad[k] += psol.dot(fe.N,k,hep.nsd);

      pnorm[ip++] += rGrad*rGrad*fe.detJxW;
      // Recovery based estimate
      pnorm[ip++] += (rGrad-gradUh)*(rGrad-gradUh)*fe.detJxW;

      if (anasol && anasol->getScalarSecSol()) {
        pnorm[ip++] += (dT-rGrad)*(dT-rGrad)*fe.detJxW;
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
    "|T^h|_res",
    "effectivity index"
  };


  const char** n = s;
  if (i > 1) {
    size_t nNrm = this->getNoFields(i);
    if (nNrm == 1 || (nNrm == 2 && anasol))
      n = res;
    else
      n = rec;
  }

  if (!prefix)
    return n[j-1];

  return prefix + std::string(" ") + n[j-1];
}
