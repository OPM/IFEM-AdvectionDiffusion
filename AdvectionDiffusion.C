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
#include "Utilities.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "WeakOperators.h"


AdvectionDiffusion::AdvectionDiffusion (unsigned short int n,
                                        AdvectionDiffusion::Stabilization s) :
  nsd(n), order(1), stab(s), Cinv(5)
{
  primsol.resize(1);

  source = 0, Uad = 0, flux = 0, reaction = 0;
}


AdvectionDiffusion::~AdvectionDiffusion()
{
  delete source;
  delete Uad;
  delete reaction;
}


int AdvectionDiffusion::getIntegrandType () const
{
  if (stab == NONE)
    return STANDARD;

  return SECOND_DERIVATIVES | ELEMENT_CORNERS;
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


/*!
  \brief Returns the characteristic element size.
  \param[in] XC The element corner coordinates
  \param[in] nsd Number of space dimensions

  \details The size is taken as the longest diagonal.
*/

static double getElementSize (const Vec3Vec& XC, int nsd)
{
  Vec3Vec D(nsd*2-2);

  // Compute the element diagonals
  if (nsd == 2 && XC.size() >= 4) {
    D[0] = XC[3] - XC[0];
    D[1] = XC[1] - XC[2];
  }
  else if (nsd == 3 && XC.size() >= 8) {
    D[0] = XC[7] - XC[0];
    D[1] = XC[6] - XC[1];
    D[2] = XC[4] - XC[3];
    D[3] = XC[2] - XC[6];
  }
  else
    return 0.0;

  // Take our element size as the longest diagonal
  double h = D.front().length();
  for (size_t i = 1; i < D.size(); i++)
    h = std::max(h,D[i].length());

  return h;
}


bool AdvectionDiffusion::evalInt (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const Vec3& X) const
{
  ElementInfo& elMat = static_cast<ElementInfo&>(elmInt);
  if (stab != NONE)
    elMat.hk = getElementSize(fe.XC,nsd);
  elMat.iEl = fe.iel;

  double f = source ? (*source)(X) : 0.0;

  if (!elMat.A.empty()) {
    // evaluate reaction field
    double react = reaction ? (*reaction)(X) : 0.0;

    // evaluate advection field
    Vec3 U;
    if (Uad)
      U = (*Uad)(X);

    WeakOperators::Laplacian(elMat.A[0], fe, props.getDiffusivity());
    WeakOperators::Mass(elMat.A[0], fe, react);
    WeakOperators::Advection(elMat.A[0], fe, U, 1.0);

    // loop over test functions (i) and basis functions (j)
    for (size_t i = 1; i <= fe.N.size(); ++i) {
      for (size_t j = 1; j <= fe.N.size(); ++j) {
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

        if (stab == GLS) {
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

        if (stab == MS) {
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
    }
  }

  // Integrate source, if defined
  if (source)
    elMat.b.front().add(fe.N,f*fe.detJxW);

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
  elMat.b.front().add(fe.N,T*fe.detJxW);

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


bool AdvectionDiffusion::evalSol(Vector& s, const VecFunc& asol,
                                 const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd);
  return true;
}


std::string AdvectionDiffusion::getField1Name (size_t, const char* prefix) const
{
  if (!prefix)
    return "theta";

  return prefix + std::string(" theta");
}


std::string AdvectionDiffusion::getField2Name (size_t i,
                                               const char* prefix) const
{
  if (i > 2)
    return "";

  static const char* s[3] = { "theta_x", "theta_y", "theta_z" };
  if (!prefix)
    return s[i];

  return prefix + std::string(" ") + s[i];
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
  if (asol)
    return new AdvectionDiffusionNorm(*const_cast<AdvectionDiffusion*>(this),
                                      asol->getScalarSol(),
                                      asol->getScalarSecSol());
  else
    return new AdvectionDiffusionNorm(*const_cast<AdvectionDiffusion*>(this));
}


AdvectionDiffusionNorm::AdvectionDiffusionNorm (AdvectionDiffusion& p,
                                                RealFunc* val, VecFunc* grad) :
  NormBase(p), phi(val), gradPhi(grad)
{
  nrcmp = myProblem.getNoFields(2);
  finalOp = ASM::ABS;
}


size_t AdvectionDiffusionNorm::getNoFields (int fld) const
{
  switch (fld) {
  case 0 : return 1 + prjsol.size();
  case 1 : return gradPhi ? 5 : 2;
  default: return gradPhi ? 4 : 2;
  }
}


/*!
  \brief Returns the L2 norm contribution in an integration point.
*/

static inline double L2Norm (double val)
{
  return val*val;
}

/*!
  \brief Returns the H1 semi-norm contribution in an integration point.
*/

static inline double H1Norm (const Vec3& grad)
{
  return grad*grad;
}

/*!
  \brief Returns the residual in an integration point.
*/

static inline double residualNorm (double val, const Vec3& U,
                                   const Vec3& grad, const Vec3& hess,
                                   double kappa, double f, double react)
{
  double res = -kappa*hess.sum() + U*grad + react*val - f;
  return res*res;
}


bool AdvectionDiffusionNorm::evalInt (LocalIntegral& elmInt,
                                      const FiniteElement& fe,
                                      const Vec3& X) const
{
  AdvectionDiffusion& problem = static_cast<AdvectionDiffusion&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the advection field
  Vec3 U;
  if (problem.Uad)
    U = (*problem.Uad)(X);

  double react = problem.reaction ? (*problem.reaction)(X) : 0.0;
  double f = problem.source ? (*problem.source)(X) : 0.0;
  double h = getElementSize(fe.XC,problem.nsd);
  double val = pnorm.vec.front().dot(fe.N);

  Vec3 grad;
  Vec3 hess;
  for (size_t k = 1; k <= problem.nsd; k++)
    for (size_t j = 1; j <= fe.N.size(); j++) {
      grad[k-1] += fe.dNdX(j,k)*pnorm.vec.front()(j);
      if (this->getIntegrandType() & Integrand::SECOND_DERIVATIVES)
        hess[k-1] += fe.d2NdX2(j,k,k)*pnorm.vec.front()(j);
    }

  // Intrinsic parameter to get L^2 norm error (upper bound)
  double kk = std::min((h/(sqrt(2)*sqrt(13))),(h*h/(3*sqrt(10)*problem.props.getDiffusivity())));

  int norm = 0;

  // Choice to have VMS based error
  pnorm[norm++] += kk*kk*residualNorm(val, U, grad, hess,
                                      problem.props.getDiffusivity(), f, react)*fe.detJxW;

  // or Gradient based error depends on what to select here
  //pnorm[norm++] += H1Norm(grad)*fe.detJxW;

  pnorm[norm++] += H1Norm(grad)*fe.detJxW;

  Vec3 eGrad;
  if (phi && gradPhi) {
    double eVal = (*phi)(X);
    eGrad = (*gradPhi)(X);

    // Recovery based error
    //pnorm[norm++] += H1Norm(grad-eGrad)*fe.detJxW;

    // VMS based error estimation
    pnorm[norm++] += L2Norm(val-eVal)*fe.detJxW;

    pnorm[norm++] += kk*kk*residualNorm(val, U, grad, hess,
                                        problem.props.getDiffusivity(), f, react)*fe.detJxW;

    norm++; // effectivity index
  }

  for (size_t i = 0; i < pnorm.psol.size(); i++) {
    Vec3 rGrad;
    for (size_t k = 0; k < problem.nsd; k++)
      rGrad[k] += pnorm.psol[i].dot(fe.N,k,problem.nsd);

    pnorm[norm++] += H1Norm(rGrad)*fe.detJxW;
    // Recovery based estimate
    pnorm[norm++] += H1Norm(rGrad-grad)*fe.detJxW;

    if (phi && gradPhi) {
      pnorm[norm++] += H1Norm(eGrad-rGrad)*fe.detJxW;
      norm++; // effectivity index
    }
  }

  return true;
}


bool AdvectionDiffusionNorm::finalizeElement (LocalIntegral& elmInt)
{
  ElmNorm& norm = static_cast<ElmNorm&>(elmInt);

  // no effectiviy indices without an analytic solution
  if (!phi)
    return true;

  norm[4] = norm[3]/norm[2];

  int j = 2;
  size_t skip = this->getNoFields(1);
  for (size_t i = skip; i < norm.size(); i += this->getNoFields(j++))
    norm[i+3] = norm[i+1]/norm[2];

  return true;
}


const char* AdvectionDiffusionNorm::getName (size_t i, size_t j,
                                             const char* prefix) const
{
  static const char* s[5] = {
    "VMS explicit estimate (residual)",
    "Gradient u_h",
    "Exact error in L^2 norm,  e=u-u^h",
    "VMS explicit residual estimate",
    "eff index, residual"
  };

  static const char* p[4] = {
    "|u^r|_H1",
    "|e|_H1 e=u^r-u^h",
    "|e|_H1, e=u-u^r",
    "eff index, H1"
  };

  const char** n = i > 2 ? p : s;

  if (!prefix)
    return n[j-1];

  static std::string name;
  name = prefix + std::string(" ");
  name += n[j-1];

  return name.c_str();
}


AdvectionDiffusion::WeakDirichlet::WeakDirichlet (unsigned short int n,
                                                  double CBI_, double gamma_) :
  CBI(CBI_), gamma(gamma_), Uad(NULL), nsd(n)
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

  double g=0.0;
  if (flux)
    g = (*flux)(X);

  double An = U*normal;

  // element size
  double h = getElementSize(fe.XC,nsd);
  double C = CBI*fabs(props.getDiffusivity())/h;
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
