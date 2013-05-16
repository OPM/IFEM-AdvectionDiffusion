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
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "VTF.h"


AdvectionDiffusion::AdvectionDiffusion (unsigned short int n,
                                        AdvectionDiffusion::Stabilization s) : 
  nsd(n), kappa(1.0), Pr(1.0), stab(s), Cinv(5)
{
  primsol.resize(1);

  source = 0, Uad = 0, flux = 0, reaction = 0;
}


AdvectionDiffusion::~AdvectionDiffusion()
{
  delete source, delete Uad, delete reaction;
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
  ElementInfo* result = new ElementInfo;

  if (neumann) {
    result->withLHS = false;
    result->resize(0,1);
    result->b[0].resize(nen,true);
  } else {
    result->withLHS = true;
    result->resize(1,1);
    result->A[0].resize(nen,nen,true);
    result->b[0].resize(nen,true);
  }
  
  if (stab != NONE) {
    result->eBl = 0;
    result->eBu = 0;
    result->eMs.resize(nen,nen,true);
    result->Cv.resize(nsd+1,true);
    result->eSs.resize(nen,true);
  }

  return result;
}

// define quadratic bubble functions
static Vector bubblederivative(const FiniteElement& E, const Vec3& X)
{
  Vector value(5);
  double v[2], dv[2], sv[2];

  for (size_t i = 0; i < 2; i++)
    if ((X[i] > E.XC[0][i]) && (X[i] < E.XC[3][i])) {

      v[i] = (X[i]-E.XC[0][i])*(E.XC[3][i]-X[i])/
             ((E.XC[3][i]-E.XC[0][i])*(E.XC[3][i]-E.XC[0][i]));

      dv[i] = (E.XC[0][i]+E.XC[3][i]-2*X[i])/
              ((E.XC[3][i]-E.XC[0][i])*(E.XC[3][i]-E.XC[0][i]));

      sv[i] = -2/((E.XC[3][i]-E.XC[0][i])*(E.XC[3][i]-E.XC[0][i]));
    }
    else
      v[i] = dv[i] = sv[i] = 0.0;

  value[0] = v[0]*v[1];
  value[1] = dv[0]*v[1];
  value[2] = v[0]*dv[1];
  value[3] = sv[0]*v[0];
  value[4] = v[0]*sv[1];

  return value;
}


bool AdvectionDiffusion::evalInt (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const Vec3& X) const
{
  ElementInfo& elMat = static_cast<ElementInfo&>(elmInt);
  if (stab != NONE)
    elMat.hk = getElementSize(fe.XC, nsd);
  elMat.iEl = fe.iel;

  Vector Bub;
  if (stab != NONE && nsd == 2)
    Bub = bubblederivative(fe, X);

  double f = 0.f;
  if (source)
    f = (*source)(X);

  if (!elMat.A.empty()) {
    // evaluate reaction field
    double react = 0.f;
    if (reaction)
      react = (*reaction)(X);

    // evaluate advection field
    Vec3 U;
    if (Uad)
      U = (*Uad)(X);

    // loop over test functions (i) and basis functions (j)
    for (size_t i = 1; i <= fe.N.size(); ++i) {
      for (size_t j = 1; j <= fe.N.size(); ++j) {
        double laplace = 0.0, advect = 0.0;
        for (size_t k = 1;k <= nsd; ++k) {
          laplace += fe.dNdX(i,k)*fe.dNdX(j,k);
          advect += U[k-1]*fe.dNdX(j,k);
        }
        advect *= fe.N(i);

        elMat.A[0](i,j) += (kappa*laplace+advect+react*fe.N(i)*fe.N(j))*fe.detJxW;

        if (stab != NONE && nsd == 2) {
          elMat.eBl += (kappa*(Bub(2)*Bub(2)+Bub(3)*Bub(3))+
                       (U[0]*Bub(2)+U[1]*Bub(3))*Bub(1))*fe.detJxW;
          elMat.eBu += -(Bub(1)*Bub(1))*fe.detJxW;
        }

        if (stab == SUPG) {
          double convV=0, Lu=0;
          for (size_t k = 1;k <= nsd;k++) {
            convV += U[k-1]*fe.dNdX(i,k);
            Lu += U[k-1]*fe.dNdX(j,k)-kappa*fe.d2NdX2(j,k,k);
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
            Lv += U[k-1]*fe.dNdX(i,k)-kappa*fe.d2NdX2(i,k,k);
            Lu += U[k-1]*fe.dNdX(j,k)-kappa*fe.d2NdX2(j,k,k);
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
            Lav += -U[k-1]*fe.dNdX(i,k)-kappa*fe.d2NdX2(i,k,k);
            Lu += U[k-1]*fe.dNdX(j,k)-kappa*fe.d2NdX2(j,k,k);
            elMat.Cv(k) += U[k-1]*fe.detJxW;
          }
          Lav += react*fe.N(i);
          Lu +=  react*fe.N(j);
          elMat.Cv(nsd+1) += fe.detJxW;
          elMat.eMs(i,j) -= Lav*Lu*fe.detJxW;
          
          if (source && j == 1)
            elMat.eSs(i) -= Lav*f*fe.detJxW;
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


bool AdvectionDiffusion::evalSol(Vector& s, const Vector& N,
                                 const Matrix& dNdX,
                                 const Vec3& X,
                                 const std::vector<int>& MNPC) const
{
  Vector ePhi;
  if (utl::gather(MNPC,1,primsol.front(),ePhi))
    return false;

  dNdX.multiply(ePhi,s,true);

  return true;
}


bool AdvectionDiffusion::evalSol(Vector& s, const VecFunc& asol,
                                 const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd);
  return true;
}


const char* AdvectionDiffusion::getField1Name (size_t i,
                                               const char* prefix) const
{
  if (i >= 1)
    i = 0;

  static const char* s[1] = { "theta"};
  if (!prefix)
    return s[i];

  static std::string name;
  name = prefix + std::string(" ") + s[i];

  return name.c_str();
}


const char* AdvectionDiffusion::getField2Name (size_t i,
                                               const char* prefix) const
{
  if (i > 2)
    return 0;

  static const char* s[3] = { "theta_x", "theta_y", "theta_z" };
  if (!prefix)
    return s[i];

  static std::string name;
  name = prefix + std::string(" ") + s[i];

  return name.c_str();
}


bool AdvectionDiffusion::finalizeElement (LocalIntegral& A,
                                          const TimeDomain&, size_t)
{
  if (stab != NONE) {
    ElementInfo& E = static_cast<ElementInfo&>(A);

    // Compute stabilization parameter
    double tau = E.getTau(kappa, Cinv);
    
    E.eMs *=  tau;
    E.eSs *=  tau;
   
    tauE(E.iEl) = tau;
    
    // Add stabilization terms
    E.A[0] += E.eMs;
    E.b[0] += E.eSs;
  }

  return true;
}


double AdvectionDiffusion::ElementInfo::getTau(double kappa, double Cinv) const
{
  // eqns (211)-(214)
  double mk = std::min(1.0/3,2.0/Cinv);

  // find mean element advection velocity
  double vel=0.f;
  for (size_t k=1;k<=Cv.size()-1;++k)
    vel += pow(Cv(k)/Cv(Cv.size()),2.0);

  vel = sqrt(vel);

  double Pek = vel*hk/(2*kappa);
  double Xi = std::min(mk*Pek,1.0);

  return hk/(2*vel)*Xi;
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


double AdvectionDiffusion::getElementSize(const std::vector<Vec3>& XC, int nsd)
{
  double h;
  if (nsd == 2) {
    // Maximum diagonal
    h = (XC[3]-XC[0]).length();
    h = std::max(h, (XC[1]-XC[2]).length());
  } else {
    // Compute the element diagonals
    Vec3 D1(XC[7]-XC[0]);
    Vec3 D2(XC[6]-XC[1]);
    Vec3 D3(XC[4]-XC[3]);
    Vec3 D4(XC[2]-XC[6]);

    // and take our size as the longest one
    h = std::max(D1.length(), D2.length());
    h = std::max(h, D3.length());
    h = std::max(h, D4.length());
  }

  return h;
}


class ADElmNorm : public LocalIntegral {
  public:
    ADElmNorm(LocalIntegral& n)
    {
      myNorm = static_cast<ElmNorm*>(&n);
      iEl = 0;
    }

    ADElmNorm(size_t n)
    {
      myNorm = new ElmNorm(n);
    }

    virtual ~ADElmNorm()
    {
      myNorm->destruct();
    }

    virtual const LocalIntegral* ref() const { return myNorm; }
    LocalIntegral* ref() { return myNorm; }

    int iEl;

    ElmNorm* myNorm;
};


AdvectionDiffusionNorm::AdvectionDiffusionNorm (AdvectionDiffusion& p,
                                                RealFunc* val,
                                                VecFunc* grad) :
  NormBase(p), phi(val), gradPhi(grad)
{
  nrcmp = myProblem.getNoFields(2);
  finalOp = ASM::ABS;
}


size_t AdvectionDiffusionNorm::getNoFields (int fld) const
{
  if (fld == 0)
    return 1+prjsol.size();

  if (fld == 1)
    return gradPhi?5:2;

  return gradPhi?4:2;
}


LocalIntegral* AdvectionDiffusionNorm::getLocalIntegral (size_t nen, size_t iEl,
                                                         bool neumann) const
{
  LocalIntegral* result = NormBase::getLocalIntegral(nen,iEl,neumann);
  if (result) {
    ADElmNorm* norm = new ADElmNorm(static_cast<ElmNorm&>(*result));
    norm->iEl = iEl;
    return norm;
  }
  
  // Element norms are not requested, so allocate an inter object instead that
  // will delete itself when invoking the destruct method.
  size_t nNorms=0;
  for (size_t j=0;j<getNoFields(0);++j)
    nNorms += getNoFields(1+j);

  ADElmNorm* r2 = new ADElmNorm(nNorms);
  r2->iEl = iEl;
  return r2;
}


/*!
 \brief Returns the energy norm contribution in an integration point
*/

static inline double energyNorm(double val, const Vec3& U,
                                const Vec3& grad, double kappa,
                                double react)
{
  return kappa*grad*grad+(U*grad+react)*val;
}


/*!
 \brief Returns the H1 semi-norm contribution in an integration point
*/

static inline double H1Norm(const Vec3& grad)
{
  return grad*grad;
}

/*!
 \brief Returns the residual in an integration point
*/

static inline double residualNorm(double val, const Vec3& U,
                                  const Vec3& grad, const Vec3& hess,
                                  double kappa, double f, double react)
{
  double res = -kappa*hess.sum()+U*grad+react*val-f;
  return res*res;
}


bool AdvectionDiffusionNorm::initElement (const std::vector<int>& MNPC,
                                          const Vec3& Xc, size_t nPt,
                                          LocalIntegral& elmInt)
{
  ElmNorm& norm = *static_cast<ElmNorm*>(static_cast<ADElmNorm&>(elmInt).ref());
  return initProjection(MNPC,norm) &&
         NormBase::initElement(MNPC,Xc,nPt,norm);
}


bool AdvectionDiffusionNorm::evalInt (LocalIntegral& elmInt,
                                      const FiniteElement& fe,
                                      const Vec3& X) const
{
  AdvectionDiffusion& problem = static_cast<AdvectionDiffusion&>(myProblem);
  ElmNorm& pnorm = *(static_cast<ElmNorm*>(static_cast<ADElmNorm&>(elmInt).ref()));

  // evaluate the advection field
  Vec3 U;
  if (problem.Uad)
    U = (*problem.Uad)(X);
  double react=0;
  if (problem.reaction)
    react = (*problem.reaction)(X);
  double f=0;
  if (problem.source)
    f = (*problem.source)(X);

  double h=0;
  if (!fe.XC.empty())
    h = problem.getElementSize(fe.XC, problem.nsd);

  int norm=0;

  double val = pnorm.vec.front().dot(fe.N);
  Vec3 grad;
  Vec3 hess;
  for (size_t k=1; k <= problem.nsd; ++k) {
    for (size_t j=1; j <= fe.N.size(); ++j) {
      grad[k-1] += fe.dNdX(j,k)*pnorm.vec.front()(j);
      if (getIntegrandType() & Integrand::SECOND_DERIVATIVES)
        hess[k-1] += fe.d2NdX2(j,k,k)*pnorm.vec.front()(j);
    }
  }

  pnorm[norm++] += h*h*residualNorm(val, U, grad, hess,
                                    problem.kappa, f, react)*fe.detJxW;
  pnorm[norm++] += H1Norm(grad)*fe.detJxW;

  double eVal;
  Vec3 eGrad;
  if (phi && gradPhi) {
    eVal = (*phi)(X);
    eGrad = (*gradPhi)(X);

    pnorm[norm++] += H1Norm(grad-eGrad)*fe.detJxW;
    pnorm[norm++] += energyNorm(val-eVal, U, grad-eGrad,
                                problem.kappa, react)*fe.detJxW;
    norm++; // effectivity index
  }

  for (size_t i=0;i<pnorm.psol.size();++i) {
    Vec3 rGrad;
    for (size_t k=0;k<problem.nsd;++k)
      rGrad[k] += pnorm.psol[i].dot(fe.N,k,problem.nsd);

    pnorm[norm++] += H1Norm(rGrad)*fe.detJxW;

    pnorm[norm++] += H1Norm(grad-rGrad)*fe.detJxW;

    if (phi && gradPhi) {
      pnorm[norm++] += H1Norm(eGrad-rGrad)*fe.detJxW;
      norm++; // effectivity index
    }
  }

  return true;
}


bool AdvectionDiffusionNorm::finalizeElement (LocalIntegral& elmInt,
                                              const TimeDomain&, size_t)
{
  ADElmNorm& pnorm = static_cast<ADElmNorm&>(elmInt);
//  AdvectionDiffusion& problem = static_cast<AdvectionDiffusion&>(myProblem);
  ElmNorm& norm = const_cast<ElmNorm&>(static_cast<const ElmNorm&>((*pnorm.ref())));
  // norm[1] *= fabs(problem.getElementTau(pnorm.iEl));

  // no effectiviy indices without an analytic solution
  if (!phi)
    return true;

  size_t skip = getNoFields(1);
  norm[4] = norm[0]/norm[1];
  int j=2;
  for (size_t i = skip; i < norm.size(); i += getNoFields(j++))
    norm[i+3] = norm[i+1]/norm[2];

  return true;
}


const char* AdvectionDiffusionNorm::getName (size_t i, size_t j,
                                             const char* prefix) const
{
  static const char* s[5] = {
    "R(u^h)",
    "Gradient u_h",
    "|e|_H1, e=u-u^h",
    "a(e,e)^0.5, e=u-u^h",
    "eff index, residual"
  };

  static const char* p[4] = {
    "|u^r|_H1",
    "|e|_H1 e=u^r-u^h",
    "|e|_H1, e'=u-u^r",
    "eff index, H1"
  };

  const char** n = i > 1 ? p : s;

  if (!prefix)
    return n[j-1];

  static std::string name;
  name = prefix + std::string(" ");
  name += n[j-1];

  return name.c_str();
}

AdvectionDiffusion::WeakDirichlet::WeakDirichlet (unsigned short int n,
                                                  double CBI_,
                                                  double gamma_) :
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

  result->withLHS = true;
  result->resize(1,1);
  result->A[0].resize(nen,nen);
  result->b[0].resize(nen);

  return result;
}


bool AdvectionDiffusion::WeakDirichlet::initElementBou(const std::vector<int>& MNPC,
                                                       LocalIntegral& A)
{
  size_t nvec   = primsol.size();

  int ierr = 0;
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
  double h = AdvectionDiffusion::getElementSize(fe.XC, nsd);
  double C = CBI*fabs(kappa)/h;

  // loop over test functions (i) and basis functions (j)
  for (size_t i = 1; i <= fe.N.size(); ++i) {
    double addI = 0.0;
    for (int k=1; k <= nsd; ++k)
      addI += fe.dNdX(i,k)*normal[k-1];
    addI *= gamma*kappa;
    for (size_t j = 1; j <= fe.N.size(); ++j) {
      // adjoint
      double addJ = 0.0;
      for (int k=1; k <= nsd; ++k)
        addJ += fe.dNdX(j,k)*normal[k-1];
      elMat.A[0](i,j) += (-kappa*addJ+An*fe.N(j))*fe.N(i)*fe.detJxW;

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
