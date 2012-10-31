// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusion.h
//!
//! \date Jun 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Advection-Diffusion problems.
//!
//==============================================================================

#ifndef _ADVECTIONDIFFUSION_H
#define _ADVECTIONDIFFUSION_H

#include "IntegrandBase.h"
#include "ElmMats.h"
#include "Vec3.h"


/*!
  \brief Class representing the integrand of the Advection-Diffusion problem.
  \details The blah. 
*/

class AdvectionDiffusion : public IntegrandBase
{
public:
  enum Stabilization { NONE, SUPG, GLS, MS };

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  AdvectionDiffusion(unsigned short int n = 3,
                     Stabilization stab=NONE);

  class ElementInfo : public ElmMats
  {
    public:
      //! \brief Default constructor
      ElementInfo() {}
      
      //! \brief Empty destructor
      virtual ~ElementInfo() {}

      //! \brief Return the stabilization parameter
      double getTau(double kappa, double Cinv) const;

      Matrix eMs; //!< Stabilized matrix
      Vector eSs; //!< Stabilized vector
      Vector Cv;  //!< velocity + area
      double eBl; //!< a(b_1^e, b_2^e)-bilinear term using quadratic bubble as trial and test functions
      double eBu; //!< b_1^e*\int_K b_2^e d\Omega -where new stabilization can be defined as \tau=eBu/eBl (assuming element residual is constant)
    
      double hk;  //!< element size
      size_t iEl; //!< element index
  };

  //! \brief Empty destructor.
  virtual ~AdvectionDiffusion();

  //! \brief Defines the source function
  void setSource(RealFunc* src) { source = src; }

  //! \brief Set the Cinv stabilization parameter
  void setCinv(double Cinv_) { Cinv = Cinv_; }

  //! \brief Set kappa
  void setKappa(double kappa_) { kappa = kappa_; }

  //! \brief Set the stabilization type
  void setStabilization(Stabilization s) { stab = s; }

  //! \brief Defines the advection field 
  void setAdvectionField(VecFunc* U) { Uad = U; }

  //! \brief Defines the flux function
  void setFlux(RealFunc* f) { flux = f; }

  //! \brief Defines the reaction field
  void setReactionField(RealFunc* f) { reaction = f; }

  //! \brief Set the global number of elements
  void setElements(size_t els)
  {
    tauE.resize(els);
  }

  void setPrandtlNumber(double Pr_)
  {
    Pr = Pr_;
  }

  //! \brief Returns a previously calculated tau value for the given element.
  //! \brief param[in] el The element number
  //!\details Used with norm calculations
  double getElementTau(size_t el)
  {
    if (el > tauE.size())
      return 0;

    return tauE(el);
  }

  virtual void advanceStep() {}

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents.
  //! It can also be used to implement multiple integration point loops within
  //! the same element, provided the necessary integration point values are
  //! stored internally in the object during the first integration loop.
  virtual bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t = 0);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite Element quantities
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt,  const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the analytical secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (tensor field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s,
                       const VecFunc& asol, const Vec3& X) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 1) const { return fld > 1 ? nsd : 1; }
  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  virtual const char* getField1Name(size_t, const char* prefix = 0) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField2Name(size_t i, const char* prefix = 0) const;

  virtual int getIntegrandType() const
  { 
    if (stab == NONE)
      return Integrand::STANDARD;

    return Integrand::SECOND_DERIVATIVES |
           Integrand::ELEMENT_CORNERS;
  }

  //! \brief Returns characteristic element size
  //! \param XC The element corner coordinates
  //! \details The size is taken as the longest diagonal 
  double getElementSize(const std::vector<Vec3>& XC) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const;

  VecFunc* Uad; //!< Pointer to advection field
  RealFunc* reaction; //<!< Pointer to the reaction field
  RealFunc* source; //!< Pointer to source
  unsigned short int nsd; //!< Number of space dimensions (1, 2 or, 3)
  double kappa; //!< diffusion coefficient
  double Pr;     //!< Prandtl number
protected:
  RealFunc* flux; //<!< Pointer to the flux field
  Vector tauE;     //<!< Stored tau values - need for norm integration

  Stabilization stab; //!< The type of stabilization used

  double Cinv; //!< stabilization parameter
};


/*!
  \brief Class representing the integrand of Advection-Diffusion energy norms.
*/

class AdvectionDiffusionNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Advection-Diffusion problem to evaluate norms for
  AdvectionDiffusionNorm(AdvectionDiffusion& p, RealFunc* u=0,
                                                VecFunc*  du=0);
  //! \brief Empty destructor.
  virtual ~AdvectionDiffusionNorm() {}

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return false; }

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] Xc Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const Vec3& Xc, size_t nPt,
                           LocalIntegral& elmInt);

  //! \brief Returns the number of norm quantities.
  //! \param[in] fld If 1, the number of exact norms of the exact solution,
  //                 else the number of norms for FE fields
  virtual size_t getNoFields(int fld=0) const;

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X) const;


  //! \brief Finalizes the element norms after the numerical integration.
  //! \details This method is used to compute effectivity indices.
  //! \param elmInt The local integral object to receive the contributions
  virtual bool finalizeElement(LocalIntegral& elmInt, const TimeDomain&,size_t);

  virtual int getIntegrandType() const
  {
    return myProblem.getIntegrandType();
  }

  const char* getName(size_t i, size_t j, const char* prefix);

  bool hasElementContributions(size_t i, size_t j)
  {
    return true;
  }
protected:
  RealFunc* phi;
  VecFunc*  gradPhi;
};

#endif
