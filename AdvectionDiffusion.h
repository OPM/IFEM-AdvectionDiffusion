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

#ifndef _ADVECTION_DIFFUSION_H
#define _ADVECTION_DIFFUSION_H

#include "IntegrandBase.h"
#include "ElmMats.h"


/*!
  \brief Class representing the integrand of the Advection-Diffusion problem.
*/

class AdvectionDiffusion : public IntegrandBase
{
public:
  //! \brief Enum defining the available stabilization methods.
  enum Stabilization { NONE, SUPG, GLS, MS };

  //! \brief Class representing the weak Dirichlet integrand.
  class WeakDirichlet : public IntegrandBase
  {
  public:
    //! \brief Default constructor.
    //! \param[in] n Number of spatial dimensions
    //! \param[in] CBI_ Model constant
    //! \param[in] gamma_ Adjoint factor
    WeakDirichlet(unsigned short int n, double CBI_ = 4, double gamma_ = 1.0);
    //! \brief The destructor deletes the advection field.
    virtual ~WeakDirichlet();

    //! \brief Returns that this integrand has no interior contributions.
    virtual bool hasInteriorTerms() const { return false; }

    //! \brief Defines which FE quantities are needed by the integrand.
    virtual int getIntegrandType() const { return ELEMENT_CORNERS; }

    //! \brief Returns a local integral contribution object for given element.
    //! \param[in] nen Number of nodes on element
    virtual LocalIntegral* getLocalIntegral(size_t nen, size_t, bool) const;

    //! \brief Initializes current element for boundary integration.
    //! \param[in] MNPC Matrix of nodal point correspondance for current element
    //! \param elmInt Local integral for element
    virtual bool initElementBou(const std::vector<int>& MNPC,
                                LocalIntegral& elmInt);

    //! \brief Evaluates the integrand at a boundary point.
    //! \param elmInt The local integral object to receive the contributions
    //! \param[in] fe Finite element data of current integration point
    //! \param[in] X Cartesian coordinates of current integration point
    //! \param[in] normal Boundary normal vector at current integration point
    virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                         const Vec3& X, const Vec3& normal) const;

    //! \brief Defines the advection field.
    void setAdvectionField(VecFunc* U) { Uad = U; }

    //! \brief Defines the flux function.
    void setFlux(RealFunc* f) { flux = f; }

    //! \brief Defines kappa.
    void setKappa(double kappa_) { kappa = kappa_; }
    //! \brief Defines the kappa function.
    void setKappa(ScalarFunc* kappa_) { kFunc = kappa_; }

  protected:
    const double CBI;   //!< Model constant
    const double gamma; //!< Adjoint factor
    VecFunc*     Uad;   //!< Pointer to advection field
    RealFunc*    flux;  //!< Pointer to the flux field
    unsigned short int nsd; //!< Number of space dimensions (1, 2 or, 3)
    double       kappa; //!< Diffusion coefficient
    ScalarFunc*  kFunc; //!< Diffusion coefficient function
  };

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] stab Stabilization option
  AdvectionDiffusion(unsigned short int n = 3, Stabilization stab = NONE);

  //! \brief Class representing the Advection-Diffusion element matrices.
  class ElementInfo : public ElmMats
  {
  public:
    //! \brief Default constructor.
    ElementInfo(bool lhs = true) : ElmMats(lhs) {}
    //! \brief Empty destructor.
    virtual ~ElementInfo() {}

    //! \brief Returns the stabilization parameter.
    double getTau(double kappa, double Cinv, int p) const;

    Matrix eMs; //!< Stabilized matrix
    Vector eSs; //!< Stabilized vector
    Vector Cv;  //!< velocity + area

    double hk;  //!< element size
    size_t iEl; //!< element index
  };

  //! \brief The destructor deletes the advection field.
  virtual ~AdvectionDiffusion();

  //! \brief Defines the source function.
  void setSource(RealFunc* src) { source = src; }

  //! \brief Defines the Cinv stabilization parameter.
  void setCinv(double Cinv_) { Cinv = Cinv_; }
  //! \brief Returns the current Cinv value.
  double getCinv() const { return Cinv; }

  //! \brief Defines kappa.
  void setKappa(double kappa_) { kappa = kappa_; }
  //! \brief Defines the kappa function.
  void setKappa(ScalarFunc* kappa_) { kFunc = kappa_; }
  //! \brief Returns current kappa.
  double getKappa() const { return kappa; }

  //! \brief Defines the stabilization type.
  void setStabilization(Stabilization s) { stab = s; }

  //! \brief Obtain the current stabilization type.
  Stabilization getStabilization() const { return stab; }

  //! \brief Defines the advection field.
  void setAdvectionField(VecFunc* U) { Uad = U; }

  //! \brief Defines the flux function.
  void setFlux(RealFunc* f) { flux = f; }

  //! \brief Defines the reaction field
  void setReactionField(RealFunc* f) { reaction = f; }

  //! \brief Defines the global number of elements.
  void setElements(size_t els) { tauE.resize(els); }

  //! \brief Sets the basis order.
  void setOrder(int p) { order = p; }

  //! \brief Defines the Prandtl number.
  void setPrandtlNumber(double Pr_) { Pr = Pr_; }

  //! \brief Returns the current Prandtl number.
  double getPrandtlNumber() const { return Pr; }

  //! \brief Returns a previously calculated tau value for the given element.
  //! \brief param[in] e The element number
  //! \details Used with norm calculations
  double getElementTau(size_t e) const { return e > tauE.size() ? 0 : tauE(e); }

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const;

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  using IntegrandBase::finalizeElement;
  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents.
  virtual bool finalizeElement(LocalIntegral& elmInt);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe,
                       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the analytical secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] asol The analytical solution field (tensor field)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const VecFunc& asol, const Vec3& X) const;

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

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const;

  //! \brief Advances the integrand one time step forward.
  virtual void advanceStep() {}

protected:
  VecFunc*  Uad;      //!< Pointer to advection field
  RealFunc* reaction; //!< Pointer to the reaction field
  RealFunc* source;   //!< Pointer to source field
  RealFunc* flux;     //!< Pointer to the flux field
  unsigned short int nsd; //!< Number of space dimensions (1, 2 or, 3)
  double kappa; //!< diffusion coefficient
  ScalarFunc* kFunc;  //!< Diffusion coefficient function
  double Pr;    //!< Prandtl number
  Vector tauE;  //!< Stored tau values - need for norm integration
  int order;    //!< Basis order

  Stabilization stab; //!< The type of stabilization used
  double        Cinv; //!< Stabilization parameter

  friend class AdvectionDiffusionNorm;
};


/*!
  \brief Class representing the integrand of Advection-Diffusion energy norms.
*/

class AdvectionDiffusionNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Advection-Diffusion problem to evaluate norms for
  //! \param[in] u Analytical solution field
  //! \param[in] du Analytical gradient field
  AdvectionDiffusionNorm(AdvectionDiffusion& p, RealFunc* u=0, VecFunc* du=0);
  //! \brief Empty destructor.
  virtual ~AdvectionDiffusionNorm() {}

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using NormBase::finalizeElement;
  //! \brief Finalizes the element norms after the numerical integration.
  //! \details This method is used to compute effectivity indices.
  virtual bool finalizeElement(LocalIntegral& elmInt);

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group = 0) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  virtual const char* getName(size_t i, size_t j, const char* prefix) const;

protected:
  RealFunc* phi;     //!< Analytical solution field
  VecFunc*  gradPhi; //!< Analytical gradient field
};

#endif
