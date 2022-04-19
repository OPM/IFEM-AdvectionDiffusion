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

#include "ADFluidProperties.h"

#include "ElmMats.h"
#include "EqualOrderOperators.h"
#include "Functions.h"
#include "IntegrandBase.h"
#include "MatVec.h"
#include "SIMenums.h"
#include "TimeIntUtils.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>


class AnaSol;
class Fields;
class FiniteElement;
class LocalIntegral;
class Vec3;
class VecFunc;


/*!
  \brief Class representing the integrand of the Advection-Diffusion problem.
*/

class AdvectionDiffusion : public IntegrandBase
{
public:
  using WeakOps     = EqualOrderOperators::Weak;     //!< Convenience renaming
  using ResidualOps = EqualOrderOperators::Residual; //!< Convenience renaming

  //! \brief Enum defining the available stabilization methods.
  enum Stabilization { NONE, SUPG, GLS, MS };

  //! \brief Class representing the Robin boundary conditions.
  class Robin : public IntegrandBase
  {
  public:
    //! \brief Default constructor.
    //! \param[in] n Number of spatial dimensions
    //! \param[in] itg Main integrand instance
    Robin(unsigned short int n, const AdvectionDiffusion& itg);
    //! \brief Empty destructor.
    virtual ~Robin() {}

    //! \brief Returns that this integrand has no interior contributions.
    bool hasInteriorTerms() const override { return false; }

    using IntegrandBase::getLocalIntegral;
    //! \brief Returns a local integral contribution object for given element.
    //! \param[in] nen Number of nodes on element
    //! \param[in] iEl Element number
    LocalIntegral* getLocalIntegral(size_t nen, size_t iEl, bool) const override
    {
      return integrand.getLocalIntegral(nen, iEl, false);
    }

    using IntegrandBase::evalBou;
    //! \brief Evaluates the integrand at a boundary point.
    //! \param elmInt The local integral object to receive the contributions
    //! \param[in] fe Finite element data of current integration point
    //! \param[in] X Cartesian coordinates of current integration point
    //! \param[in] normal Boundary normal vector at current integration point
    bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                 const Vec3& X, const Vec3& normal) const override;

    //! \brief Set coefficient and constant.
    //! \param f The coefficient function
    void setAlpha(const VecFunc* f) { alpha = f; g = nullptr; }

    //! \brief Set flux.
    //! \param f The flux function
    void setFlux(const RealFunc* f) { alpha = nullptr; g = f; }

  protected:
    const VecFunc* alpha = nullptr; //!< Coefficient - alpha(1) * u + du/dn = alpha(2)
    const RealFunc* g = nullptr; //!< Coefficient - u + du/dn = g
    const AdvectionDiffusion& integrand; //!< Main integrand instance
  };

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] s Stabilization option
  AdvectionDiffusion(unsigned short int n = 3, Stabilization s = NONE);

  //! \brief Class representing the Advection-Diffusion element matrices.
  class ElementInfo : public ElmMats
  {
  public:
    //! \brief Default constructor.
    explicit ElementInfo(bool lhs = true) : ElmMats(lhs), hk(0.0), iEl(0) {}
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
  void setSource(RealFunc* src) { source.reset(src); }

  //! \brief Defines the Cinv stabilization parameter.
  void setCinv(double Cinv_) { Cinv = Cinv_; }
  //! \brief Returns the current Cinv value.
  double getCinv() const { return Cinv; }

  //! \brief Defines the Cbar element size parameter.
  void setCbar(double Cbar_) { Cbar = Cbar_; }

  //! \brief Returns the current Cbar value.
  double getCbar() const { return Cbar; }

  //! \brief Defines whether or not to use modified element size in residual estimator.
  void setModified(bool use) { useModified = use; }

  //! \brief Returns whether or not to use modified element size in residual estimator.
  bool useModifiedElmSize() const { return useModified; }

  //! \brief Defines the stabilization type.
  void setStabilization(Stabilization s) { stab = s; }

  //! \brief Obtain the current stabilization type.
  Stabilization getStabilization() const { return stab; }

  //! \brief Defines the advection field.
  void setAdvectionField(VecFunc* U) { Uad.reset(U); }

  //! \brief Defines the flux function.
  void setFlux(RealFunc* f) { flux = f; }

  //! \brief Defines the reaction field.
  void setReactionField(RealFunc* f) { reaction.reset(f); }

  //! \brief Defines the global number of elements.
  void setElements(size_t els) { tauE.resize(els); }

  //! \brief Sets the basis order.
  void setOrder(int p) { order = p; }

  //! \brief Returns a previously calculated tau value for the given element.
  //! \brief param[in] e The element number
  //! \details Used with norm calculations
  double getElementTau(size_t e) const { return e > tauE.size() ? 0 : tauE(e); }

  //! \brief Defines which FE quantities are needed by the integrand.
  int getIntegrandType() const override;

  //! \brief Set whether we evaluate residual norm.
  void setResidualNorm(bool on) { residualNorm = on; }

  //! \brief True if we evaluate residual norm.
  bool doResidualNorm() const { return residualNorm; }

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                  bool neumann) const override;

  using IntegrandBase::finalizeElement;
  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents.
  bool finalizeElement(LocalIntegral& elmInt) override;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  using IntegrandBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] eV Element-level primary solution vectors
  //! \param[in] fe Finite element data at current point
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3&) const override;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld) const override { return fld > 1 ? nsd : 1; }
  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  std::string getField1Name(size_t, const char* prefix) const override;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField2Name(size_t i, const char* prefix) const override;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol) const override;

  //! \brief Advances the integrand one time step forward.
  virtual void advanceStep() {}

  //! \brief Returns a reference to the fluid properties.
  AD::FluidProperties& getFluidProperties() { return props; }
  //! \brief Returns a const reference to the fluid properties.
  const AD::FluidProperties& getFluidProperties() const { return props; }

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  void setMode(SIM::SolutionMode mode) override;

  //! \brief Set time scaling factor from ODE solver.
  //! \param scale Time scaling factor
  void setTimeScale(double scale) { timeScale = scale; }

  //! \brief Set advection formulation.
  //! \param form Formulation to use.
  void setAdvectionForm(WeakOperators::ConvectionForm form)
  { advForm = form; };

  //! \brief Returns used advection form.
  WeakOperators::ConvectionForm getAdvForm() const
  { return advForm; }

  //! \brief Returns time integration method used.
  virtual TimeIntegration::Method getTimeMethod() const
  { return TimeIntegration::NONE; }

  //! \brief Registers where we can inject a mixed-basis vector field.
  void setNamedFields(const std::string&, Fields*) override;


protected:
  std::unique_ptr<VecFunc>  Uad;      //!< Pointer to advection field
  std::unique_ptr<RealFunc> reaction; //!< Pointer to the reaction field
  std::unique_ptr<RealFunc> source;   //!< Pointer to source field
  RealFunc* flux;     //!< Pointer to the flux field
  double timeScale = 1.0; //!< Time scale factor

  Vector tauE;  //!< Stored tau values - need for norm integration
  int    order; //!< Basis order

  AD::FluidProperties props; //!< Fluid properties

  Stabilization stab; //!< The type of stabilization used
  double        Cinv; //!< Stabilization parameter
  double        Cbar = 0.0; //!< Used in element size calculations
  bool  residualNorm; //!< If \e true, we will evaluate residual norm
  bool useModified = false; //!< If \e true use modified element size in residual estimate
  WeakOperators::ConvectionForm advForm = WeakOperators::CONVECTIVE; //!< Advection formulation to use

  Vectors velocity; //!< The advecting velocity field
  std::array<std::unique_ptr<Fields>,2> uFields; //!< Externally provided velocity fields

  friend class AdvectionDiffusionNorm;
};


/*!
  \brief Class representing the integrand of Advection-Diffusion energy norms.
*/

class AdvectionDiffusionNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The heat equation problem to evaluate norms for
  //! \param[in] a The analytical aolution (optional)
  explicit AdvectionDiffusionNorm(AdvectionDiffusion& p, AnaSol* a = nullptr);
  //! \brief Empty destructor.
  virtual ~AdvectionDiffusionNorm() {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  size_t getNoFields(int group = 0) const override;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  std::string getName(size_t i, size_t j, const char* prefix) const override;

  using NormBase::finalizeElement;
  //! \brief Finalizes the element norms after the numerical integration.
  //! \details This method is used to compute effectivity indices.
  //! \param elmInt The local integral object to receive the contributions
  bool finalizeElement(LocalIntegral&) override;

  //! \brief Defines which FE quantities are needed by the integrand.
  int getIntegrandType() const override;

private:
  AnaSol* anasol; //!< Analytical solution
};

#endif
