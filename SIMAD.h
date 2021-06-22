// $Id$
//==============================================================================
//!
//! \file SIMAD.h
//!
//! \date Jun 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for Advection-Diffusion problems.
//!
//==============================================================================

#ifndef _SIM_AD_H
#define _SIM_AD_H

#include "SIMMultiPatchModelGen.h"
#include "SIMsolution.h"
#include "SIMconfigure.h"
#include "SIMenums.h"
#ifdef HAS_PETSC
#include "PETScSupport.h"
#endif


class DataExporter;
class SIMoutput;
class TimeStep;


/*!
  \brief Driver class for an Advection-diffusion simulator.
  \details The class incapsulates data and methods for solving a
  Advection-Diffusion problem using NURBS-based finite elements.
*/

template<class Dim, class Integrand>
class SIMAD : public SIMMultiPatchModelGen<Dim>, public SIMsolution
{
public:
  //! \brief Setup properties.
  struct SetupProps
  {
    bool shareGrid = false; //!< True to share grid with some other simulator
    Integrand* integrand = nullptr; //!< Integrand to use
    SIMoutput* share = nullptr; //!< Simulator to share grid with
    bool standalone = false; //!< Simulator runs standalone
  };

  //! \brief Default constructor.
  //! \param[in] ad Integrand for advection-diffusion problem
  //! \param[in] alone Integrand is used stand-alone (controls time stepping)
  SIMAD(Integrand& ad, bool alone = false);

  //! \brief Constructs from given properties.
  explicit SIMAD(const SetupProps& props) :
    SIMAD(*props.integrand, props.standalone)
  {}

  //! \brief The destructor clears the VTF-file pointer, unless stand-alone.
  //! \details This is needed when the VTF-file is assumed to be owned by
  //! another SIM-object and is deleted by the SIMbase destructor of that one.
  virtual ~SIMAD();

  //! \brief Initializes the property containers of the model.
  void clearProperties() override;

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem) override;

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "AdvectionDiffusion"; }

  //! \brief Initialize simulator.
  bool init(const TimeStep&);

  //! \brief Preprocessing performed before the FEM model generation.
  //! \details This method is reimplemented to couple the weak Dirichlet
  //! integrand to the Robin property codes.
  void preprocessA() override;

  //! \brief Defines the global number of elements.
  bool preprocessB() override;

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  bool saveModel(char* fileName, int& geoBlk, int&);

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep&);

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp, bool = false);

  //! \brief No solution postprocessing.
  bool postSolve(const TimeStep&) { return true; }

  //! \brief Evaluates and prints out solution norms.
  void printFinalNorms(const TimeStep& tp);

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock);

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(SerializeMap& data) const override;

  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const SerializeMap& data) override;

  //! \brief Returns a reference to current solution vector.
  Vector& getSolution() { return solution.front(); }
  //! \brief Returns a const reference to current solution vector.
  //! \details If set, returns external vector for \a idx=0.
  const Vector& getSolution(int idx) const override
  {
    if (idx == 0 && extsol)
      return *extsol;

    return solution[idx];
  }

  //! \brief Register fields for simulation result export.
  void registerFields(DataExporter& exporter, const std::string& prefix="");

  //! \brief Set context to read from input file
  void setContext(int ctx);

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  //!
  //! \details Since this solver is linear, this is just a normal solve.
  SIM::ConvStatus solveIteration(TimeStep& tp);

#ifdef HAS_PETSC
  //! \brief Sets MPI communicator for the linear equation solvers.
  void setCommunicator(const MPI_Comm* comm) { Dim::adm.setCommunicator(comm); }
#endif

  //! \brief Sets the externally provided solution vector (adaptive simulation).
  void setSol(const Vector* sol) { extsol = sol; }

  //! \brief Prints a summary of the calculated solution to std::cout.
  //! \param[in] solution The solution vector
  //! \param[in] printSol Print solution only if size is less than this value
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  void printSolutionSummary(const Vector& solution, int printSol, const char*,
                            std::streamsize outPrec = 0) override
  {
    this->SIMMultiPatchModelGen<Dim>::printSolutionSummary(solution, printSol,
                                                           "temperature ", outPrec);
  }

  //! \brief Set time scaling factor from ODE solver.
  //! \param scale Time scaling factor
  void setTimeScale(double scale) { AD.setTimeScale(scale); }

  //! \brief Prints integrated solution norms to the log stream.
  //! \param[in] norms The norm values
  void printNorms (const Vectors& norms, size_t) const override;

  //! \brief Returns a reference to the element norm matrix.
  Matrix& getElmNorms() { return eNorm; }

  //! \brief Returns a reference to the projection vectors.
  Vectors& getProjections() { return projs; }

protected:
  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  bool initNeumann(size_t propInd) override;

  //! \brief Prints solution norms to the log stream.
  //! \param[in] gNorm The norm values
  void printExactNorms(const Vector& gNorm) const;

  //! \brief Prints a norm group to the log stream.
  void printNormGroup(const Vector& rNorm, const Vector& fNorm,
                      const std::string& name) const override;

private:
  Integrand& AD; //!< Problem integrand definition
  typename Integrand::Robin robinBC; //!< Robin integrand definition
  Vectors projs; //!< Projection vectors
  Matrix eNorm; //!< Element norms

  const Vector* extsol = nullptr; //!< Solution vector for adaptive simulators
  bool standalone = false; //!< If \e true, this simulator owns the VTF object
  std::string inputContext; //!< Input context
  int aCode[2] = {0}; //!< Analytical BC code (used by destructor)
};


//! \brief Partial specialization for configurator
template<class Dim, class Integrand>
struct SolverConfigurator< SIMAD<Dim,Integrand> > {
  //! \brief Configure a SIMAD instance.
  //! \param ad The SIMAD instance to configure
  //! \param[in] props Configuration properties
  //! \param[in] infile The input file to read
  int setup(SIMAD<Dim,Integrand>& ad,
            const typename SIMAD<Dim,Integrand>::SetupProps& props,
            char* infile);
};

#endif
