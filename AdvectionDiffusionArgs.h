// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionArgs.h
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preparsing of input files for the Advection-Diffusion application.
//!
//==============================================================================

#include "AppCommon.h"
#include "Integrand.h"
#include "TimeIntUtils.h"

/*! \brief Struct holding application parameters.
 */

class AdvectionDiffusionArgs : public SIM::AppXMLInputBase
{
public:
  bool adap = false; //!< True to run an adaptive simulator
  double errTol = 0.0; //!< Error tolerance for embedded time stepping
  TimeIntegration::Method timeMethod = TimeIntegration::NONE; //!< Time integration method
  int integrandType = Integrand::STANDARD; //!< Integrand formulation

protected:
  //! \brief Parse an element from the input file
  bool parse(const TiXmlElement* elem) override;
};
