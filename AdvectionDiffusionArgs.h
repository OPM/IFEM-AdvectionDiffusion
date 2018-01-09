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

#ifndef _ADVECTION_DIFFUSION_ARGS_H
#define _ADVECTION_DIFFUSION_ARGS_H

#include "SIMargsBase.h"
#include "Integrand.h"
#include "TimeIntUtils.h"


/*!
  \brief Class holding application parameters.
*/

class AdvectionDiffusionArgs : public SIMargsBase
{
public:
  double errTol = 1e-6; //!< Error tolerance for embedded time stepping
  TimeIntegration::Method timeMethod = TimeIntegration::NONE; //!< Time integration method
  int integrandType = Integrand::STANDARD; //!< Integrand formulation

  //! \brief Default constructor.
  AdvectionDiffusionArgs() : SIMargsBase("advectiondiffusion") {}

  //! \brief Parses a command-line argument.
  virtual bool parseArg(const char* argv);

protected:
  //! \brief Parse an element from the input file
  bool parse(const TiXmlElement* elem) override;
};

#endif
