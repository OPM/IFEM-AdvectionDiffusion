// $Id$
//==============================================================================
//!
//! \file AdvectionDiffusionArgs.C
//!
//! \date Mar 9 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preparsing of input files for the Advection-Diffusion application.
//!
//==============================================================================

#include "AdvectionDiffusionArgs.h"
#include "Utilities.h"
#include "tinyxml.h"


bool AdvectionDiffusionArgs::parseArg (const char* argv)
{
  TimeIntegration::Method tmp;
  if (argv[0] != '-')
    return false;
  else if ((tmp = TimeIntegration::get(argv+1)) > TimeIntegration::NONE)
    timeMethod = tmp;
  else if (!strcmp(argv,"-supg"))
    integrandType = Integrand::SECOND_DERIVATIVES | Integrand::G_MATRIX;
  else
    return this->SIMargsBase::parseArg(argv);

  return true;
}


bool AdvectionDiffusionArgs::parse (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"timestepping")) {
    std::string type;
    utl::getAttribute(elem,"tol",errTol);
    if (utl::getAttribute(elem,"type",type))
      timeMethod = TimeIntegration::get(type);
  }

  return this->SIMargsBase::parse(elem);
}
