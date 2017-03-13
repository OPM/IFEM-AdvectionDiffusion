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


bool AdvectionDiffusionArgs::parse(const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"advectiondiffusion"))
    utl::getAttribute(elem,"adaptive",adap);

  if (!strcasecmp(elem->Value(),"timestepping")) {
    utl::getAttribute(elem, "tol", errTol);

    std::string type;
    utl::getAttribute(elem,"type",type);
    if (type == "be")
      timeMethod = TimeIntegration::BE;
    else if (type == "bdf2")
      timeMethod = TimeIntegration::BDF2;
    else if (type == "cn")
      timeMethod = TimeIntegration::THETA;
    else if (type == "euler")
      timeMethod = TimeIntegration::EULER;
    else if (type == "heun")
      timeMethod = TimeIntegration::HEUN;
    else if (type == "heuneuler")
      timeMethod = TimeIntegration::HEUNEULER;
    else if (type == "bs")
      timeMethod = TimeIntegration::BOGACKISHAMPINE;
    else if (type == "fehlberg")
      timeMethod = TimeIntegration::FEHLBERG;
    else if (type == "rk3")
      timeMethod = TimeIntegration::RK3;
    else if (type == "rk4")
      timeMethod = TimeIntegration::RK4;
  }

  return SIM::AppXMLInputBase::parse(elem);
}
