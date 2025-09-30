//==============================================================================
//!
//! \file TestFluidProperties.C
//!
//! \date Oct 8 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for fluid properties for advection diffusion problems.
//!
//==============================================================================

#include "ADFluidProperties.h"
#include <tinyxml2.h>

#include "Catch2Support.h"


TEST_CASE("TestADFluidProperties.Physical")
{
  tinyxml2::XMLDocument doc;
  doc.Parse("<fluidproperties rho=\"2.0\" kappa=\"3.0\" C=\"4.0\"/>");
  AD::FluidProperties props;
  props.parse(doc.RootElement());

  REQUIRE_THAT(props.getMassDensity(), WithinRel(2.0));
  REQUIRE_THAT(props.getDiffusivity(), WithinRel(3.0));
  REQUIRE_THAT(props.getHeatCapacity(), WithinRel(4.0));
}


TEST_CASE("TestADFluidProperties.PrRa")
{
  tinyxml2::XMLDocument doc;
  doc.Parse("<fluidproperties type=\"Pr/Ra\" Pr=\"2.0\" Ra=\"3.0\"/>");
  AD::FluidProperties props;
  props.parse(doc.RootElement());

  REQUIRE_THAT(props.getPrandtlNumber(), WithinRel(2.0));
  REQUIRE_THAT(props.getRayleighNumber(), WithinRel(3.0));
}
