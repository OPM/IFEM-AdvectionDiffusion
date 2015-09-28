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

#include "FluidProperties.h"
#include <tinyxml.h>

#include "gtest/gtest.h"

TEST(TestFluidProperties, ADPhysical)
{
  TiXmlDocument doc;
  doc.Parse("<fluidproperties rho=\"2.0\" kappa=\"3.0\" C=\"4.0\"/>");
  AD::FluidProperties props;
  props.parse(doc.RootElement());

  EXPECT_FLOAT_EQ(props.getMassDensity(), 2.0);
  EXPECT_FLOAT_EQ(props.getDiffusivity(), 3.0);
  EXPECT_FLOAT_EQ(props.getHeatCapacity(), 4.0);
}


TEST(TestFluidProperties, PrRa)
{
  TiXmlDocument doc;
  doc.Parse("<fluidproperties type=\"Pr/Ra\" Pr=\"2.0\" Ra=\"3.0\"/>");
  AD::FluidProperties props;
  props.parse(doc.RootElement());

  EXPECT_FLOAT_EQ(props.getPrandtlNumber(), 2.0);
  EXPECT_FLOAT_EQ(props.getRayleighNumber(), 3.0);
}
