//==============================================================================
//!
//! \file TestSIMAD.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for solution driver for Advection-Diffusion problems.
//!
//==============================================================================

#include "AdvectionDiffusionBDF.h"
#include "SIMAD.h"
#include "SIM2D.h"

#include "gtest/gtest.h"

TEST(TestSIMAD, Parse)
{
  SIMAD<SIM2D> sim(new AdvectionDiffusionBDF(2, 2, 0), true);
  EXPECT_TRUE(sim.read("Lshape.xinp"));
}
