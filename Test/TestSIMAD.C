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
  AdvectionDiffusionBDF integrand(2, 2, 0);
  SIMAD<SIM2D,AdvectionDiffusionBDF> sim(integrand, true);
  EXPECT_TRUE(sim.read("Lshape.xinp"));

  const AdvectionDiffusionBDF& ad =
        static_cast<const AdvectionDiffusionBDF&>(*sim.getProblem());

  ASSERT_FLOAT_EQ(ad.getCinv(), 1.0);
  ASSERT_FLOAT_EQ(ad.getFluidProperties().getDiffusivity(), 1e-6);
  ASSERT_FLOAT_EQ(ad.getFluidProperties().getPrandtlNumber(), 0.5);
  EXPECT_EQ(ad.getStabilization(), AdvectionDiffusion::MS);
}
