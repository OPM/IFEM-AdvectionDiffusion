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

#include "AdvectionDiffusion.h"
#include "AdvectionDiffusionBDF.h"
#include "ADFluidProperties.h"
#include "SIMAD.h"

#include "SIM2D.h"
#include "TimeIntUtils.h"

#include "Catch2Support.h"


TEST_CASE("TestSIMAD.Parse")
{
  AdvectionDiffusionBDF integrand(2, TimeIntegration::BDF2, 0);
  SIMAD<SIM2D,AdvectionDiffusionBDF> sim(integrand, true);
  REQUIRE(sim.read("Lshape.xinp"));

  const AdvectionDiffusionBDF& ad =
        static_cast<const AdvectionDiffusionBDF&>(*sim.getProblem());

  REQUIRE_THAT(ad.getCinv(), WithinRel(1.0));
  REQUIRE_THAT(ad.getFluidProperties().getDiffusivity(), WithinRel(1e-6));
  REQUIRE_THAT(ad.getFluidProperties().getPrandtlNumber(), WithinRel(0.5));
  REQUIRE(ad.getStabilization() == AdvectionDiffusion::MS);
}
