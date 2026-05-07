cmake_minimum_required(VERSION 3.22)

project(ADLib LANGUAGES C CXX)

ifem_add_library(
  NAME
    ADR
  ALIAS
    AppLibrary::AdvectionDiffusion
  SOURCES
    AdvectionDiffusion.C
    AdvectionDiffusionBDF.C
    AdvectionDiffusionExplicit.C
    AdvectionDiffusionImplicit.C
    AdvectionDiffusionSource.C
    ADFluidProperties.C
    SIMAD.C
  HEADERS
    AdvectionDiffusion.h
    AdvectionDiffusionBDF.h
    AdvectionDiffusionExplicit.h
    AdvectionDiffusionImplicit.h
    AdvectionDiffusionSource.h
    ADFluidProperties.h
    SIMAD.h
  LIBRARIES
    AppLibrary::Common
    IFEM
)

# Unit tests
ifem_add_test_app(
  SOURCES
    ${PROJECT_SOURCE_DIR}/Test/TestAdvectionDiffusionSource.C
    ${PROJECT_SOURCE_DIR}/Test/TestFluidProperties.C
    ${PROJECT_SOURCE_DIR}/Test/TestSIMAD.C
    WORKDIR
      ${PROJECT_SOURCE_DIR}/../../Test
    NAME
      AdvectionDiffusion
    LIBRARIES
      AppLibrary::AdvectionDiffusion
)

set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
