project(CommonAD)

cmake_minimum_required(VERSION 2.6)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

include_directories(../../../Common)

set(AD_SOURCES AdvectionDiffusion.C
               AdvectionDiffusionBDF.C
               AdvectionDiffusionExplicit.C
               AdvectionDiffusionImplicit.C
               AdvectionDiffusionSource.C
               ADFluidProperties.C
               AdvectionDiffusion.h
               AdvectionDiffusionBDF.h
               AdvectionDiffusionExplicit.h
               AdvectionDiffusionImplicit.h
               AdvectionDiffusionSource.h
               ADFluidProperties.h
               SIMAD.C
               SIMAD.h)

add_library(CommonAD STATIC ${AD_SOURCES})

list(APPEND CHECK_SOURCES ${AD_SOURCES})

set(CHECK_SOURCES ${CHECK_SOURCES} PARENT_SCOPE)
