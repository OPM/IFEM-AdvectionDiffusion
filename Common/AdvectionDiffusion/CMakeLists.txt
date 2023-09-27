cmake_minimum_required(VERSION 3.5)

project(CommonAD)

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

foreach(src ${AD_SOURCES})
  list(APPEND CHECK_SOURCES ${PROJECT_SOURCE_DIR}/${src})
endforeach()

# Unit tests
set(TEST_SOURCES ${PROJECT_SOURCE_DIR}/Test/TestAdvectionDiffusionSource.C
                 ${PROJECT_SOURCE_DIR}/Test/TestFluidProperties.C
                 ${PROJECT_SOURCE_DIR}/Test/TestSIMAD.C)
IFEM_add_test_app("${TEST_SOURCES}"
                  ${PROJECT_SOURCE_DIR}/../../Test
                  AdvectionDiffusion 0
                  CommonAD IFEMAppCommon ${IFEM_LIBRARIES})

set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
set(CHECK_SOURCES ${CHECK_SOURCES} PARENT_SCOPE)
