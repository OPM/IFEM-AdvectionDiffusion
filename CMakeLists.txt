PROJECT(AdvectionDiffusion)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

# Required packages
SET(IFEM_PRECONFIGURED ${IFEM_CONFIGURED})
IF (NOT IFEM_CONFIGURED)
  FIND_PACKAGE(IFEM REQUIRED)
ENDIF(NOT IFEM_CONFIGURED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

INCLUDE_DIRECTORIES(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR})

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

IF(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ENDIF(NOT WIN32)

IF (NOT IFEM_PRECONFIGURED)
  ADD_SUBDIRECTORY(../Common CommonIFEM.dir)
ENDIF(NOT IFEM_PRECONFIGURED)

add_library(CommonAD STATIC AdvectionDiffusion.C
                            AdvectionDiffusionBDF.C
                            AdvectionDiffusionExplicit.C
                            FluidProperties.C)

ADD_EXECUTABLE(AdvectionDiffusion main_AdvectionDiffusion.C)

TARGET_LINK_LIBRARIES(AdvectionDiffusion CommonAD ${IFEM_LIBRARIES} CommonIFEM)

# Installation
INSTALL(TARGETS AdvectionDiffusion DESTINATION bin)

# Testing
ENABLE_TESTING()

# Generate regtest script with correct paths
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/../Common/scripts/regtest.sh.in
               regtest.sh)

IF(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
ELSE(IFEM_USE_PARALLEL_PETSC)
  FILE(GLOB AD_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                         ${PROJECT_SOURCE_DIR}/Test/*.reg)
  FOREACH(TESTFILE ${AD_TESTFILES})
    IFEM_ADD_TEST(${TESTFILE} AdvectionDiffusion)
  ENDFOREACH()
ENDIF(IFEM_USE_PARALLEL_PETSC)
list(APPEND TEST_APPS AdvectionDiffusion)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  AdvectionDiffusion
                  CommonAD ${IFEM_LIBRARIES} CommonIFEM)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
