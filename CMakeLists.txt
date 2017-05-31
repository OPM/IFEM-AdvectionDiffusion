project(AdvectionDiffusion)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

# Required packages
IF (NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
ENDIF(NOT IFEM_CONFIGURED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

include_directories(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

set(AD_SOURCES AdvectionDiffusion.C
               AdvectionDiffusionArgs.C
               AdvectionDiffusionBDF.C
               AdvectionDiffusionExplicit.C
               ADFluidProperties.C)

add_library(CommonAD STATIC ${AD_SOURCES})

add_executable(AdvectionDiffusion main_AdvectionDiffusion.C)
list(APPEND CHECK_SOURCES ${AD_SOURCES} main_AdvectionDiffusion.C)

target_link_libraries(AdvectionDiffusion CommonAD IFEMAppCommon ${IFEM_LIBRARIES})

# Installation
install(TARGETS AdvectionDiffusion DESTINATION bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}")
add_doc_target(AdvectionDiffusion AdvectionDiffusion)

# Testing
enable_testing()
include(IFEMTesting)

if(MPI_FOUND)
  if(HDF5_FOUND AND CEREAL_FOUND)
    ifem_add_restart_test(MPI/Square-ad-restart.reg AdvectionDiffusion 5 4)
  endif()
endif()
if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  set(TESTFILES    Lshape.reg
                   Square-abd1-ad-bdf2.reg
                   Square-abd1-ad-be.reg
                   Square-abd1-ad-bs.reg
                   Square-abd1-ad-cn.reg
                   Square-abd1-ad-euler.reg
                   Square-abd1-ad-heuneuler.reg
                   Square-abd1-ad-heun.reg
                   Square-abd1-ad-rk3.reg
                   Square-abd1-ad-rk4.reg
                   Square-abd2-ad-bdf2.reg
                   Square-abd2-ad-be.reg
                   Square-abd2-ad-bs.reg
                   Square-abd2-ad-cn.reg
                   Square-abd2-ad-euler.reg
                   Square-abd2-ad-heuneuler.reg
                   Square-abd2-ad-heun.reg
                   Square-abd2-ad-rk3.reg
                   Square-abd2-ad-rk4.reg
                   Square-ad-RaPr.reg
                   Square-ad.reg)
  if(NOT __LINUX_COMPILER_PGI) # These tests are extremely compiler sensitive, disable.
    list(APPEND TESTFILES Square-abd1-ad-fehlberg.reg
                          Square-abd2-ad-fehlberg.reg)
  endif()
  foreach(TESTFILE ${TESTFILES})
    ifem_add_test(${TESTFILE} AdvectionDiffusion)
  endforeach()
  if(HDF5_FOUND AND CEREAL_FOUND)
    ifem_add_restart_test(Square-abd2-ad-restart.reg AdvectionDiffusion 5)
  endif()
endif()
list(APPEND TEST_APPS AdvectionDiffusion)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  AdvectionDiffusion
                  CommonAD IFEMAppCommon ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
