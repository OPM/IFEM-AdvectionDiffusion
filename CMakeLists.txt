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

if(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

add_library(CommonAD STATIC AdvectionDiffusion.C
                            AdvectionDiffusionBDF.C
                            AdvectionDiffusionExplicit.C
                            FluidProperties.C)

add_executable(AdvectionDiffusion main_AdvectionDiffusion.C)

target_link_libraries(AdvectionDiffusion CommonAD IFEMAppCommon ${IFEM_LIBRARIES})

# Installation
install(TARGETS AdvectionDiffusion DESTINATION bin)

# For generating the doxy
if(IFEM_AS_SUBMODULE OR IFEM_INTREE_BUILD)
  set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}")
  configure_file(${IFEM_PATH}/doc/Doxyfile.in Doxyfile)
else()
  configure_file(doc/Doxyfile.in Doxyfile)
endif()
configure_file(doc/AdvectionDiffusion.dox.in AdvectionDiffusion.dox)
add_custom_target(doc doxygen ${PROJECT_BINARY_DIR}/Doxyfile
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  COMMENT "Generating API documentation" VERBATIM)

if(IFEM_INSTALL_DOXY)
  install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")" COMPONENT doc)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html DESTINATION ${CMAKE_INSTALL_DOCDIR}/Apps/AdvectionDiffusion
          COMPONENT doc
          PATTERN *.md5 EXCLUDE
          PATTERN *.map EXCLUDE)
endif(IFEM_INSTALL_DOXY)

# Testing
enable_testing()
include(IFEMTesting)

# Generate regtest script with correct paths
configure_file(${IFEM_REGTEST_SCRIPT} regtest.sh)

if(MPI_FOUND)
  # Add parallel tests here
else()
  file(GLOB AD_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                         ${PROJECT_SOURCE_DIR}/Test/*.reg)
  foreach(TESTFILE ${AD_TESTFILES})
    ifem_add_test(${TESTFILE} AdvectionDiffusion)
  endforeach()
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
