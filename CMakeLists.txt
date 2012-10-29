PROJECT(AdvectionDiffusion)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Required defines
SET(CMAKE_CXX_FLAGS "-DReal=double")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DINDEX_CHECK=2")
IF(VERBOSE_DEBUG GREATER 0)
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DINT_DEBUG=${VERBOSE_DEBUG}")
ENDIF(VERBOSE_DEBUG GREATER 0)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

ENABLE_TESTING()

IF(NOT IFEM_BUILD_TYPE)
  SET(IFEM_BUILD_TYPE ${CMAKE_BUILD_TYPE})
ENDIF(NOT IFEM_BUILD_TYPE)

IF(${CMAKE_BUILD_TYPE} MATCHES "Nopt")
  SET(CMAKE_BUILD_TYPE Debug)
ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Nomp")
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(${CMAKE_BUILD_TYPE} MATCHES "Nopt")

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

# Required packages
FIND_PACKAGE(IFEM REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

INCLUDE_DIRECTORIES(${IFEM_INCLUDES} ../Stokes/Common/Helpers 
                    ${PROJECT_SOURCE_DIR}/../LinearElasticity
                    ../Stokes/Common/NavierStokes)

SET(EXECUTABLE_OUTPUT_PATH bin)

IF(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ENDIF(NOT WIN32)

ADD_SUBDIRECTORY(../Stokes/Common/NavierStokes CommonNavierStokes.dir)

ADD_EXECUTABLE(AdvectionDiffusion AdvectionDiffusion.C 
                                  AdvectionDiffusionBDF.C
                                  main_AdvectionDiffusion.C)
TARGET_LINK_LIBRARIES(AdvectionDiffusion ${IFEM_LIBRARIES} CommonNavierStokes)

ADD_CUSTOM_TARGET(doc doxygen 
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  COMMENT "Generating API documentation" VERBATIM)

ADD_CUSTOM_TARGET(check-commits scripts/check-patchseries.sh
                  "${IFEM_BUILD_TYPE}" "Apps/AdvectionDiffusion" 
                                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../../
                                COMMENT "Checking commits" VERBATIM)

FILE(GLOB_RECURSE AD_TESTFILES "${PROJECT_SOURCE_DIR}/Test/*.reg")
FOREACH(TESTFILE ${AD_TESTFILES})
  ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/../../scripts/regtest.sh "${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/AdvectionDiffusion" "${TESTFILE}")
ENDFOREACH()

IF(NOT IFEM_INSTALL_PREFIX)
  SET(IFEM_INSTALL_PREFIX $ENV{HOME})
ENDIF(NOT IFEM_INSTALL_PREFIX)

INSTALL(TARGETS AdvectionDiffusion DESTINATION ${IFEM_INSTALL_PREFIX}/bin)
