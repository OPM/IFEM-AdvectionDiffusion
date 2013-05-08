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

INCLUDE_DIRECTORIES(${IFEM_INCLUDES} ../Common)

SET(EXECUTABLE_OUTPUT_PATH bin)

IF(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ENDIF(NOT WIN32)

ADD_SUBDIRECTORY(../Common CommonIFEM.dir)

ADD_EXECUTABLE(AdvectionDiffusion AdvectionDiffusion.C 
                                  AdvectionDiffusionBDF.C
                                  AdvectionDiffusionExplicit.C
                                  main_AdvectionDiffusion.C)
TARGET_LINK_LIBRARIES(AdvectionDiffusion ${IFEM_LIBRARIES} CommonIFEM)

ADD_CUSTOM_TARGET(check-commits scripts/check-patchseries.sh
                  "${IFEM_BUILD_TYPE}" "Apps/AdvectionDiffusion" 
                                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../../
                                COMMENT "Checking commits" VERBATIM)

IF(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
ELSE(IFEM_USE_PARALLEL_PETSC)
  FILE(GLOB AD_TESTFILES "${PROJECT_SOURCE_DIR}/Test/*.reg")
  FOREACH(TESTFILE ${AD_TESTFILES})
    ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/../../scripts/regtest.sh "${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/AdvectionDiffusion" "${TESTFILE}")
  ENDFOREACH()
ENDIF(IFEM_USE_PARALLEL_PETSC)

INSTALL(TARGETS AdvectionDiffusion DESTINATION bin)
