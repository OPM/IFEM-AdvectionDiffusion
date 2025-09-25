cmake_minimum_required(VERSION 3.22)

project(AdvectionDiffusion LANGUAGES C CXX)

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

if(NOT TARGET CommonAD)
  add_subdirectory(Common/AdvectionDiffusion)
endif()

ifem_add_application(
  NAME
    AdvectionDiffusion
  SOURCES
    main_AdvectionDiffusion.C
    AdvectionDiffusionArgs.C
  HEADERS
    AdvectionDiffusionArgs.h
  LIBRARIES
    CommonAD
)

enable_testing()

ifem_add_hdf5_test(
  TARGET
    AdvectionDiffusion
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    MPI/Square-2-LR-bdf2.hreg
  PARALLEL
    2
)

ifem_add_hdf5_test(
  TARGET
    AdvectionDiffusion
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Square-ad-adap-rec.hreg
)

ifem_add_hdf5_test(
  TARGET
    AdvectionDiffusion
  TEST_FILES
    Square-ad.hreg
)

ifem_add_vtf_test(
  TARGET
    AdvectionDiffusion
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Square-ad-adap-rec.vreg
    Square-ad-adap-varkappa.vreg
)

ifem_add_vtf_test(
  TARGET
    AdvectionDiffusion
  TEST_FILES
    Square-ad.vreg
    Square-robin.vreg
)

ifem_add_restart_test(
  TARGET
    AdvectionDiffusion
  TEST_FILES
    MPI/Square-ad-restart.reg
  RESTART_LEVEL
    5
  PARALLEL
    4
)

ifem_add_restart_test(
  TARGET
    AdvectionDiffusion
  TEST_FILES
    Square-abd2-ad-restart.reg
  RESTART_LEVEL
    5
)

ifem_add_regression_test(
  TARGET
    AdvectionDiffusion
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    MPI/Square-2-LR-bdf2.reg
  PARALLEL
    2
)

ifem_add_regression_test(
  TARGET
    AdvectionDiffusion
  CONDITIONS
    LRSpline_FOUND
    Zoltan_FOUND
  TEST_FILES
    MPI/Square-ad-adap.reg
  PARALLEL
    4
)

ifem_add_regression_test(
  TARGET
    AdvectionDiffusion
  TEST_FILES
    Cube-robin.reg
    Lshape.reg
    Square-abd1-ad-bdf2.reg
    Square-abd1-ad-be.reg
    Square-abd1-ad-bs.reg
    Square-abd1-ad-cn.reg
    Square-abd1-ad-euler.reg
    Square-abd1-ad-heuneuler.reg
    Square-abd1-ad-heun.reg
    Square-abd1-ad-rk3.reg
    Square-abd1-ad-rk4.reg
    Square-abd2-ad-ab1.reg
    Square-abd2-ad-ab2.reg
    Square-abd2-ad-ab3.reg
    Square-abd2-ad-ab4.reg
    Square-abd2-ad-ab5.reg
    Square-abd2-ad-am1.reg
    Square-abd2-ad-am2.reg
    Square-abd2-ad-am3.reg
    Square-abd2-ad-am4.reg
    Square-abd2-ad-bdf2.reg
    Square-abd2-ad-be.reg
    Square-abd2-ad-bs.reg
    Square-abd2-ad-cn.reg
    Square-abd2-ad-euler.reg
    Square-abd2-ad-heuneuler.reg
    Square-abd2-ad-heun.reg
    Square-abd2-ad-rk3.reg
    Square-abd2-ad-rk4.reg
    Square-ad.reg
    Square-ad-varkappa.reg
    Square-ad-varkappa-autodiff.reg
    Square-robin.reg
)

ifem_add_regression_test(
  TARGET
    AdvectionDiffusion
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Square-2-LR-bdf2.reg
    Square-ad-adap.reg
    Square-ad-adap-anasource.reg
    Square-ad-adap-rec.reg
    Square-ad-adap-rec-noana.reg
    Square-ad-adap-res.reg
    Square-ad-adap-varkappa.reg
    Square-ad-adap-res-noana.reg
)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(TARGET AdvectionDiffusion)

# Installation
include(GNUInstallDirs)
install(TARGETS AdvectionDiffusion DESTINATION bin)
