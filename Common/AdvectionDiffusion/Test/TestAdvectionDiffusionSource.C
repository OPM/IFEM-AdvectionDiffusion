//==============================================================================
//!
//! \file TestAdvectionDiffusionSource.C
//!
//! \date Sep 28 2023
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for Advection-Diffusion source functions.
//!
//==============================================================================

#include "AdvectionDiffusionSource.h"

#include "ADFluidProperties.h"

#include "AnaSol.h"

#include "ExprFunctions.h"
#include "Functions.h"
#include "Vec3.h"
#include <tinyxml.h>

#include "gtest/gtest.h"


TEST(TestAdvectionDiffusionComponentSource, Parse2D)
{
  TiXmlDocument doc;
  doc.Parse("<source type=\"components\">"
            "<laplacian type=\"expression\">"
            "  -4*sin(2*x)*sin(y) | -sin(2*x)*sin(y)"
            "</laplacian>"
            "<temperature_grad type=\"expression\">"
            "  2*cos(2*x)*sin(y) | sin(2*x)*cos(y)"
            "</temperature_grad>"
            "<velocity type=\"expression\">"
            " x*y | sin(x)*cos(y)"
            "</velocity>"
            "</source>");

  AD::FluidProperties props;
  AD::AdvectionDiffusionSource source(doc.RootElement(),props);

  auto Tx  = [](const Vec3& X) { return 2*cos(2*X.x)*sin(X.y); };
  auto Txx = [](const Vec3& X) { return -4*sin(2*X.x)*sin(X.y); };
  auto Ty  = [](const Vec3& X) { return sin(2*X.x)*cos(X.y); };
  auto Tyy = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y); };
  auto u   = [](const Vec3& X) { return X.x*X.y; };
  auto v   = [](const Vec3& X) { return sin(X.x)*cos(X.y); };

  for (const double x : {1.0, 2.0})
    for (const double y : {3.0, 4.0}) {
      Vec3 X(x,y);
      EXPECT_DOUBLE_EQ(source(X), -1.0*(Txx(X) + Tyy(X)) + u(X)*Tx(X) + v(X)*Ty(X));
    }
}


TEST(TestAdvectionDiffusionComponentSource, Parse2DReaction)
{
  TiXmlDocument doc;
  doc.Parse("<source type=\"components\">"
            "<laplacian type=\"expression\">"
            "  -4*sin(2*x)*sin(y) | -sin(2*x)*sin(y)"
            "</laplacian>"
            "<temperature_grad type=\"expression\">"
            "  2*cos(2*x)*sin(y) | sin(2*x)*cos(y)"
            "</temperature_grad>"
            "<velocity type=\"expression\">"
            " x*y | sin(x)*cos(y)"
            "</velocity>"
            "<temperature type=\"expression\">sin(2*x)*sin(y)</temperature>"
            "</source>");

  AD::FluidProperties props;
  ConstFunc react(0.5);
  AD::AdvectionDiffusionSource source(doc.RootElement(), props, &react);

  auto T   = [](const Vec3& X) { return 0.5*sin(2*X.x)*sin(X.y); };
  auto Tx  = [](const Vec3& X) { return 2*cos(2*X.x)*sin(X.y); };
  auto Txx = [](const Vec3& X) { return -4*sin(2*X.x)*sin(X.y); };
  auto Ty  = [](const Vec3& X) { return sin(2*X.x)*cos(X.y); };
  auto Tyy = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y); };
  auto u   = [](const Vec3& X) { return X.x*X.y; };
  auto v   = [](const Vec3& X) { return sin(X.x)*cos(X.y); };

  for (const double x : {1.0, 2.0})
    for (const double y : {3.0, 4.0}) {
      Vec3 X(x,y);
      EXPECT_DOUBLE_EQ(source(X), -1.0*(Txx(X) + Tyy(X)) + u(X)*Tx(X) + v(X)*Ty(X) + T(X));
    }
}


TEST(TestAdvectionDiffusionComponentSource, Parse3D)
{
  TiXmlDocument doc;
  doc.Parse("<source type=\"components\">"
            "<laplacian type=\"expression\">"
            "  -4*sin(2*x)*sin(y)*cos(z) | -sin(2*x)*sin(y)*cos(z) | -sin(2*x)*sin(y)*cos(z)"
            "</laplacian>"
            "<temperature_grad type=\"expression\">"
            "  2*cos(2*x)*sin(y)*cos(z) | sin(2*x)*cos(y)*cos(z) | -sin(2*x)*sin(y)*sin(z)"
            "</temperature_grad>"
            "<velocity type=\"expression\">"
            " x*y*exp(z) | sin(x)*cos(y)*sin(z) | x*y*z"
            "</velocity>"
            "</source>");

  AD::FluidProperties props;
  AD::AdvectionDiffusionSource source(doc.RootElement(),props);

  auto Tx  = [](const Vec3& X) { return 2*cos(2*X.x)*sin(X.y)*cos(X.z); };
  auto Txx = [](const Vec3& X) { return -4*sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto Ty  = [](const Vec3& X) { return sin(2*X.x)*cos(X.y)*cos(X.z); };
  auto Tyy = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto Tz  = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*sin(X.z); };
  auto Tzz = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto u   = [](const Vec3& X) { return X.x*X.y*exp(X.z); };
  auto v   = [](const Vec3& X) { return sin(X.x)*cos(X.y)*sin(X.z); };
  auto w   = [](const Vec3& X) { return X.x*X.y*X.z; };

  for (const double x : {1.0, 2.0})
    for (const double y : {3.0, 4.0})
      for (const double z : {0.5, 0.6}) {
        Vec3 X(x,y,z);
        EXPECT_DOUBLE_EQ(source(X), -1.0*(Txx(X) + Tyy(X) + Tzz(X)) +
                                        u(X)*Tx(X) + v(X)*Ty(X) + w(X)*Tz(X));
      }
}


TEST(TestAdvectionDiffusionComponentSource, Parse3DReaction)
{
  TiXmlDocument doc;
  doc.Parse("<source type=\"components\">"
            "<laplacian type=\"expression\">"
            "  -4*sin(2*x)*sin(y)*cos(z) | -sin(2*x)*sin(y)*cos(z) | -sin(2*x)*sin(y)*cos(z)"
            "</laplacian>"
            "<temperature type=\"expression\">sin(2*x)*sin(y)*cos(z)</temperature>"
            "<temperature_grad type=\"expression\">"
            "  2*cos(2*x)*sin(y)*cos(z) | sin(2*x)*cos(y)*cos(z) | -sin(2*x)*sin(y)*sin(z)"
            "</temperature_grad>"
            "<velocity type=\"expression\">"
            " x*y*exp(z) | sin(x)*cos(y)*sin(z) | x*y*z"
            "</velocity>"
            "</source>");

  AD::FluidProperties props;
  EvalFunction react("x*y*z");
  AD::AdvectionDiffusionSource source(doc.RootElement(),props,&react);

  auto T   = [](const Vec3& X) { return X.x*X.y*X.z*sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto Tx  = [](const Vec3& X) { return 2*cos(2*X.x)*sin(X.y)*cos(X.z); };
  auto Txx = [](const Vec3& X) { return -4*sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto Ty  = [](const Vec3& X) { return sin(2*X.x)*cos(X.y)*cos(X.z); };
  auto Tyy = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto Tz  = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*sin(X.z); };
  auto Tzz = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto u   = [](const Vec3& X) { return X.x*X.y*exp(X.z); };
  auto v   = [](const Vec3& X) { return sin(X.x)*cos(X.y)*sin(X.z); };
  auto w   = [](const Vec3& X) { return X.x*X.y*X.z; };

  for (const double x : {1.0, 2.0})
    for (const double y : {3.0, 4.0})
      for (const double z : {0.5, 0.6}) {
        Vec3 X(x,y,z);
        EXPECT_NEAR(source(X), -1.0*(Txx(X) + Tyy(X) + Tzz(X)) +
                               u(X)*Tx(X) + v(X)*Ty(X) + w(X)*Tz(X) + T(X), 1e-12);
      }
}


TEST(TestAdvectionDiffusionAnaSolSource, Parse2D)
{
  TiXmlDocument doc;
  doc.Parse("<anasol type=\"expression\" autodiff=\"true\">"
            "<primary>sin(2*x)*sin(y)</primary>"
            "</anasol>");

  AD::FluidProperties props;
  AnaSol aSol(doc.RootElement(), true);
  aSol.setupSecondarySolutions();
  VecFuncExpr U("x*y | sin(x)*cos(y)");
  EXPECT_TRUE(aSol.getScalarSol() != nullptr);
  EXPECT_TRUE(aSol.getScalarSecSol() != nullptr);
  AD::AdvectionDiffusionAnaSolSource source(aSol,U,props,nullptr,true);

  auto Tx  = [](const Vec3& X) { return 2*cos(2*X.x)*sin(X.y); };
  auto Txx = [](const Vec3& X) { return -4*sin(2*X.x)*sin(X.y); };
  auto Ty  = [](const Vec3& X) { return sin(2*X.x)*cos(X.y); };
  auto Tyy = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y); };
  auto u   = [](const Vec3& X) { return X.x*X.y; };
  auto v   = [](const Vec3& X) { return sin(X.x)*cos(X.y); };

  for (const double x : {1.0, 2.0})
    for (const double y : {3.0, 4.0}) {
      Vec3 X(x,y);
      EXPECT_DOUBLE_EQ(source(X), -1.0*(Txx(X) + Tyy(X)) + u(X)*Tx(X) + v(X)*Ty(X));
    }
}


TEST(TestAdvectionDiffusionAnaSolSource, Parse2DReaction)
{
  TiXmlDocument doc;
  doc.Parse("<anasol type=\"expression\" autodiff=\"true\">"
            "<primary>sin(2*x)*sin(y)</primary>"
            "</anasol>");

  AD::FluidProperties props;
  AnaSol aSol(doc.RootElement(), true);
  aSol.setupSecondarySolutions();
  VecFuncExpr U("x*y | sin(x)*cos(y)");
  ConstFunc react(0.5);
  EXPECT_TRUE(aSol.getScalarSol() != nullptr);
  EXPECT_TRUE(aSol.getScalarSecSol() != nullptr);
  AD::AdvectionDiffusionAnaSolSource source(aSol,U,props,&react,true);

  auto T   = [](const Vec3& X) { return 0.5*sin(2*X.x)*sin(X.y); };
  auto Tx  = [](const Vec3& X) { return 2*cos(2*X.x)*sin(X.y); };
  auto Txx = [](const Vec3& X) { return -4*sin(2*X.x)*sin(X.y); };
  auto Ty  = [](const Vec3& X) { return sin(2*X.x)*cos(X.y); };
  auto Tyy = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y); };
  auto u   = [](const Vec3& X) { return X.x*X.y; };
  auto v   = [](const Vec3& X) { return sin(X.x)*cos(X.y); };

  for (const double x : {1.0, 2.0})
    for (const double y : {3.0, 4.0}) {
      Vec3 X(x,y);
      EXPECT_DOUBLE_EQ(source(X), -1.0*(Txx(X) + Tyy(X)) + u(X)*Tx(X) + v(X)*Ty(X) + T(X));
    }
}


TEST(TestAdvectionDiffusionAnaSolSource, Parse2DUnsteady)
{
  TiXmlDocument doc;
  doc.Parse("<anasol type=\"expression\" autodiff=\"true\">"
            "<primary>sin(2*x)*sin(y)*sin(t)</primary>"
            "</anasol>");

  AD::FluidProperties props;
  AnaSol aSol(doc.RootElement(), true);
  aSol.setupSecondarySolutions();
  VecFuncExpr U("x*y*cos(t) | sin(x)*cos(y)*exp(t)");
  EXPECT_TRUE(aSol.getScalarSol() != nullptr);
  EXPECT_TRUE(aSol.getScalarSecSol() != nullptr);
  AD::AdvectionDiffusionAnaSolSource source(aSol,U,props,nullptr,false);

  auto Tt  = [](const Vec4& X) { return sin(2*X.x)*sin(X.y)*cos(X.t); };
  auto Tx  = [](const Vec4& X) { return 2*cos(2*X.x)*sin(X.y)*sin(X.t); };
  auto Txx = [](const Vec4& X) { return -4*sin(2*X.x)*sin(X.y)*sin(X.t); };
  auto Ty  = [](const Vec4& X) { return sin(2*X.x)*cos(X.y)*sin(X.t); };
  auto Tyy = [](const Vec4& X) { return -sin(2*X.x)*sin(X.y)*sin(X.t); };
  auto u   = [](const Vec4& X) { return X.x*X.y*cos(X.t); };
  auto v   = [](const Vec4& X) { return sin(X.x)*cos(X.y)*exp(X.t); };

  for (const double t : {0.1, 0.2})
    for (const double x : {1.0, 2.0})
      for (const double y : {3.0, 4.0}) {
        Vec4 X(x,y,0,t);
          EXPECT_DOUBLE_EQ(source(X), -1.0*(Txx(X) + Tyy(X)) +
                                      u(X)*Tx(X) + v(X)*Ty(X) + Tt(X));
      }
}


TEST(TestAdvectionDiffusionAnaSolSource, Parse3D)
{
  TiXmlDocument doc;
  doc.Parse("<anasol type=\"expression\" autodiff=\"true\">"
            "<primary>sin(2*x)*sin(y)*cos(z)</primary>"
            "</anasol>");

  AD::FluidProperties props;
  AnaSol aSol(doc.RootElement(), true);
  aSol.setupSecondarySolutions();
  VecFuncExpr U("x*y*exp(z) | sin(x)*cos(y)*sin(z) | x*y*z");
  EXPECT_TRUE(aSol.getScalarSol() != nullptr);
  EXPECT_TRUE(aSol.getScalarSecSol() != nullptr);
  AD::AdvectionDiffusionAnaSolSource source(aSol,U,props,nullptr,true);

  auto Tx  = [](const Vec3& X) { return 2*cos(2*X.x)*sin(X.y)*cos(X.z); };
  auto Txx = [](const Vec3& X) { return -4*sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto Ty  = [](const Vec3& X) { return sin(2*X.x)*cos(X.y)*cos(X.z); };
  auto Tyy = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto Tz  = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*sin(X.z); };
  auto Tzz = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto u   = [](const Vec3& X) { return X.x*X.y*exp(X.z); };
  auto v   = [](const Vec3& X) { return sin(X.x)*cos(X.y)*sin(X.z); };
  auto w   = [](const Vec3& X) { return X.x*X.y*X.z; };

  for (const double x : {1.0, 2.0})
    for (const double y : {3.0, 4.0})
      for (const double z : {0.5, 0.6}) {
        Vec3 X(x,y,z);
        EXPECT_DOUBLE_EQ(source(X), -1.0*(Txx(X) + Tyy(X) + Tzz(X)) +
                                        u(X)*Tx(X) + v(X)*Ty(X) + w(X)*Tz(X));
      }
}


TEST(TestAdvectionDiffusionAnaSolSource, Parse3DReaction)
{
  TiXmlDocument doc;
  doc.Parse("<anasol type=\"expression\" autodiff=\"true\">"
            "<primary>sin(2*x)*sin(y)*cos(z)</primary>"
            "</anasol>");

  AD::FluidProperties props;
  AnaSol aSol(doc.RootElement(), true);
  aSol.setupSecondarySolutions();
  VecFuncExpr U("x*y*exp(z) | sin(x)*cos(y)*sin(z) | x*y*z");
  EvalFunction react("x*y*z");
  EXPECT_TRUE(aSol.getScalarSol() != nullptr);
  EXPECT_TRUE(aSol.getScalarSecSol() != nullptr);
  AD::AdvectionDiffusionAnaSolSource source(aSol,U,props,&react,true);

  auto T   = [](const Vec3& X) { return X.x*X.y*X.z*sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto Tx  = [](const Vec3& X) { return 2*cos(2*X.x)*sin(X.y)*cos(X.z); };
  auto Txx = [](const Vec3& X) { return -4*sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto Ty  = [](const Vec3& X) { return sin(2*X.x)*cos(X.y)*cos(X.z); };
  auto Tyy = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto Tz  = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*sin(X.z); };
  auto Tzz = [](const Vec3& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z); };
  auto u   = [](const Vec3& X) { return X.x*X.y*exp(X.z); };
  auto v   = [](const Vec3& X) { return sin(X.x)*cos(X.y)*sin(X.z); };
  auto w   = [](const Vec3& X) { return X.x*X.y*X.z; };

  for (const double x : {1.0, 2.0})
    for (const double y : {3.0, 4.0})
      for (const double z : {0.5, 0.6}) {
        Vec3 X(x,y,z);
        EXPECT_NEAR(source(X), -1.0*(Txx(X) + Tyy(X) + Tzz(X)) +
                               u(X)*Tx(X) + v(X)*Ty(X) + w(X)*Tz(X) + T(X), 1e-12);
      }
}


TEST(TestAdvectionDiffusionAnaSolSource, Parse3DUnsteady)
{
  TiXmlDocument doc;
  doc.Parse("<anasol type=\"expression\" autodiff=\"true\">"
            "<primary>sin(2*x)*sin(y)*cos(z)*cos(t)</primary>"
            "</anasol>");

  AD::FluidProperties props;
  AnaSol aSol(doc.RootElement(), true);
  aSol.setupSecondarySolutions();
  VecFuncExpr U("x*y*exp(z)*t | sin(x)*cos(y)*sin(z)*t*t | x*y*z*sin(t)");
  EXPECT_TRUE(aSol.getScalarSol() != nullptr);
  EXPECT_TRUE(aSol.getScalarSecSol() != nullptr);
  AD::AdvectionDiffusionAnaSolSource source(aSol,U,props,nullptr,false);

  auto Tt  = [](const Vec4& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z)*sin(X.t); };
  auto Tx  = [](const Vec4& X) { return 2*cos(2*X.x)*sin(X.y)*cos(X.z)*cos(X.t); };
  auto Txx = [](const Vec4& X) { return -4*sin(2*X.x)*sin(X.y)*cos(X.z)*cos(X.t); };
  auto Ty  = [](const Vec4& X) { return sin(2*X.x)*cos(X.y)*cos(X.z)*cos(X.t); };
  auto Tyy = [](const Vec4& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z)*cos(X.t); };
  auto Tz  = [](const Vec4& X) { return -sin(2*X.x)*sin(X.y)*sin(X.z)*cos(X.t); };
  auto Tzz = [](const Vec4& X) { return -sin(2*X.x)*sin(X.y)*cos(X.z)*cos(X.t); };
  auto u   = [](const Vec4& X) { return X.x*X.y*exp(X.z)*X.t; };
  auto v   = [](const Vec4& X) { return sin(X.x)*cos(X.y)*sin(X.z)*X.t*X.t; };
  auto w   = [](const Vec4& X) { return X.x*X.y*X.z*sin(X.t); };

  for (const double t : {0.4, 0.8})
    for (const double x : {1.0, 2.0})
      for (const double y : {3.0, 4.0})
        for (const double z : {0.5, 0.6}) {
          Vec4 X(x,y,z,t);
          EXPECT_DOUBLE_EQ(source(X), -1.0*(Txx(X) + Tyy(X) + Tzz(X)) +
                                          u(X)*Tx(X) + v(X)*Ty(X) + w(X)*Tz(X) + Tt(X));
        }
}
