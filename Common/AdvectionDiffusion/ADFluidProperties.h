// $Id$
//==============================================================================
//!
//! \file ADFluidProperties.h
//!
//! \date Aug 11 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Advection-Diffusion fluid properties.
//!
//==============================================================================

#ifndef _AD_FLUID_PROPERTIES_H
#define _AD_FLUID_PROPERTIES_H

#include <memory>
#include <string>

class RealFunc;
class Vec3;
namespace tinyxml2 { class XMLElement; }

namespace AD {

/*!
  \brief Class representing fluid properties for advection-diffusion problems.
*/

class FluidProperties
{
public:
  //! \brief Enumeration of available scalings
  enum Scaling {
    PHYSICAL, //!< physical dimensions.
    PR_RA     //!< Prandtl/Rayleigh scaling.
  };

  //! \brief Empty constructor.
  FluidProperties();

  //! \brief Destructor.
  ~FluidProperties();

  //! \brief Parses material parameters from an XML element.
  void parse(const tinyxml2::XMLElement*);

  //! \brief Prints out fluid properties to the log stream.
  void printLog() const;

  //! \brief Returns the mass density.
  double getMassDensity() const { return rho; }

  //! \brief Returns thermal diffusivity.
  double getDiffusivity() const { return kappa; }

  //! \brief Returns heat capacity.
  double getHeatCapacity() const { return C; }

  //! \brief Returns Prandtl number.
  double getPrandtlNumber() const { return Pr; }

  //! \brief Return Rayleigh number.
  double getRayleighNumber() const { return Ra; }

  //! \brief Returns the thermal expansion coefficient.
  double getExpansion() const { return alpha; }

  //! \brief Return mass and advection constant according to scaling.
  double getMassAdvectionConstant() const;

  //! \brief Returns the diffusion constant acccording to scaling.
  double getDiffusionConstant(const Vec3& X) const;

  //! \brief Returns the reaction constant acccording to scaling.
  double getReactionConstant(const Vec3& X) const;

  //! \brief Get thermal expansion term
  //! \param[in] T Temperature
  double getThermalExpansion(double T) const;

  //! \brief Returns kappa function.
  const RealFunc* kappaFunc() const
  { return kappaF.get(); }

  //! \brief Returns a new copy of kappa function.
  RealFunc* newKappaFunc() const;

protected:
  // Physical properties (constant)
  double rho;   //!< Mass density
  double kappa; //!< Thermal diffusivity
  double alpha; //!< Thermal expansion coefficient
  double C;     //!< Heat capacity

  //! \brief Struct holding a function definition.
  struct FuncDef {
    std::string type; //!< Type of the function
    std::string def; //!< Textual definition of the function
  };

  double Ra; //!< Rayleigh number
  FuncDef RaFdef; //!< Definition of Rayleigh number as a function
  std::unique_ptr<RealFunc> RaF; //!< Rayleigh number as a function
  FuncDef kappaFdef; //!< Definition of kappa as a function
  std::unique_ptr<RealFunc> kappaF; //!< Scalar thermal diffusivity as a function
  double Pr; //!< Prandtl number

  Scaling scaling; //!< Scaling to use.
};

}

#endif
