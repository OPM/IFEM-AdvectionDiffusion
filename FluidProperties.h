// $Id$
//==============================================================================
//!
//! \file FluidProperties.h
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

class TiXmlElement;


/*!
  \brief Class representing fluid properties for advection-diffusion problems.
*/

namespace AD {

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

  //! \brief Empty destructor.
  ~FluidProperties() {}

  //! \brief Parses material parementers from an XML element.
  void parse(const TiXmlElement*);

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

  //! \brief Return mass and advection constant according to scaling.
  double getMassAdvectionConstant() const;

  //! \brief Returns the diffusion constant acccording to scaling.
  double getDiffusionConstant() const;

  //! \brief Returns the reaction constant acccording to scaling.
  double getReactionConstant() const;
protected:
  // Physical properties (constant)
  double rho;   //!< Mass density.
  double kappa; //!< Thermal diffusivity.
  double C;     //!< Heat capacity.

  double Ra; //!< Rayleigh number.
  double Pr; //!< Prandtl number.

  Scaling scaling; //!< Scaling to use.
};

}

#endif
