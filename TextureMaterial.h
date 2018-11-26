// $Id$
//==============================================================================
//!
//! \file TextureMaterial.h
//!
//! \date Sep 13 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Isotropic linear elastic material model defined through a texture.
//!
//==============================================================================

#ifndef _ISOTROPIC_TEXTURE_MAT_H
#define _ISOTROPIC_TEXTURE_MAT_H

#include "ADFluidProperties.h"
#include <array>
#include <vector>
#include <map>


/*!
  \brief Class representing an isotropic material model with a texture.
*/

class TextureMaterial : public AD::FluidProperties
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param[in] ps If \e true, assume plane stress in 2D
  //! \param[in] ax If \e true, assume 3D axi-symmetric material
  TextureMaterial() : AD::FluidProperties() {}
  //! \brief Empty destructor.
  virtual ~TextureMaterial() = default;

  //! \brief Parses material parameters from an XML element.
  void parse(const TiXmlElement* elem) override;

  //! \brief Prints out material parameters to the log stream.
  void printLog() const override;

  //! \brief Returns the diffusion constant acccording to scaling.
  double getDiffusionConstant(const FiniteElement& fe) const override;

  //! \brief Returns the reaction constant acccording to scaling.
  double getReactionConstant(const FiniteElement& fe) const override;

  //! \brief Returns thermal diffusivity.
  double getDiffusivity(const FiniteElement& fe) const override;

private:
  //! \brief Locates the appropriate material as indicated by texture.
  const AD::FluidProperties* findMaterial(const FiniteElement& fe) const;

protected:
  typedef std::pair<double,double> Doubles; //!< Convenience type
  typedef std::array<double,4>     rgba;    //!< RGB color code

  //! Material for different texture regions
  std::map< Doubles,AD::FluidProperties > materials;
  //! Raw image texture information describing the spatial material variation
  std::vector< std::vector<rgba> > textureData;
};

#endif
