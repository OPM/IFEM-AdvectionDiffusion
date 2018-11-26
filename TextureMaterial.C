// $Id$
//==============================================================================
//!
//! \file TextureMaterial.C
//!
//! \date Sep 13 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Isotropic linear elastic material model defined through a texture.
//!
//==============================================================================

#include "TextureMaterial.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


void TextureMaterial::parse (const TiXmlElement* elem)
{
  std::string textureFile;
  utl::getAttribute(elem, "file", textureFile);

  int width, height, nrChannels;
  unsigned char* image = stbi_load(textureFile.c_str(),
                                   &width, &height, &nrChannels, 0);
  if (!image) {
    std::cerr << "File not found: " << textureFile << std::endl;
    return;
  }

  textureData.resize(width,std::vector<rgba>(height,{0.0,0.0,0.0}));
  const unsigned char* data = image;
  for (int j = height-1; j >= 0; j--)
    for (int i = 0; i < width; i++)
      for (int c = 0; c < nrChannels; c++)
        textureData[i][j][c] = double(*data++) / 255.0;

  free(image);

  Doubles     range;
  AD::FluidProperties mat;
  const TiXmlElement* child = elem->FirstChildElement("range");
  for (; child; child = child->NextSiblingElement("range"))
  {
    utl::getAttribute(child,"min",range.first);
    utl::getAttribute(child,"max",range.second);
    IFEM::cout << (materials.empty() ? "\n\t" : "\t");
    mat.parse(child);
    materials[range] = mat;
  }
}


void TextureMaterial::printLog () const
{
  for (const std::pair<Doubles,AD::FluidProperties>& mat : materials)
  {
    IFEM::cout <<"Material with range ["
               << mat.first.first <<","<< mat.first.second <<"]:\n";
    mat.second.printLog();
  }
}


const AD::FluidProperties* TextureMaterial::findMaterial (const FiniteElement& fe) const
{
  if (textureData.empty())
    return nullptr;

  int nrow = textureData.size();
  int ncol = textureData.front().size();
  int i = fe.u * (nrow-1);
  int j = fe.v * (ncol-1);
  if (i < 0 || i >= nrow || j < 0 || j >= ncol)
  {
    std::cerr <<" *** Texture index out of bounds "<< i <<" "<< j << std::endl;
    return nullptr;
  }

  double I = textureData[i][j].front();
  auto mat = std::find_if(materials.begin(),materials.end(),
                          [I](const std::pair<Doubles,AD::FluidProperties>& a)
                          {
                            return a.first.first <= I && I <= a.first.second;
                          });

  return mat == materials.end() ? nullptr : &mat->second;
}

double TextureMaterial::getDiffusionConstant(const FiniteElement& fe) const {
  const AD::FluidProperties* mat = this->findMaterial(fe);
  return mat ? mat->getDiffusionConstant(fe) : 1.0;
}

double TextureMaterial::getDiffusivity(const FiniteElement& fe) const {
  const AD::FluidProperties* mat = this->findMaterial(fe);
  return mat ? mat->getDiffusivity(fe) : 1.0;
}

double TextureMaterial::getReactionConstant(const FiniteElement& fe) const {
  const AD::FluidProperties* mat = this->findMaterial(fe);
  return mat ? mat->getReactionConstant(fe) : 1.0;
}

