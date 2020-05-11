//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef FACE_H
#define FACE_H

//! \file      Face.h
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.1
//! \date      June 5 2019

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include "../Maths/Coord.h"
#include "../Errors.h"

class Face;

#include "Element.h"

class Face
{
public:
  Face();
  virtual ~Face();

  //Accesseurs
  const Coord& getNormal() const { return m_normal; };
  const Coord& getTangent() const { return m_tangent; };
  const Coord& getBinormal() const { return m_binormal; };
  const double& getSurface() const { return m_surface; };
  const Coord& getPos() const { return m_position; };

  virtual void setSurface(const double &surface){ Errors::errorMessage("setSurface not available for requested face"); };
  virtual void initializeAutres(const double &surface, const Coord &normal, const Coord &tangent, const Coord &binormal){ Errors::errorMessage("initializeAutres not available for requested face"); }
  virtual void setPos(const double &X, const double &Y, const double &Z) { Errors::errorMessage("setPos not available for requested face"); };
  virtual void setNormal(const double &X, const double &Y, const double &Z) { Errors::errorMessage("setNormal not available for requested face"); };
  virtual void setTangent(const double &X, const double &Y, const double &Z) { Errors::errorMessage("setTangent not available for requested face"); };
  virtual void setBinormal(const double &X, const double &Y, const double &Z) { Errors::errorMessage("setBinormal not available for requested face"); };
  virtual void setSize(const double &sizeX, const double &sizeY, const double &sizeZ) { Errors::errorMessage("setSize not available for requested face"); };
  virtual void setSize(const Coord &size) { Errors::errorMessage("setSize not available for requested face"); };

  Coord vecteur(Element *e);   /*!< Cree vecteur entre center face et center d un element */
  double distance(Element *e); /*!< Calcul de la distance a un center d element */

  virtual void printInfo() const{ Errors::errorMessage("AfficheInfos not available for requested face"); };

  virtual const double& getSizeX() { Errors::errorMessage("getSizeX not available for requested face"); return Errors::defaultDouble; };
  virtual const double& getSizeY() { Errors::errorMessage("getSizeY not available for requested face"); return Errors::defaultDouble; };
  virtual const double& getSizeZ() { Errors::errorMessage("getSizeZ not available for requested face"); return Errors::defaultDouble; };
  virtual const Coord& getSize() { Errors::errorMessage("getSize not available for requested face"); return Coord::defaultCoord; };

  //Pour methode AMR
  virtual Face* creerNouvelleFace() { Errors::errorMessage("creerNouvelleFace not available for requested face"); return 0; };

protected:

  Coord m_position;     /*!< Position du center de la face */
  double m_surface;     /*!< 1.0 pour element 0D, longueur pour element 1D, surface pour element 2D */
  Coord m_normal;
  Coord m_tangent;
  Coord m_binormal;

};

#endif // FACE_H
