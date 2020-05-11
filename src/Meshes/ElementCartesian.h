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

#ifndef ELEMENTCARTESIAN_H
#define ELEMENTCARTESIAN_H

//! \file      ElementCartesian.h
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include "Element.h"

class ElementCartesian : public Element
{
public:
  ElementCartesian();
  virtual ~ElementCartesian();

  virtual void setVolume(const double &volume);
  virtual void setLCFL(const double &lCFL);
  virtual void setPos(const double &X, const double &Y, const double &Z);
  virtual void setPos(const Coord &pos);
  virtual void setPosX(const double &X);
  virtual void setPosY(const double &Y);
  virtual void setPosZ(const double &Z);
  virtual void setSize(const double &sizeX, const double &sizeY, const double &sizeZ);
  virtual void setSize(const Coord &size);

  virtual const double& getSizeX() { return m_size.getX(); };
  virtual const double& getSizeY() { return m_size.getY(); };
  virtual const double& getSizeZ() { return m_size.getZ(); };
  virtual const Coord& getSize() { return m_size; };

  //Pour methode AMR
  virtual void creerElementChild();
  virtual Element* getElementChild(const int &numberChild);
  virtual Element* getElementChildBack();
  virtual void finalizeElementsChildren();

protected:

  Coord m_size;     //!< dimensions of cartesian cell

  //Attributs pour methode AMR
  std::vector<ElementCartesian*> m_elementsChildren;      /*!< Vector d'elements enfants */

private:
};

#endif // ELEMENTCARTESIAN_H
