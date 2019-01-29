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

#ifndef FACECARTESIAN_H
#define FACECARTESIAN_H

//! \file      FaceCartesian.h
//! \author    F. Petitpas, K.Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "Face.h"
class FaceCartesian : public Face
{
public:
  FaceCartesian();
  virtual ~FaceCartesian();

  virtual void setSurface(const double &surface);
  virtual void initializeAutres(const double &surface, const Coord &normal, const Coord &tangent, const Coord &binormal);
  virtual void setPos(const double &X, const double &Y, const double &Z);
  virtual void setNormal(const double &X, const double &Y, const double &Z);
  virtual void setTangent(const double &X, const double &Y, const double &Z);
  virtual void setBinormal(const double &X, const double &Y, const double &Z);
  virtual void setSize(const double &sizeX, const double &sizeY, const double &sizeZ);
  virtual void setSize(const Coord &size);

  virtual double getSizeX();
  virtual double getSizeY();
  virtual double getSizeZ();
  virtual Coord getSize();

  //Pour methode AMR
  virtual Face* creerNouvelleFace();

protected:

  Coord m_size;    //!< dimensions of cartesian face

};

#endif // FACECARTESIAN_H