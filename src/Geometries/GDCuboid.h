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

#ifndef GDCUBOID_H
#define GDCUBOID_H

//! \file      GDCuboid.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 19 2017

#include "GeometricalDomain.h"

//! \class     GDCuboid
//! \brief     Class describing a cuboid geometrical domain
class GDCuboid :
  public GeometricalDomain
{
public:
  //! \brief     Geometrical constructor from a XML format reading
  //! \details   Reading data from XML file under the following format:
  //!        ex :	<dataCuboid lAxisX="0.3" lAxisY="0.2" lAxisZ="0.4">
  //!               <posInferiorVertex x = "0.4" y = "0.5" z = "0." />
  //!             </dataCuboid>
  //! \param     vecPhases      Phases vector variables to copy in geometrical domain
  //! \param     mixture        Mixture variables to copy in geometrical domain
  //! \param     vecTransports  Transports vector varaiables to copy in geometrical domain
  //! \param     element        XML element to read for geometrical properties
  //! \param     physicalEntity physical entity number relative to mesh generation (see mesh tool)
  //! \param     fileName       String name of readed XML file
  GDCuboid(std::string name, std::vector<Phase*> vecPhases, Mixture *mixture, std::vector<Transport> vecTransports, tinyxml2::XMLElement *element, const int &physicalEntity, std::string fileName = "Fichier Inconnu");
  virtual ~GDCuboid();

  virtual bool belong(Coord &posElement, const int &lvl) const;
private:
  Coord m_posXmYmZm;            //!< Coordinates of the bottom left front corner (minimum positions in X, Y, Z)
  double m_lX, m_lY, m_lZ;      //!< Width along X, Y and Z axis
};

#endif //GDCUBOID_H