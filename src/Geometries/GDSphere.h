//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//      Official webSite: https://code-mphi.github.io/ECOGEN/
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

#ifndef GDSPHERE_H
#define GDSPHERE_H

#include "GeometricalDomain.h"

//! \class     GDSphere
//! \brief     Class describing a sphere geometrical domain
class GDSphere :
  public GeometricalDomain
{
public:
  //! \brief     Geometrical constructor from a XML format reading
  //! \details   Reading data from XML file under the following format:
  //!             ex :  <dataSphere radius = "0.5">
  //!                     <center x = "1." y = "0.5" z = "0.5"/>
  //!                   </dataSphere>
  //! \param     vecPhases      Phases vector variables to copy in geometrical domain
  //! \param     mixture        Mixture variables to copy in geometrical domain
  //! \param     vecTransports  Transports vector varaiables to copy in geometrical domain
  //! \param     element        XML element to read for geometrical properties
  //! \param     physicalEntity physical entity number relative to mesh generation (see mesh tool)
  //! \param     fileName       String name of readed XML file
  GDSphere(std::string name, std::vector<Phase*> vecPhases, Mixture* mixture, std::vector<Transport> vecTransports, tinyxml2::XMLElement* element, const int& physicalEntity, std::string fileName = "Fichier Inconnu");
  virtual ~GDSphere();

  virtual bool belong(Coord& posElement, const int& /*lvl*/) const;
private:
  Coord m_centerPos;         //!< Sphere position center
  double m_radius;           //!< Sphere radius
};

#endif //GDSPHERE_H