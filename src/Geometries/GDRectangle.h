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

#ifndef GDRECTANGLE_H
#define GDRECTANGLE_H

//! \file      GDRectangle.h
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 19 2017

#include "GeometricalDomain.h"

//! \class     GDRectangle
//! \brief     Class describing a rectangle geometrical domain
class GDRectangle :
  public GeometricalDomain
{
public:
  //! \brief     Geometrical constructor from a XML format reading
  //! \details   Reading data from XML file under the following format:
  //!           ex : <dataRectangle axis1 = "x" axis2 = "y" lAxis1 = "0.3" lAxis2 = "0.2">
  //!                  <posInferiorVertex x = "0.4" y = "0.5" z = "0."/>
  //!                </dataRectangle>
  //! \param     vecPhases      Phases vector variables to copy in geometrical domain
  //! \param     mixture        Mixture variables to copy in geometrical domain
  //! \param     vecTransports  Transports vector varaiables to copy in geometrical domain
  //! \param     element        XML element to read for geometrical properties
  //! \param     physicalEntity physical entity number relative to mesh generation (see mesh tool)
  //! \param     fileName       String name of readed XML file
  GDRectangle(std::string name, std::vector<Phase*> vecPhases, Mixture *mixture, std::vector<Transport> vecTransports, tinyxml2::XMLElement *element, const int &physicalEntity, std::string fileName = "Fichier Inconnu");
  virtual ~GDRectangle();

  virtual bool belong(Coord &posElement, const int &lvl) const;
private:
  Coord m_posLeftBottom;       //!< Coordinates of left bottom corner (minimum positions in X, Y, Z)
  Axis m_axis1, m_axis2;          //!< Axes defining rectangle plane
  double m_lAxis1, m_lAxis2;     //!< Width along axes 1 and 2
};

#endif //GDRECTANGLE_H