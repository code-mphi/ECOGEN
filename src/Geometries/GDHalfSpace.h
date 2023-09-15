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

#ifndef GDHALFSPACE_H
#define GDHALFSPACE_H

#include "GeometricalDomain.h"

//! \class     GDHalfSpace
//! \brief     Class describing a half space geometrical domain
class GDHalfSpace : public GeometricalDomain
{
  public:
    //! \brief     Geometrical constructor from a XML format reading
    //! \details   Reading data from XML file under the following format:
    //!            ex : <dataHalfSpace axis="x" origin="0.5" direction="positive"/>
    //! \param     vecPhases      Phases vector variables to copy in geometrical domain
    //! \param     mixture        Mixture variables to copy in geometrical domain
    //! \param     vecTransports  Transports vector varaiables to copy in geometrical domain
    //! \param     element        XML element to read for geometrical properties
    //! \param     physicalEntity physical entity number relative to mesh generation (see mesh tool)
    //! \param     fileName       String name of readed XML file
    GDHalfSpace(std::string name, std::vector<Phase*> vecPhases, Mixture* mixture, std::vector<Transport> vecTransports, tinyxml2::XMLElement* element, const int& physicalEntity, std::string fileName="Unknown file");
    virtual ~GDHalfSpace();

    virtual bool belong(Coord& posElement, const int& /*lvl*/) const;
    virtual void fillIn(Cell* cell) const;

  private:
    double m_position;  //!< Origin of the half space along axis
    Axis m_axis;          //!< Axis orthogonal to the origin plane of the half space
    int m_direction;    //!< Direction along axis (positive or negative)
};

#endif //GDHALFSPACE_H