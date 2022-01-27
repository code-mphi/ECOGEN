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

#ifndef SOURCENUMGRAVITY_H
#define SOURCENUMGRAVITY_H

#include "SourceNum.h"

//! \class     SourceGravity
//! \brief     Class for gravity source terms
class SourceNumGravity : public SourceNum
{
  public:
    //! \brief     SourceNumGravity constructor from a XML format reading
    //! \details   Reading data from XML file under the following format:
    //!            ex: <gravity x="0." y="-9.81" z="0." / >
    //! \param     element          XML element to read for source term
    //! \param     fileName         string name of readed XML file
    SourceNumGravity(tinyxml2::XMLElement* element, int order, int physicalEntity, std::string fileName = "Unknown file");
    virtual ~SourceNumGravity();

    virtual void prepSourceTerms(Cell* /*cell*/, const int& i = 0);

  private:
    Coord m_g;       //! Gravity acceleration vector
};

#endif // SOURCENUMGRAVITY_H
