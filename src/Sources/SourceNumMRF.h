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

#ifndef SOURCENUMMRF_H
#define SOURCENUMMRF_H

#include "SourceNum.h"

//! \class     SourceNumMRF
//! \brief     Moving Reference Frame method for rotating flows
class SourceNumMRF : public SourceNum
{
public:
  //! \brief     SourceNumMRF constructor from a XML format reading
  //! \details   Reading data from XML file under the following format:
  //!            ex: <dataMRF omega="1.d3"/>
  //! \param     element          XML element to read for source term
  //! \param     fileName         string name of readed XML file
  SourceNumMRF(tinyxml2::XMLElement* element, int order, int physicalEntity, std::string fileName = "Unknown file");
  virtual ~SourceNumMRF();

  virtual void prepSourceTerms(Cell* cell, const int& i=0 );
  virtual void sourceEvolution(const double& time);

  virtual Coord computeAbsVelocity(const Coord& relVelocity, const Coord& position);

private:
  Coord m_omega;    //!Angular velocity
  double m_tf;      //!Optional final time to increase linearly omega
  double m_incr;    //!To increment angular velocity from zero to omega
};

#endif //SOURCENUMMRF_H
