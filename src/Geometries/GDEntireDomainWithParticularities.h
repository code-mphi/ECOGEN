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

#ifndef GDENTIREDOMAINWITHPARTICULARITIES_H
#define GDENTIREDOMAINWITHPARTICULARITIES_H

//! \file      GDEntireDomainWithParticularities.h
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      January 3 2018

#include "GeometricalDomain.h"

//! \class     GDEntireDomainWithParticularities
//! \brief     Class describing geometrical domain for the entire geometry with particularities involved. Need to be recompiled after each changement of the particularities.
class GDEntireDomainWithParticularities :
  public GeometricalDomain
{
public:
  //! \brief     Geometrical constructor for the entire geometrical domain with particularities involved
  //! \param     vecPhases      Phases vector variables to copy in geometrical domain
  //! \param     mixture        Mixture variables to copy in geometrical domain
  //! \param     vecTransports  Transports vector varaiables to copy in geometrical domain
  //! \param     physicalEntity physical entity number relative to mesh generation (see mesh tool)
  GDEntireDomainWithParticularities(std::string name, std::vector<Phase*> vecPhases, Mixture *mixture, std::vector<Transport> vecTransports, const int &physicalEntity);
  virtual ~GDEntireDomainWithParticularities();

  virtual bool belong(Coord &posElement, const int &lvl) const;
  virtual void fillIn(Cell *cell, const int &numberPhases, const int &numberTransports) const;
};

#endif //GDENTIREDOMAINWITHPARTICULARITIES_H