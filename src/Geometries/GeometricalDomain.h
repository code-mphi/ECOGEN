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

#ifndef GEOMETRICALDOMAIN_H
#define GEOMETRICALDOMAIN_H

#include <string>
#include <vector>
#include <cmath>
#include "../Maths/Coord.h"
#include "../Models/Phase.h"
#include "../libTierces/tinyxml2.h"
#include "../Errors.h"
#include "../Tools.h"

class GeometricalDomain; //pre-declaration of GeometricalDomain class needed for Cell.h inclusion
#include "../Order1/Cell.h"

//! \class     GeometricalDomain
//! \brief     General class for geometrical domain
//! \details   This is a pure virtual class: can not be instantiated
class GeometricalDomain
{
public:
  //! \brief     Generic geometrical constructor
  //! \param     vecPhases         Phases vector variables to copy in geometrical domain
  //! \param     mixture           Mixture variables to copy in geometrical domain
  //! \param     vecTransports     Transports vector varaiables to copy in geometrical domain
  //! \param     physicalEntity physical entity number relative to mesh generation (see mesh tool)
  GeometricalDomain(std::string name, std::vector<Phase*> vecPhases, Mixture* mixture, std::vector<Transport> vecTransports, const int& physicalEntity);
  virtual ~GeometricalDomain();

  //! \brief     Method to verify inclusion of a vertex in geometrical domain
  //! \param     posElement        Point coordinates
  //! \param     lvl               Level of the cell
  //! \return    True if the vertex belongs to geometrical domain
  virtual bool belong(Coord& /*posElement*/, const int& /*lvl*/) const = 0;
  //! \brief     Method to fill in the cell data with the ones of the corresponding domain
  //! \param     cell              Cell
  //! \param     numberPhases      Number of phases
  //! \param     numberTransports  Number of transport
  virtual void fillIn(Cell* cell, const int& numberPhases, const int& numberTransports) const;

  const std::string& getName() { return m_name; };

protected:
  std::string m_name;           //!< Geometrical domain name
  int m_numberPhases;           //!< Phases number
  int m_numberTransports;       //!< Transport equations number
  Phase** m_vecPhases;          //!< Phases variable vector
  Mixture* m_mixture;           //!< Mixture variables
  Transport* m_vecTransports;   //!< Transport variables vector
  int m_physicalEntity;         //!< Physical entity number (see mesh software and mesh file)
};

#endif // GEOMETRICALDOMAIN_H

