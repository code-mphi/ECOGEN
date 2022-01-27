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

#include "GeometricalDomain.h"

//******************************************************************

GeometricalDomain::GeometricalDomain(std::string name, std::vector<Phase*> vecPhases, Mixture* mixture, std::vector<Transport> vecTransports, const int& physicalEntity) : m_name(name), m_physicalEntity(physicalEntity)
{
  int numbPhases = vecPhases.size();
  int numbTransports = vecTransports.size();
  m_vecPhases = new Phase*[numbPhases];
  for (int k = 0; k < numbPhases; k++){
    vecPhases[k]->allocateAndCopyPhase(&m_vecPhases[k]);
  }
  mixture->allocateAndCopyMixture(&m_mixture);
  if (numbTransports > 0) { m_vecTransports = new Transport[numbTransports]; }
  for (int k = 0; k < numbTransports; k++) {
    m_vecTransports[k].setValue(vecTransports[k].getValue());
  }
}

//******************************************************************

GeometricalDomain::~GeometricalDomain()
{
  for (int k = 0; k < numberPhases; k++) {
    delete m_vecPhases[k];
  }
  delete[] m_vecPhases;
  delete m_mixture;
  if (numberTransports != 0) delete[] m_vecTransports;
}

//******************************************************************

void GeometricalDomain::fillIn(Cell* cell) const
{
  //Test if the cell belongs to the geometrical domain
  bool belongs(true);
  if (cell->getElement() != 0) {
    Coord coord(cell->getPosition());
    if (!this->belong(coord, cell->getLvl())) { belongs = false; }
    //Test if the cell belongs to physical mesh entity (for unstructured meshes)
    if (cell->getElement()->getAppartenancePhysique() > 0 && m_physicalEntity > 0) {
      if(cell->getElement()->getAppartenancePhysique() != m_physicalEntity) { belongs = false; }
    }
  }

  if (belongs) {
    for (int k = 0; k < numberPhases; k++) { cell->copyPhase(k, m_vecPhases[k]); }
    cell->copyMixture(m_mixture);
    for (int k = 0; k < numberTransports; k++) { cell->setTransport(m_vecTransports[k].getValue(), k); }
  }
}