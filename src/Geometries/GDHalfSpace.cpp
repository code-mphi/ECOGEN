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

#include "GDHalfSpace.h"

using namespace tinyxml2;

//***************************************************************

GDHalfSpace::GDHalfSpace(std::string name, std::vector<Phase*> vecPhases, Mixture* mixture, std::vector<Transport> vecTransports, XMLElement* element, const int& physicalEntity, std::string fileName) :
  GeometricalDomain(name, vecPhases, mixture, vecTransports, physicalEntity)
{
  XMLElement* sousElement(element->FirstChildElement("dataHalfSpace"));
  if (sousElement == NULL) throw ErrorXMLElement("dataHalfSpace", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  //Origin
  error = sousElement->QueryDoubleAttribute("origin", &m_position);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("origin", fileName, __FILE__, __LINE__);
  //Axis
  std::string axis(sousElement->Attribute("axis"));
  Tools::uppercase(axis);
  if      (axis == "X"){ m_axis = X; }
  else if (axis == "Y"){ m_axis = Y; }
  else if (axis == "Z"){ m_axis = Z; }
  else { throw ErrorXMLAttribut("axis", fileName, __FILE__, __LINE__); }
  //Direction
  std::string direction(sousElement->Attribute("direction"));
  Tools::uppercase(direction);
  if      (direction == "POSITIVE"){ m_direction = 1; }
  else if (direction == "NEGATIVE"){ m_direction = -1; }
  else { throw ErrorXMLAttribut("direction", fileName, __FILE__, __LINE__); }

  //To uncomment only for special test cases
  //4. Random velocity perturbations: O(1e−4 u_s)
  //---------------------------------------------
  //Inialize random number generator
  // int rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // srand(rank);
}

//***************************************************************

GDHalfSpace::~GDHalfSpace(){}

//***************************************************************

bool GDHalfSpace::belong(Coord& posElement, const int& /*lvl*/) const
{
  bool result(false);
  switch (m_axis)
  {
  case X:
    if (m_direction >= 0){if (posElement.getX() >= m_position) result = true;}
    else{ if (posElement.getX() <= m_position) result = true; }
    break;
  case Y:
    if (m_direction >= 0){ if (posElement.getY() >= m_position) result = true; }
    else{ if (posElement.getY() <= m_position) result = true; }
    break;
  case Z:
    if (m_direction >= 0){ if (posElement.getZ() >= m_position) result = true; }
    else{ if (posElement.getZ() <= m_position) result = true; }
    break;
  default:
    result = false; //false if not belongs to
    break;
  }
  return result;
}

//***************************************************************

void GDHalfSpace::fillIn(Cell* cell) const
{
  //As basic fillIn: Test if the cell belongs to the geometrical domain
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
    if(m_physicalEntity == -1){ cell->setWall(true); }
    else{ cell->setWall(false); }

    //To uncomment only for special test cases
    //4. Random velocity perturbations: O(1e−4 u_s)
    //---------------------------------------------
    // Coord perturbedVelocity(cell->getMixture()->getVelocity());
    // perturbedVelocity.setX(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getX());
    // perturbedVelocity.setY(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getY());
    // perturbedVelocity.setZ(static_cast<double>(rand() % 2001 - 1000)/1.e3 * 1.e-3*151.821433232719 + perturbedVelocity.getZ());
    // cell->getMixture()->setVelocity(perturbedVelocity);
  }
}

//***************************************************************