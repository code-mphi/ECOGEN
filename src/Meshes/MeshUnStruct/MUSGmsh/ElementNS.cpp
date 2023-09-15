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

#include "ElementNS.h"

//***********************************************************************

ElementNS::ElementNS(){}

//***********************************************************************

ElementNS::ElementNS(const int& typeGmsh, const int& numberNodes, const int& numberFaces, const int& typeVTK) :
m_typeGmsh(typeGmsh),
m_typeVTK(typeVTK),
m_numberNodes(numberNodes),
m_numberFaces(numberFaces),
m_isFantome(false),
m_isCommunicant(false),
m_otherCPU(0)
{
  m_numNodes = new int[numberNodes];
}

//***********************************************************************

ElementNS::~ElementNS()
{
  delete[] m_numNodes;
  if (m_otherCPU != 0) delete[] m_otherCPU;
}

//***********************************************************************

void ElementNS::construitElement(const int* numNodes, const Coord* nodes, const int numberEntitePhysique, const int numberEntiteGeometrique, int& indexElement)
{
  m_index = indexElement;
  //Attribution des number de node vis a vis du array de node global et compute position du Centre de l'element
  m_position = 0.;
  for (int i = 0; i < m_numberNodes; i++){
    m_numNodes[i] = numNodes[i];
    m_position += nodes[i];
  }
  m_position /= static_cast<double>(m_numberNodes);

  m_appartenancePhysique = numberEntitePhysique;
  m_appartenanceGeometrique = numberEntiteGeometrique;

  //Calculs des proprietes de l'element
  this->computeVolume(nodes);
  this->computeLCFL(nodes);
}

//***********************************************************************

void ElementNS::construitElementParallele(const Coord* nodes)
{
  Coord* nodeslocal = new Coord[m_numberNodes];
  for (int i = 0; i < m_numberNodes; i++){ nodeslocal[i] = nodes[m_numNodes[i]]; }

  //Attribution des number de node vis a vis du array de node global et compute position du Centre de l'element
  for (int i = 0; i < m_numberNodes; i++)
  {
    m_position += nodeslocal[i];
  }
  m_position /= static_cast<double>(m_numberNodes);

  //Calculs des proprietes de l'element
  this->computeVolume(nodeslocal);
  this->computeLCFL(nodeslocal);
}

//***********************************************************************

void ElementNS::setIndex(int& index)
{
  m_index = index;
}

//***********************************************************************

void ElementNS::setAppartenancePhysique(int& appartenancePhysique)
{
  m_appartenancePhysique = appartenancePhysique;
}

//***********************************************************************

void ElementNS::setNumNode(int* numNodes)
{
  for (int i = 0; i < m_numberNodes; i++){ m_numNodes[i] = numNodes[i]; }
}

//***********************************************************************

void ElementNS::setNumNode(int& node, int& numNode)
{
  m_numNodes[node] = numNode;
}

//***********************************************************************

void ElementNS::setIsFantome(bool isFantome)
{
  m_isFantome = isFantome;
}

//***********************************************************************

void ElementNS::setIsCommunicant(bool isCommunicant)
{
  m_isCommunicant = isCommunicant;
}

//***********************************************************************

void ElementNS::setAppartenanceCPU(const int* numCPU, const int& numberCPU)
{
  m_CPU = numCPU[0] - 1;
  m_numberOtherCPU = numberCPU - 1;
  m_otherCPU = new int[m_numberOtherCPU];
  for (int i = 1; i < numberCPU; i++)
  {
    m_otherCPU[i - 1] = -numCPU[i] - 1;
  }
}

//***********************************************************************

void ElementNS::removeCPUOthers(std::vector<int>& numCPU)
{
  //if (numCPU >= m_numberOtherCPU){ Errors::errorMessage("Probleme dans removeCPUOthers"); }
  //Copie des anciens
  int* otherCPUTemp = new int[m_numberOtherCPU];
  for (int i = 0; i < m_numberOtherCPU; i++)
  {
    otherCPUTemp[i] = m_otherCPU[i];
  }

  //reperage des CPU a remove
  bool *removeCPU = new bool[m_numberOtherCPU];
  for (int i = 0; i < m_numberOtherCPU; i++)
  {
    removeCPU[i] = false;
    for (unsigned int p = 0; p < numCPU.size(); p++)
    {
      if (i == numCPU[p]){ removeCPU[i] = true; break; }
    }
  }

  //Construction nouveau
  delete[] m_otherCPU;
  m_otherCPU = new int[m_numberOtherCPU - numCPU.size()];
  int indexNouveau(0);
  for (int i = 0; i < m_numberOtherCPU; i++)
  {
    if (!removeCPU[i]){ m_otherCPU[indexNouveau++] = otherCPUTemp[i]; }
  }
  m_numberOtherCPU = indexNouveau;

  delete[] otherCPUTemp;
  delete[] removeCPU;
}

//***********************************************************************

const int& ElementNS::getAutreCPU(const int& autreCPU) const
{
  if (sizeof(m_otherCPU) <= (unsigned int)autreCPU)
  {
    Errors::errorMessage("probleme de dimension dans m_otherCPU");
  }
  return m_otherCPU[autreCPU];
}

//***********************************************************************

void ElementNS::printInfo() const
{
  std::cout << "-------------" << std::endl;
  //cout << "Infos Element" << endl;
  //for (int i = 0; i < m_numberNodes; i++)
  //  cout << " " << m_numNodes[i];
  //cout << endl;
  std::cout << "center : " << m_position.getX() << " " << m_position.getY() << " " << m_position.getZ() << std::endl;
  //cout << " volume : " << m_volume << endl;
  //cout << " lCfl : " << m_lCFL << endl;
  //cout << " CPU n : " << m_CPU << endl;
  //for (int i = 0; i < m_numberOtherCPU; i++)
  //  cout << " autre CPU : " << m_otherCPU[i] << endl;
}

//***********************************************************************