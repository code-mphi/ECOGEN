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

//! \file      ElementNS.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 20 2017

#include "ElementNS.h"

using namespace std;

//***********************************************************************

ElementNS::ElementNS(){}

//***********************************************************************

ElementNS::ElementNS(const int &typeGmsh, const int &numberNoeuds, const int &numberFaces, const int &typeVTK) :
m_typeGmsh(typeGmsh),
m_numberNoeuds(numberNoeuds),
m_numberFaces(numberFaces),
m_typeVTK(typeVTK),
m_isFantome(false),
m_isCommunicant(false)
{
  m_numNoeuds = new int[numberNoeuds];
}

//***********************************************************************

ElementNS::~ElementNS()
{
  delete[] m_numNoeuds;
  delete[] m_autresCPU;
}

//***********************************************************************

void ElementNS::construitElement(const int *numNoeuds, const Coord *noeuds, const int numberEntitePhysique, const int numberEntiteGeometrique, int &indexElement)
{
  m_index = indexElement;
  //Attribution des number de noeud vis a vis du tableau de noeud global et compute position du Centre de l'element
  m_position = 0.;
  for (int i = 0; i < m_numberNoeuds; i++){
    m_numNoeuds[i] = numNoeuds[i];
    m_position += noeuds[i];
  }
  m_position /= static_cast<double>(m_numberNoeuds);

  m_appartenancePhysique = numberEntitePhysique;
  m_appartenanceGeometrique = numberEntiteGeometrique;

  //Calculs des proprietes de l'element
  this->computeVolume(noeuds);
  this->computeLCFL(noeuds);
}

//***********************************************************************

void ElementNS::construitElementParallele(const Coord *noeuds)
{
  Coord *noeudslocal = new Coord[m_numberNoeuds];
  for (int i = 0; i < m_numberNoeuds; i++){ noeudslocal[i] = noeuds[m_numNoeuds[i]]; }

  //Attribution des number de noeud vis a vis du tableau de noeud global et compute position du Centre de l'element
  for (int i = 0; i < m_numberNoeuds; i++)
  {
    m_position += noeudslocal[i];
  }
  m_position /= static_cast<double>(m_numberNoeuds);

  //Calculs des proprietes de l'element
  this->computeVolume(noeudslocal);
  this->computeLCFL(noeudslocal);
}

//***********************************************************************

void ElementNS::setIndex(int &index)
{
  m_index = index;
}

//***********************************************************************

void ElementNS::setAppartenancePhysique(int &appartenancePhysique)
{
  m_appartenancePhysique = appartenancePhysique;
}

//***********************************************************************

void ElementNS::setNumNoeud(int *numNoeuds)
{
  for (int i = 0; i < m_numberNoeuds; i++){ m_numNoeuds[i] = numNoeuds[i]; }
}

//***********************************************************************

void ElementNS::setNumNoeud(int &noeud, int &numNoeud)
{
  m_numNoeuds[noeud] = numNoeud;
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

void ElementNS::setAppartenanceCPU(const int *numCPU, const int &numberCPU)
{
  m_CPU = numCPU[0] - 1;
  m_numberautresCPU = numberCPU - 1;
  m_autresCPU = new int[m_numberautresCPU];
  for (int i = 1; i < numberCPU; i++)
  {
    m_autresCPU[i - 1] = -numCPU[i] - 1;
  }
}

//***********************************************************************

void ElementNS::enleveCPUAutres(vector<int> &numCPU)
{
  //if (numCPU >= m_numberautresCPU){ Errors::errorMessage("Probleme dans enleveCPUAutres"); }
  //Copie des anciens
  int *autresCPUTemp = new int[m_numberautresCPU];
  for (int i = 0; i < m_numberautresCPU; i++)
  {
    autresCPUTemp[i] = m_autresCPU[i];
  }

  //reperage des CPU a enlever
  bool *enleveCPU = new bool[m_numberautresCPU];
  for (int i = 0; i < m_numberautresCPU; i++)
  {
    enleveCPU[i] = false;
    for (unsigned int p = 0; p < numCPU.size(); p++)
    {
      if (i == numCPU[p]){ enleveCPU[i] = true; break; }
    }
  }

  //Construction nouveau
  delete[] m_autresCPU;
  m_autresCPU = new int[m_numberautresCPU - numCPU.size()];
  int indexNouveau(0);
  for (int i = 0; i < m_numberautresCPU; i++)
  {
    if (!enleveCPU[i]){ m_autresCPU[indexNouveau++] = autresCPUTemp[i]; }
  }
  m_numberautresCPU = indexNouveau;

  delete[] autresCPUTemp;
  delete[] enleveCPU;
}

//***********************************************************************

int ElementNS::getIndex() const
{
  return m_index;
}

//***********************************************************************

int ElementNS::getNumberNoeuds() const
{
  return m_numberNoeuds;
}

//***********************************************************************

int ElementNS::getNumberFaces() const
{
  return m_numberFaces;
}

//***********************************************************************

int ElementNS::getTypeGmsh() const
{
  return m_typeGmsh;
}

//***********************************************************************

int ElementNS::getTypeVTK() const
{
  return m_typeVTK;
}

//***********************************************************************

int ElementNS::getNumNoeud(int &noeud) const
{
  return m_numNoeuds[noeud];
}

//***********************************************************************

int ElementNS::getAppartenancePhysique() const
{
  return m_appartenancePhysique;
}

//***********************************************************************

int ElementNS::getAppartenanceGeometrique() const
{
  return m_appartenanceGeometrique;
}

//***********************************************************************

int ElementNS::getCPU() const
{
  return m_CPU;
}

//***********************************************************************

int ElementNS::getNumberAutresCPU() const
{
  return m_numberautresCPU;
}

//***********************************************************************

int ElementNS::getAutreCPU(const int &autreCPU) const
{
  if (sizeof(m_autresCPU) <= autreCPU)
  {
    Errors::errorMessage("probleme de dimension dans m_autresCPU");
  }
  return m_autresCPU[autreCPU];
}

//***********************************************************************

bool ElementNS::isFantome() const
{
  return m_isFantome;
}

//***********************************************************************

bool ElementNS::isCommunicant() const
{
  return m_isCommunicant;
}

//***********************************************************************

void ElementNS::printInfo() const
{
  cout << "-------------" << endl;
  //cout << "Infos Element" << endl;
  //for (int i = 0; i < m_numberNoeuds; i++)
  //  cout << " " << m_numNoeuds[i];
  //cout << endl;
  cout << "center : " << m_position.getX() << " " << m_position.getY() << " " << m_position.getZ() << endl;
  //cout << " volume : " << m_volume << endl;
  //cout << " lCfl : " << m_lCFL << endl;
  //cout << " CPU n : " << m_CPU << endl;
  //for (int i = 0; i < m_numberautresCPU; i++)
  //  cout << " autre CPU : " << m_autresCPU[i] << endl;
}

//***********************************************************************