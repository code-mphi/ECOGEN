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

//! \file      BoundCondOutflow.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 13 2019

#include "BoundCondOutflow.h"

using namespace tinyxml2;

//****************************************************************************

BoundCondOutflow::BoundCondOutflow(){}

//****************************************************************************

BoundCondOutflow::BoundCondOutflow(int numPhysique, XMLElement *element, int &numberPhases, int &numberTransports, std::vector<std::string> nameTransports, std::string fileName) :
  BoundCond(numPhysique)
{
  //Lecture de the pressure en sortie
  XMLElement *sousElement(element->FirstChildElement("dataOutflow"));
  if (sousElement == NULL) throw ErrorXMLElement("dataOutflow", fileName, __FILE__, __LINE__);
  //Recuperation des attributs
  //--------------------------
  XMLError error;
  error = sousElement->QueryDoubleAttribute("p0", &m_p0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("p0", fileName, __FILE__, __LINE__);

  //Lecture des transports
  int couleurTrouvee(0);
  m_valueTransport = new double[numberTransports];
  XMLElement *elementTransport(sousElement->FirstChildElement("transport"));
  std::string nameTransport;
  while (elementTransport != NULL)
  {
    nameTransport = elementTransport->Attribute("name");
    if (nameTransport == "") throw ErrorXMLAttribut("name", fileName, __FILE__, __LINE__);
    int e(0);
    for (e = 0; e < numberTransports; e++) {
      if (nameTransport == nameTransports[e]) { break; }
    }
    if (e != numberTransports) {
      error = elementTransport->QueryDoubleAttribute("value", &m_valueTransport[e]);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("value", fileName, __FILE__, __LINE__);
      couleurTrouvee++;
    }
    //Transport suivant
    elementTransport = elementTransport->NextSiblingElement("transport");
  }
  if (numberTransports > couleurTrouvee) throw ErrorXMLAttribut("Pas assez d equations de tansport dans CL inj", fileName, __FILE__, __LINE__);
  m_numberTransports = numberTransports;

  //Allocation pour stocker les debits
  m_debits = new double[numberPhases];
  m_numberPhases = numberPhases;
}

//****************************************************************************

BoundCondOutflow::BoundCondOutflow(double p0) : m_p0(p0){}

//****************************************************************************

BoundCondOutflow::BoundCondOutflow(const BoundCondOutflow& Source, const int lvl) : BoundCond(Source)
{
  m_numberPhases = Source.m_numberPhases;
  m_numberTransports = Source.m_numberTransports;

  m_p0 = Source.m_p0;

  m_valueTransport = new double[m_numberTransports];
  for (int k = 0; k < m_numberTransports; k++) {
    m_valueTransport[k] = Source.m_valueTransport[k];
  }
  
  m_debits = new double[m_numberPhases];
  for (int k = 0; k < m_numberPhases; k++) {
    m_debits[k] = 0.;
  }

  m_lvl = lvl;
}

//****************************************************************************

BoundCondOutflow::~BoundCondOutflow()
{
  delete[] m_valueTransport;
  delete[] m_debits;
}

//****************************************************************************

void BoundCondOutflow::creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondOutflow(*(this)));
}

//****************************************************************************

void BoundCondOutflow::solveRiemannLimite(Cell &cellLeft, const int & numberPhases, const double & dxLeft, double & dtMax)
{
  m_mod->solveRiemannOutflow(cellLeft, numberPhases, dxLeft, dtMax, m_p0, m_debits);
  for (int k = 0; k < numberPhases; k++) {
    m_debits[k] *= this->getFace()->getSurface();
    //if (1) m_debits[k] *= 3.14*2.*this->getFace()->getPos().getY();
  }
}

//****************************************************************************

void BoundCondOutflow::solveRiemannTransportLimite(Cell &cellLeft, const int & numberTransports) const
{
	m_mod->solveRiemannTransportOutflow(cellLeft, numberTransports, m_valueTransport);
}

//****************************************************************************

void BoundCondOutflow::printInfo()
{
  std::cout << m_numPhysique << std::endl;
  std::cout << m_p0 << std::endl;
}

//****************************************************************************

//double BoundCondOutflow::getDebit(int numPhase) const 
//{
//  return m_debits[numPhase];
//}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void BoundCondOutflow::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondOutflow(*this, m_lvl + 1));
}

//****************************************************************************