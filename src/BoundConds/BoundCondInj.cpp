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

//! \file      BoundCondInj.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      May 07 2020

#include "BoundCondInj.h"

using namespace tinyxml2;

//****************************************************************************

BoundCondInj::BoundCondInj(){}

//****************************************************************************

BoundCondInj::BoundCondInj(int numPhysique, XMLElement* element, int& numberPhases, int& numberTransports, std::vector<std::string> nameTransports, Eos** eos, std::string fileName) :
  BoundCond(numPhysique), m_T0(0.)
{
  m_numberPhase = numberPhases;
  m_ak0 = new double[m_numberPhase];
  m_rhok0 = new double[m_numberPhase];
  m_pk0 = new double[m_numberPhase];

  //Reading injection surface-mass-flow condition (kg/s/m²)
  //-------------------------------------------------------
  XMLElement *sousElement(element->FirstChildElement("dataInjection"));
  if (sousElement == NULL) throw ErrorXMLElement("dataInjection", fileName, __FILE__, __LINE__);
  //Attributes reading
  XMLError error;
  error = sousElement->QueryDoubleAttribute("m0", &m_m0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("m0", fileName, __FILE__, __LINE__);
  m_m0 = -m_m0; //Changement de signe car cond limite droite par defaut

  //Reading volume fraction, density and pressure of the phases
  //-----------------------------------------------------------
  if (m_numberPhase == 1) {
    m_ak0[0] = 1.;
    XMLElement* fluid(element->FirstChildElement("dataFluid"));
    if (fluid == NULL) throw ErrorXMLElement("dataFluid", fileName, __FILE__, __LINE__);
    //Attributes reading
    error = fluid->QueryDoubleAttribute("density", &m_rhok0[0]);
    if (error != XML_NO_ERROR) error = fluid->QueryDoubleAttribute("temperature", &m_T0);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("density", fileName, __FILE__, __LINE__);
    error = fluid->QueryDoubleAttribute("pressure", &m_pk0[0]);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("pressure", fileName, __FILE__, __LINE__);

    if(m_T0 !=0.) m_rhok0[0] = eos[0]->computeDensity(m_pk0[0], m_T0);

  }
  else {
    //Reading proportion of inflow fluids
    //-----------------------------------
    XMLElement* fluid(element->FirstChildElement("dataFluid"));

    for(int e = 0; e <m_numberPhase; e++){
      //Attributes reading
      error = fluid->QueryDoubleAttribute("alpha", &m_ak0[e]);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("alpha", fileName, __FILE__, __LINE__);
      error = fluid->QueryDoubleAttribute("density", &m_rhok0[e]);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("density", fileName, __FILE__, __LINE__);
      error = fluid->QueryDoubleAttribute("pressure", &m_pk0[e]);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("pressure", fileName, __FILE__, __LINE__);
      fluid = fluid->NextSiblingElement("dataFluid");
    }

    //Proportions checking
    //--------------------
    double sum(0.);
    for (int k = 0; k < m_numberPhase; k++) {
      if (m_ak0[k]<0. || m_ak0[k]>1.) throw ErrorXMLAttribut("alpha should be in [0,1]", fileName, __FILE__, __LINE__);
      sum += m_ak0[k];
    }
    if (std::fabs(sum - 1.) > 1.e-6) { throw ErrorXMLAttribut("sum of alpha should be 1", fileName, __FILE__, __LINE__); }
    else {
      for (int k = 0; k < m_numberPhase; k++) { m_ak0[k] /= sum; }
    }
  }

  //Reading of transports
  //---------------------
  if (numberTransports) {
    XMLElement *sousElement(element->FirstChildElement("dataInjection"));
    if (sousElement == NULL) throw ErrorXMLElement("dataInjection", fileName, __FILE__, __LINE__);
    XMLError error;

    int foundColors(0);
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
        foundColors++;
      }
      //Following transport
      elementTransport = elementTransport->NextSiblingElement("transport");
    }
    if (numberTransports > foundColors) throw ErrorXMLAttribut("Not enough transport equations in inj BC", fileName, __FILE__, __LINE__);
  }
  m_numberTransports = numberTransports;
}

//****************************************************************************

BoundCondInj::BoundCondInj(const BoundCondInj &Source, const int lvl) : BoundCond(Source)
{
  m_numberPhase = Source.m_numberPhase;
  m_numberTransports = Source.m_numberTransports;
  m_ak0 = new double[m_numberPhase];
  m_rhok0 = new double[m_numberPhase];
  m_pk0 = new double[m_numberPhase];

  m_m0 = Source.m_m0;

  for (int k = 0; k < m_numberPhase; k++)
  {
    m_ak0[k] = Source.m_ak0[k];
    m_rhok0[k] = Source.m_rhok0[k];
    m_pk0[k] = Source.m_pk0[k];
  }

  m_valueTransport = new double[Source.m_numberTransports];
  for (int k = 0; k < Source.m_numberTransports; k++) {
    m_valueTransport[k] = Source.m_valueTransport[k];
  }

  m_lvl = lvl;
}

//****************************************************************************

BoundCondInj::~BoundCondInj()
{
  delete[] m_ak0;
  delete[] m_rhok0;
  delete[] m_pk0;
  delete[] m_valueTransport;
}

//****************************************************************************

void BoundCondInj::creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondInj(*(this)));
}

//****************************************************************************

void BoundCondInj::solveRiemannLimite(Cell &cellLeft, const int & numberPhases, const double & dxLeft, double & dtMax)
{
  m_mod->solveRiemannInflow(cellLeft, numberPhases, dxLeft, dtMax, m_m0, m_ak0, m_rhok0, m_pk0);
}

//****************************************************************************

void BoundCondInj::solveRiemannTransportLimite(Cell &cellLeft, const int & numberTransports) const
{
	m_mod->solveRiemannTransportInflow(cellLeft, numberTransports, m_valueTransport);
}

//****************************************************************************

void BoundCondInj::printInfo()
{
  std::cout << m_numPhysique << std::endl;
  std::cout << m_m0 << std::endl;
  std::cout << m_rhok0[0] << std::endl;
}


//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void BoundCondInj::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondInj(*this, m_lvl + 1));
}

//****************************************************************************