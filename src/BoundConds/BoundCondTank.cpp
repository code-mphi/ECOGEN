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

//! \file      BoundCondTank.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 13 2019

#include "BoundCondTank.h"

using namespace tinyxml2;

//****************************************************************************

BoundCondTank::BoundCondTank(){}

//****************************************************************************

BoundCondTank::BoundCondTank(int numPhysique, XMLElement *element, int &numberPhases, int &numberTransports, std::vector<std::string> nameTransports, Eos **eos, std::string fileName) :
  BoundCond(numPhysique)
{
  m_numberPhase = numberPhases;
  m_ak0 = new double[m_numberPhase];
  m_Yk0 = new double[m_numberPhase];
  m_rhok0 = new double[m_numberPhase]; 
  
  //Reading tank pressure and temperature conditions
  //------------------------------------------------
  XMLElement *sousElement(element->FirstChildElement("dataTank"));
  if (sousElement == NULL) throw ErrorXMLElement("dataTank", fileName, __FILE__, __LINE__);
  //Attributes reading
  XMLError error;
  error = sousElement->QueryDoubleAttribute("p0", &m_p0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("p0", fileName, __FILE__, __LINE__);
  error = sousElement->QueryDoubleAttribute("T0", &m_T0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("T0", fileName, __FILE__, __LINE__);

  if (m_numberPhase == 1) {
    m_rhok0[0] = eos[0]->computeDensity(m_p0, m_T0);
    m_ak0[0] = 1.;
    m_Yk0[0] = 1.;
  }
  else {
    //Reading fluids proportion in tank
    //---------------------------------
    sousElement = element->FirstChildElement("fluidsProp");
    if (sousElement == NULL) throw ErrorXMLElement("fluidsProp", fileName, __FILE__, __LINE__);
    XMLElement* fluid(sousElement->FirstChildElement("dataFluid"));

    int nbFluids(0); std::string nameEOS;
    int presenceAlpha(0), presenceMassFrac(0);
    while (fluid != NULL)
    {
      nbFluids++;
      //EOS name searching
      nameEOS = fluid->Attribute("EOS");
      int e(0);
      for (e = 0; e < m_numberPhase; e++) {
        if (nameEOS == eos[e]->getName()) { break; }
      }
      if (e == m_numberPhase) { throw ErrorXMLEOSInconnue(nameEOS, fileName, __FILE__, __LINE__); }

      //Reading fluid proportion
      if (fluid->QueryDoubleAttribute("alpha", &m_ak0[e]) == XML_NO_ERROR) {
        if (presenceMassFrac) throw ErrorXMLAttribut("only one of following is required : alpha, massFrac", fileName, __FILE__, __LINE__);
        presenceAlpha = 1;
      }
      if (fluid->QueryDoubleAttribute("massFrac", &m_Yk0[e]) == XML_NO_ERROR) {
        if (presenceAlpha) throw ErrorXMLAttribut("only one of following is required : alpha, massFrac", fileName, __FILE__, __LINE__);
        presenceMassFrac = 1;
      }
      fluid = fluid->NextSiblingElement("dataFluid");
    }
    if (nbFluids != m_numberPhase) throw ErrorXMLEtat("Tank", fileName, __FILE__, __LINE__);

    //Proportions checking
    //--------------------
    double sum(0.);
    if (presenceAlpha) {
      for (int k = 0; k < m_numberPhase; k++) {
        if (m_ak0[k]<0. || m_ak0[k]>1.) throw ErrorXMLAttribut("alpha should be in [0,1]", fileName, __FILE__, __LINE__);
        sum += m_ak0[k];
      }
      if (std::fabs(sum - 1.) > 1.e-6) { throw ErrorXMLAttribut("sum of alpha should be 1", fileName, __FILE__, __LINE__); }
      else {
        for (int k = 0; k < m_numberPhase; k++) { m_ak0[k] /= sum; }
      }
    }
    else if (presenceMassFrac) {
      for (int k = 0; k < m_numberPhase; k++) {
        if (m_Yk0[k]<0. || m_Yk0[k]>1.) throw ErrorXMLAttribut("massFrac should be in [0,1]", fileName, __FILE__, __LINE__);
        sum += m_Yk0[k];
      }
      if (std::fabs(sum - 1.) > 1.e-6) { throw ErrorXMLAttribut("sum of massFrac should be 1", fileName, __FILE__, __LINE__); }
      else {
        for (int k = 0; k < m_numberPhase; k++) { m_Yk0[k] /= sum; }
      }
    }
    else { throw ErrorXMLAttribut("One of following is required : alpha, massFrac", fileName, __FILE__, __LINE__); }

    //Fulfill tank state (rhok0, ak0 or Yk0)
    //--------------------------------------
    for (int k = 0; k < m_numberPhase; k++) {
      m_rhok0[k] = eos[k]->computeDensity(m_p0, m_T0);
    }
    double rhoMel(0.);
    if (presenceAlpha) {
      for (int k = 0; k < m_numberPhase; k++) { rhoMel += m_ak0[k] * m_rhok0[k]; }
      for (int k = 0; k < m_numberPhase; k++) { m_Yk0[k] = m_ak0[k] * m_rhok0[k] / rhoMel; }
    }
    else {
      for (int k = 0; k < m_numberPhase; k++) { rhoMel += m_Yk0[k] / m_rhok0[k]; }
      rhoMel = 1.0 / rhoMel;
      for (int k = 0; k < m_numberPhase; k++) { m_ak0[k] = rhoMel * m_Yk0[k] / m_rhok0[k]; }
    }

  } //End proportion

  //Reading of transports
  //---------------------
  if (numberTransports) {
    XMLElement *sousElement(element->FirstChildElement("dataTank"));
    if (sousElement == NULL) throw ErrorXMLElement("dataTank", fileName, __FILE__, __LINE__);
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
    if (numberTransports > foundColors) throw ErrorXMLAttribut("Not enough transport equations in tank BC", fileName, __FILE__, __LINE__);
  }
  m_numberTransports = numberTransports;
}

//****************************************************************************

BoundCondTank::BoundCondTank(const BoundCondTank &Source, const int lvl) : BoundCond(Source)
{
  m_numberPhase = Source.m_numberPhase;
  m_numberTransports = Source.m_numberTransports;
  m_ak0 = new double[m_numberPhase];
  m_Yk0 = new double[m_numberPhase];
  m_rhok0 = new double[m_numberPhase];

  for (int k = 0; k < m_numberPhase; k++)
  {
    m_ak0[k] = Source.m_ak0[k];
    m_Yk0[k] = Source.m_Yk0[k];
    m_rhok0[k] = Source.m_rhok0[k];
  }
  m_p0 = Source.m_p0;
  m_T0 = Source.m_T0;

  m_valueTransport = new double[Source.m_numberTransports];
  for (int k = 0; k < Source.m_numberTransports; k++) {
    m_valueTransport[k] = Source.m_valueTransport[k];
  }

  m_lvl = lvl;
}

//****************************************************************************

BoundCondTank::~BoundCondTank()
{
  delete[] m_ak0;
  delete[] m_Yk0;
  delete[] m_rhok0;
  delete[] m_valueTransport;
}

//****************************************************************************

void BoundCondTank::creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondTank(*(this)));
}

//****************************************************************************

void BoundCondTank::solveRiemannLimite(Cell &cellLeft, const int & numberPhases, const double & dxLeft, double & dtMax)
{
  Coord omega(0., 0., 500.);
  m_mod->solveRiemannTank(cellLeft, numberPhases, dxLeft, dtMax, m_ak0, m_rhok0, m_p0, m_T0);
  //std::cout << m_face->getPos().getX() << " " << m_face->getPos().getY() << " " << m_face->getPos().getZ() << std::endl;
}

//****************************************************************************

void BoundCondTank::solveRiemannTransportLimite(Cell &cellLeft, const int & numberTransports) const
{
	m_mod->solveRiemannTransportTank(cellLeft, numberTransports, m_valueTransport);
}

//****************************************************************************

void BoundCondTank::printInfo()
{
  std::cout << m_numPhysique << std::endl;
  std::cout << m_rhok0[0] << std::endl;
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void BoundCondTank::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondTank(*this, m_lvl + 1));
}

//****************************************************************************
