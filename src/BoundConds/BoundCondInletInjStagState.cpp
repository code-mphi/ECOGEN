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

#include "BoundCondInletInjStagState.h"

using namespace tinyxml2;

//****************************************************************************

BoundCondInletInjStagState::BoundCondInletInjStagState(int numPhysique, XMLElement* element, const int& numbPhases, const int& numbTransports, std::vector<std::string> nameTransports, Eos** eos, std::string fileName) :
  BoundCond(numPhysique), m_T0(0.)
{
  m_ak0 = new double[numbPhases];
  m_rhok0 = new double[numbPhases];
  m_pk0 = new double[numbPhases];

  //Reading injection surface-mass-flow condition (kg/s/m�)
  //-------------------------------------------------------
  XMLElement* sousElement(element->FirstChildElement("dataInletInj"));
  if (sousElement == NULL) throw ErrorXMLElement("dataInletInj", fileName, __FILE__, __LINE__);
  //Attributes reading
  XMLError error;
  error = sousElement->QueryDoubleAttribute("m0", &m_m0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("m0", fileName, __FILE__, __LINE__);
  m_m0 = -m_m0; //Changement de signe car cond limite droite par defaut

  //Reading volume fraction, density and pressure of the phases
  //-----------------------------------------------------------
  if (numbPhases == 1) {
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

    for(int e = 0; e < numbPhases; e++){
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
    for (int k = 0; k < numbPhases; k++) {
      if (m_ak0[k]<0. || m_ak0[k]>1.) throw ErrorXMLAttribut("alpha should be in [0,1]", fileName, __FILE__, __LINE__);
      sum += m_ak0[k];
    }
    if (std::fabs(sum - 1.) > 1.e-6) { throw ErrorXMLAttribut("sum of alpha should be 1", fileName, __FILE__, __LINE__); }
    else {
      for (int k = 0; k < numbPhases; k++) { m_ak0[k] /= sum; }
    }
  }

  //Reading of transports
  //---------------------
  m_valueTransport = new double[numbTransports];
  if (numbTransports) {
    XMLElement* sousElement(element->FirstChildElement("dataInletInj"));
    if (sousElement == NULL) throw ErrorXMLElement("dataInletInj", fileName, __FILE__, __LINE__);
    XMLError error;

    int foundColors(0);
    XMLElement* elementTransport(sousElement->FirstChildElement("transport"));
    std::string nameTransport;
    while (elementTransport != NULL)
    {
      nameTransport = elementTransport->Attribute("name");
      if (nameTransport == "") throw ErrorXMLAttribut("name", fileName, __FILE__, __LINE__);
      int e(0);
      for (e = 0; e < numbTransports; e++) {
        if (nameTransport == nameTransports[e]) { break; }
      }
      if (e != numbTransports) {
        error = elementTransport->QueryDoubleAttribute("value", &m_valueTransport[e]);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("value", fileName, __FILE__, __LINE__);
        foundColors++;
      }
      //Next transport
      elementTransport = elementTransport->NextSiblingElement("transport");
    }
    if (numbTransports > foundColors) throw ErrorXMLAttribut("Not enough transport equations in inj BC", fileName, __FILE__, __LINE__);
  }
}

//****************************************************************************

BoundCondInletInjStagState::BoundCondInletInjStagState(const BoundCondInletInjStagState &Source, const int& lvl) : BoundCond(Source, lvl)
{
  m_ak0 = new double[numberPhases];
  m_rhok0 = new double[numberPhases];
  m_pk0 = new double[numberPhases];

  m_m0 = Source.m_m0;

  for (int k = 0; k < numberPhases; k++)
  {
    m_ak0[k] = Source.m_ak0[k];
    m_rhok0[k] = Source.m_rhok0[k];
    m_pk0[k] = Source.m_pk0[k];
  }

  m_valueTransport = new double[numberTransports];
  for (int k = 0; k < numberTransports; k++) {
    m_valueTransport[k] = Source.m_valueTransport[k];
  }
}

//****************************************************************************

BoundCondInletInjStagState::~BoundCondInletInjStagState()
{
  delete[] m_ak0;
  delete[] m_rhok0;
  delete[] m_pk0;
  delete[] m_valueTransport;
}

//****************************************************************************

void BoundCondInletInjStagState::createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondInletInjStagState(*(this)));
}

//****************************************************************************

void BoundCondInletInjStagState::solveRiemannBoundary(Cell& cellLeft, const double& dxLeft, double& dtMax)
{
  model->solveRiemannInletInjStagState(cellLeft, dxLeft, dtMax, m_m0, m_ak0, m_rhok0, m_pk0, m_boundData);
}

//****************************************************************************

void BoundCondInletInjStagState::solveRiemannTransportBoundary(Cell& cellLeft) const
{
	model->solveRiemannTransportInletInjStagState(cellLeft, m_valueTransport);
}

//****************************************************************************

void BoundCondInletInjStagState::printInfo()
{
  std::cout << m_numPhysique << std::endl;
  std::cout << m_m0 << std::endl;
  std::cout << m_rhok0[0] << std::endl;
}

//***************************************************************************
//******************************AMR Method***********************************
//***************************************************************************

void BoundCondInletInjStagState::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondInletInjStagState(*this, m_lvl + 1));
}

//****************************************************************************