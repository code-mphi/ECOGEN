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

#include "BoundCondInletTank.h"

using namespace tinyxml2;

//****************************************************************************

BoundCondInletTank::BoundCondInletTank(int numPhysique, XMLElement* element, const int& numbPhases, const int& numbTransports, std::vector<std::string> nameTransports, Eos** eos, std::string fileName) :
  BoundCond(numPhysique)
{
  m_ak0 = new double[numbPhases];
  m_Yk0 = new double[numbPhases];
  m_rhok0 = new double[numbPhases]; 
  
  //Reading tank pressure and temperature conditions
  //------------------------------------------------
  XMLElement* sousElement(element->FirstChildElement("dataInletTank"));
  if (sousElement == NULL) throw ErrorXMLElement("dataInletTank", fileName, __FILE__, __LINE__);
  //Attributes reading
  XMLError error;
  error = sousElement->QueryDoubleAttribute("p0", &m_p0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("p0", fileName, __FILE__, __LINE__);
  error = sousElement->QueryDoubleAttribute("T0", &m_T0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("T0", fileName, __FILE__, __LINE__);

  if (numbPhases == 1) {
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
      for (e = 0; e < numbPhases; e++) {
        if (nameEOS == eos[e]->getName()) { break; }
      }
      if (e == numbPhases) { throw ErrorXMLEOSUnknown(nameEOS, fileName, __FILE__, __LINE__); }

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
    if (nbFluids != numbPhases) throw ErrorXMLEtat("Tank", fileName, __FILE__, __LINE__);

    //Proportions checking
    //--------------------
    double sum(0.);
    if (presenceAlpha) {
      for (int k = 0; k < numbPhases; k++) {
        if (m_ak0[k]<0. || m_ak0[k]>1.) throw ErrorXMLAttribut("alpha should be in [0,1]", fileName, __FILE__, __LINE__);
        sum += m_ak0[k];
      }
      if (std::fabs(sum - 1.) > 1.e-6) { throw ErrorXMLAttribut("sum of alpha should be 1", fileName, __FILE__, __LINE__); }
      else {
        for (int k = 0; k < numbPhases; k++) { m_ak0[k] /= sum; }
      }
    }
    else if (presenceMassFrac) {
      for (int k = 0; k < numbPhases; k++) {
        if (m_Yk0[k]<0. || m_Yk0[k]>1.) throw ErrorXMLAttribut("massFrac should be in [0,1]", fileName, __FILE__, __LINE__);
        sum += m_Yk0[k];
      }
      if (std::fabs(sum - 1.) > 1.e-6) { throw ErrorXMLAttribut("sum of massFrac should be 1", fileName, __FILE__, __LINE__); }
      else {
        for (int k = 0; k < numbPhases; k++) { m_Yk0[k] /= sum; }
      }
    }
    else { throw ErrorXMLAttribut("One of following is required : alpha, massFrac", fileName, __FILE__, __LINE__); }

    //Fulfill tank state (rhok0, ak0 or Yk0)
    //--------------------------------------
    for (int k = 0; k < numbPhases; k++) {
      m_rhok0[k] = eos[k]->computeDensity(m_p0, m_T0);
    }
    double rhoMel(0.);
    if (presenceAlpha) {
      for (int k = 0; k < numbPhases; k++) { rhoMel += m_ak0[k] * m_rhok0[k]; }
      for (int k = 0; k < numbPhases; k++) { m_Yk0[k] = m_ak0[k] * m_rhok0[k] / rhoMel; }
    }
    else {
      for (int k = 0; k < numbPhases; k++) { rhoMel += m_Yk0[k] / m_rhok0[k]; }
      rhoMel = 1.0 / rhoMel;
      for (int k = 0; k < numbPhases; k++) { m_ak0[k] = rhoMel * m_Yk0[k] / m_rhok0[k]; }
    }

  } //End proportion

  //Reading of transports
  //---------------------
  m_valueTransport = new double[numbTransports];
  if (numbTransports) {
    XMLElement* sousElement(element->FirstChildElement("dataInletTank"));
    if (sousElement == NULL) throw ErrorXMLElement("dataInletTank", fileName, __FILE__, __LINE__);
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
    if (numbTransports > foundColors) throw ErrorXMLAttribut("Not enough transport equations in tank BC", fileName, __FILE__, __LINE__);
  }
}

//****************************************************************************

BoundCondInletTank::BoundCondInletTank(const BoundCondInletTank &Source, const int& lvl) : BoundCond(Source, lvl)
{
  m_ak0 = new double[numberPhases];
  m_Yk0 = new double[numberPhases];
  m_rhok0 = new double[numberPhases];

  for (int k = 0; k < numberPhases; k++)
  {
    m_ak0[k] = Source.m_ak0[k];
    m_Yk0[k] = Source.m_Yk0[k];
    m_rhok0[k] = Source.m_rhok0[k];
  }
  m_p0 = Source.m_p0;
  m_T0 = Source.m_T0;

  m_valueTransport = new double[numberTransports];
  for (int k = 0; k < numberTransports; k++) {
    m_valueTransport[k] = Source.m_valueTransport[k];
  }
}

//****************************************************************************

BoundCondInletTank::~BoundCondInletTank()
{
  delete[] m_ak0;
  delete[] m_Yk0;
  delete[] m_rhok0;
  delete[] m_valueTransport;
}

//****************************************************************************

void BoundCondInletTank::createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondInletTank(*(this)));
}

//****************************************************************************

void BoundCondInletTank::solveRiemannBoundary(Cell& cellLeft, const double& dxLeft, double& dtMax)
{
  model->solveRiemannInletTank(cellLeft, dxLeft, dtMax, m_ak0, m_rhok0, m_p0, m_T0, m_boundData);
}

//****************************************************************************

void BoundCondInletTank::solveRiemannTransportBoundary(Cell& cellLeft) const
{
	model->solveRiemannTransportInletTank(cellLeft, m_valueTransport);
}

//****************************************************************************

void BoundCondInletTank::printInfo()
{
  std::cout << m_numPhysique << std::endl;
  std::cout << m_rhok0[0] << std::endl;
}

//****************************************************************************
//******************************AMR Method************************************
//****************************************************************************

void BoundCondInletTank::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondInletTank(*this, m_lvl + 1));
}

//****************************************************************************
