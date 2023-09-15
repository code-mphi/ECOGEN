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

#include "PhasePUEq.h"
#include "../../Eos/Eos.h"

using namespace tinyxml2;

//***************************************************************************

PhasePUEq::PhasePUEq() : PhaseUEq() {}

//***************************************************************************

PhasePUEq::PhasePUEq(XMLElement* material, Eos* eos, const double& pressure, std::string fileName) : PhaseUEq()
{
  m_pressure = pressure;
  m_eos = eos;

  XMLElement* sousElement(material->FirstChildElement("dataFluid"));
  if (sousElement == NULL) throw ErrorXMLElement("dataFluid", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  //alpha
  error = sousElement->QueryDoubleAttribute("alpha", &m_alpha);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("alpha", fileName, __FILE__, __LINE__);

  //Thermodynamic data reading
  int presenceDensity(0), presenceTemperature(0);
  if (sousElement->QueryDoubleAttribute("density", &m_density) == XML_NO_ERROR) presenceDensity = 1;
  if (sousElement->QueryDoubleAttribute("temperature", &m_temperature) == XML_NO_ERROR) presenceTemperature = 1;

  //Attribute error gestion
  if (presenceDensity + presenceTemperature == 2) throw ErrorXMLAttribut("only one of following is required : density, temperature", fileName, __FILE__, __LINE__);

  //Thermodynamic reconstruction if needed
  if (presenceTemperature) m_density = m_eos->computeDensity(m_pressure, m_temperature);
  if (presenceDensity) m_temperature = m_eos->computeTemperature(m_density,m_pressure);

  m_energy = m_eos->computeEnergy(m_density, m_pressure);
  m_soundSpeed = m_eos->computeSoundSpeed(m_density, m_pressure);
}

//***************************************************************************

PhasePUEq::~PhasePUEq(){}

//***************************************************************************

void PhasePUEq::allocateAndCopyPhase(Phase** vecPhase)
{
  *vecPhase = new PhasePUEq(*this);
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhasePUEq::returnScalar(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    if (m_alpha < 1.e-20) { return 0.; }
    else { return m_alpha; }
    break;
  case 2:
    if (m_density < 1.e-20) { return 1.e-20; }
    else { return m_density; }
    break;
  case 3:
    if (m_temperature < 1.e-20) { return 1.e-20; }
    else { return m_temperature; }
    break;
  case 4:
    if (m_Y < 1.e-20) { return 0.; }
    else { return m_Y; }
    break;
  default:
    return 0.; break;
  }
}

//***************************************************************************

std::string PhasePUEq::returnNameScalar(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Alpha"; break;
  case 2:
    return "Density"; break;
  case 3:
    return "Temperature"; break;
  case 4:
    return "Mass fraction"; break;
  default:
    return "NoName"; break;
  }
}

//****************************************************************************
//************************* READING FROM FILE ********************************
//****************************************************************************

void PhasePUEq::setScalar(const int& numVar, const double& value)
{
  switch (numVar)
  {
  case 1:
    m_alpha = value; break;
  case 2:
    m_density = value; break;
  case 3:
    m_temperature = value; break;
  case 4:
    m_Y = value; break;
  default:
    Errors::errorMessage("numVar not found in PhaseUEq::setScalar"); break;
  }
}

//***************************************************************************