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

//! \file      PhaseKapila.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "PhaseKapila.h"
#include "../../Eos/Eos.h"

using namespace tinyxml2;

//***************************************************************************

PhaseKapila::PhaseKapila() :m_alpha(1.0), m_density(0.), m_pressure(0.), m_eos(0), m_energie(0.), m_soundSpeed(0.) {}

//***************************************************************************

PhaseKapila::PhaseKapila(XMLElement *material, Eos *eos, const double &pressure, std::string fileName) : m_eos(eos), m_energie(0.), m_soundSpeed(0.), m_pressure(pressure)
{
  XMLElement *sousElement(material->FirstChildElement("dataFluid"));
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
}

//***************************************************************************

PhaseKapila::~PhaseKapila(){}

//***************************************************************************

void PhaseKapila::allocateAndCopyPhase(Phase **vecPhase)
{
  *vecPhase = new PhaseKapila(*this);
}

//***************************************************************************

void PhaseKapila::copyPhase(Phase &phase)
{
  m_alpha = phase.getAlpha();
  m_density = phase.getDensity();
  m_pressure = phase.getPressure();
  m_eos = phase.getEos();
  m_energie = phase.getEnergy();
  m_soundSpeed = phase.getSoundSpeed();
}

//***************************************************************************

void PhaseKapila::extendedCalculusPhase(const Coord &velocity)
{
  m_temperature = m_eos->computeTemperature(m_density, m_pressure);
  m_energie = m_eos->computeEnergy(m_density, m_pressure);
  m_soundSpeed = m_eos->computeSoundSpeed(m_density, m_pressure);
}

//****************************************************************************

void PhaseKapila::computeMassFraction(const double &density)
{
  m_Y = m_alpha*m_density / std::max(density, epsilonAlphaNull);
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhaseKapila::returnScalar(const int &numVar) const
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

std::string PhaseKapila::returnNameScalar(const int &numVar) const
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

void PhaseKapila::setScalar(const int &numVar, const double &value)
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
    Errors::errorMessage("numVar not found in PhaseKapila::setScalar"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseKapila::numberOfTransmittedVariables() const
{
  //3 variables + number EOS
  return 4;
}

//***************************************************************************

void PhaseKapila::fillBuffer(double *buffer, int &counter) const
{
  buffer[++counter] = m_alpha;
  buffer[++counter] = m_density;
  buffer[++counter] = m_pressure;
  buffer[++counter] = static_cast<double>(m_eos->getNumber());
}

//***************************************************************************

void PhaseKapila::fillBuffer(std::vector<double> &dataToSend) const
{
  dataToSend.push_back(m_alpha);
  dataToSend.push_back(m_density);
  dataToSend.push_back(m_pressure);
  dataToSend.push_back(static_cast<double>(m_eos->getNumber()));
}

//***************************************************************************

void PhaseKapila::getBuffer(double *buffer, int &counter, Eos **eos)
{
  m_alpha = buffer[++counter];
  m_density = buffer[++counter];
  m_pressure = buffer[++counter];
  m_eos = eos[static_cast<int>(buffer[++counter])];
}

//***************************************************************************

void PhaseKapila::getBuffer(std::vector<double> &dataToReceive, int &counter, Eos **eos)
{
  m_alpha = dataToReceive[counter++];
  m_density = dataToReceive[counter++];
  m_pressure = dataToReceive[counter++];
  m_eos = eos[static_cast<int>(std::round(dataToReceive[counter++]))];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void PhaseKapila::computeSlopesPhase(const Phase &sLeft, const Phase &sRight, const double &distance)
{
  m_alpha = (sRight.getAlpha() - sLeft.getAlpha()) / distance;
  m_density = (sRight.getDensity() - sLeft.getDensity()) / distance;
  m_pressure = (sRight.getPressure() - sLeft.getPressure()) / distance;
}

//***************************************************************************

void PhaseKapila::setToZero()
{
  m_alpha = 0.; m_density = 0.; m_pressure = 0.;
}

//***************************************************************************

void PhaseKapila::extrapolate(const Phase &slope, const double &distance)
{
  m_alpha += slope.getAlpha() * distance;
  m_density += slope.getDensity() * distance;
  m_pressure += slope.getPressure() * distance;
}

//***************************************************************************

void PhaseKapila::limitSlopes(const Phase &slopeGauche, const Phase &slopeDroite, Limiter &globalLimiter, Limiter &volumeFractionLimiter)
{
  m_alpha = volumeFractionLimiter.limiteSlope(slopeGauche.getAlpha(), slopeDroite.getAlpha());
  m_density = globalLimiter.limiteSlope(slopeGauche.getDensity(), slopeDroite.getDensity());
  m_pressure = globalLimiter.limiteSlope(slopeGauche.getPressure(), slopeDroite.getPressure());
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int PhaseKapila::numberOfTransmittedSlopes() const
{
	return 3;
}

//***************************************************************************

void PhaseKapila::fillBufferSlopes(double *buffer, int &counter) const
{
	buffer[++counter] = m_alpha;
	buffer[++counter] = m_density;
	buffer[++counter] = m_pressure;
}

//***************************************************************************

void PhaseKapila::getBufferSlopes(double *buffer, int &counter)
{
	m_alpha = buffer[++counter];
	m_density = buffer[++counter];
	m_pressure = buffer[++counter];
}

//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseKapila::verifyPhase(const std::string &message) const
{
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    if (m_alpha < 0.) errors.push_back(Errors(message + "too small alpha in verifyPhase"));
    if (m_alpha > 1.) errors.push_back(Errors(message + "too big alpha in verifyPhase"));
    if (m_density < 0.) errors.push_back(Errors(message + "too small density in verifyPhase"));
    m_eos->verifyPressure(m_pressure, message);
  }
  else { // alpha = 0 is desactivated (alpha != 0)
    if (m_alpha <= 1e-15) errors.push_back(Errors(message + "too small alpha in verifyPhase"));
    if (m_alpha >= 1. - 1e-15) errors.push_back(Errors(message + "too big alpha in verifyPhase"));
    if (m_density <= 1.e-15) errors.push_back(Errors(message + "too small density in verifyPhase"));
  }
}

//***************************************************************************

void PhaseKapila::verifyAndCorrectPhase()
{
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    if (m_alpha < 0.) m_alpha = 0.;
    if (m_alpha > 1.) m_alpha = 1.;
    if (m_density <= 1.e-15) m_density = 1.e-15;
  }
  else { // alpha = 0 is desactivated (alpha != 0)
    if (m_alpha < 1e-15) m_alpha = 1e-14;
    if (m_alpha > 1. - 1e-15) m_alpha = 1. - 1e-14;
    if (m_density < 1.e-15) m_density = 1.e-14;
  }
  m_eos->verifyAndModifyPressure(m_pressure);
}

//****************************************************************************
//**************************** DATA ACCESSORS ********************************
//****************************************************************************

void PhaseKapila::setAlpha(double alpha) { m_alpha = alpha; }

//***************************************************************************

void PhaseKapila::setDensity(double density) { m_density = density; }

//***************************************************************************

void PhaseKapila::setPressure(double pressure) { m_pressure = pressure; }

//***************************************************************************

void PhaseKapila::setEos(Eos *eos) { m_eos = eos; }

//***************************************************************************

void PhaseKapila::setEnergy(double energie) { m_energie = energie; }

//***************************************************************************

void PhaseKapila::setSoundSpeed(double soundSpeed) { m_soundSpeed = soundSpeed; }

//****************************************************************************
//****************************** OPERATORS ***********************************
//****************************************************************************

void PhaseKapila::changeSign()
{
  m_alpha = -m_alpha;
  m_density = -m_density;
  m_pressure = -m_pressure;
}

//***************************************************************************

void PhaseKapila::multiplyAndAdd(const Phase &slopesPhasesTemp, const double &coeff)
{
  m_alpha += slopesPhasesTemp.getAlpha()*coeff;
  m_density += slopesPhasesTemp.getDensity()*coeff;
  m_pressure += slopesPhasesTemp.getPressure()*coeff;
}

//***************************************************************************

void PhaseKapila::divide(const double &coeff)
{
  m_alpha /= coeff;
  m_density /= coeff;
  m_pressure /= coeff;
}
