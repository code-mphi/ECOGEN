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

//! \file      PhaseMultiP.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      June 5 2017

#include "PhaseMultiP.h"
#include "../../Eos/Eos.h"
#include <fstream>

using namespace std;
using namespace tinyxml2;

//***************************************************************************

PhaseMultiP::PhaseMultiP() :m_alpha(1.0), m_density(0.), m_pressure(0.), m_eos(0), m_energie(0.), m_totalEnergy(0.), m_soundSpeed(0.) {}

//***************************************************************************

PhaseMultiP::PhaseMultiP(XMLElement *material, Eos *eos, string fileName) : m_eos(eos), m_energie(0.), m_totalEnergy(0.), m_soundSpeed(0.)
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
  int presenceDensity(0), presencePressure(0), presenceTemperature(0);
  if (sousElement->QueryDoubleAttribute("density", &m_density) == XML_NO_ERROR) presenceDensity = 1;
  if (sousElement->QueryDoubleAttribute("pressure", &m_pressure) == XML_NO_ERROR) presencePressure = 1;
  if (sousElement->QueryDoubleAttribute("temperature", &m_temperature) == XML_NO_ERROR) presenceTemperature = 1;

  //Attribute error gestion
  if (presenceDensity + presencePressure + presenceTemperature != 2) throw ErrorXMLAttribut("only two of following is required : density, pressure, temperature", fileName, __FILE__, __LINE__);

  //Thermodynamic reconstruction if needed
  if (presencePressure&&presenceTemperature) m_density = m_eos->computeDensity(m_pressure, m_temperature);
  if (presenceDensity&&presenceTemperature) throw ErrorXMLAttribut("impossible to initialize MultiP phase with density and temperature", fileName, __FILE__, __LINE__);
}

//***************************************************************************

PhaseMultiP::~PhaseMultiP(){}

//***************************************************************************

void PhaseMultiP::allocateAndCopyPhase(Phase **vecPhase)
{
  *vecPhase = new PhaseMultiP(*this);
}

//***************************************************************************

void PhaseMultiP::copyPhase(Phase &phase)
{
  m_alpha = phase.getAlpha();
  m_density = phase.getDensity();
  m_pressure = phase.getPressure();
  m_eos = phase.getEos();
  m_energie = phase.getEnergy();
  m_soundSpeed = phase.getSoundSpeed();
  m_totalEnergy = phase.getTotalEnergy();
}

//***************************************************************************

void PhaseMultiP::extendedCalculusPhase(const Coord &velocity)
{
  m_temperature = m_eos->computeTemperature(m_density, m_pressure);
  m_energie = m_eos->computeEnergy(m_density, m_pressure);
  m_soundSpeed = m_eos->computeSoundSpeed(m_density, m_pressure);
  m_totalEnergy = m_energie + 0.5*velocity.squaredNorm();
}

//****************************************************************************

void PhaseMultiP::computeMassFraction(const double &density)
{
  m_Y = m_alpha*m_density / density;
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhaseMultiP::returnScalar(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_alpha; break;
  case 2:
    return m_density; break;
  case 3:
    return m_temperature; break;
  case 4:
    return m_Y; break;
  case 5:
    return m_pressure; break;
  default:
    return 0.; break;
  }
}

//***************************************************************************

string PhaseMultiP::returnNameScalar(const int &numVar) const
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
  case 5:
    return "Pressure"; break;
  default:
    return "NoName"; break;
  }
}

//****************************************************************************
//************************* READING FROM FILE ********************************
//****************************************************************************

void PhaseMultiP::setScalar(const int &numVar, const double &value)
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
  case 5:
    m_pressure = value; break;
  default:
    Errors::errorMessage("numVar not found in PhaseMultiP::setScalar"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseMultiP::numberOfTransmittedVariables() const
{
  //3 variables + number EOS
  return 4;
}

//***************************************************************************

void PhaseMultiP::fillBuffer(double *buffer, int &counter) const
{
  buffer[++counter] = m_alpha;
  buffer[++counter] = m_density;
  buffer[++counter] = m_pressure;
  buffer[++counter] = static_cast<double>(m_eos->getNumber());
}

//***************************************************************************

void PhaseMultiP::getBuffer(double *buffer, int &counter, Eos **eos)
{
  m_alpha = buffer[++counter];
  m_density = buffer[++counter];
  m_pressure = buffer[++counter];
  m_eos = eos[static_cast<int>(buffer[++counter])];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void PhaseMultiP::computeSlopesPhase(const Phase &sLeft, const Phase &sRight, const double &distance)
{
  m_alpha = (sRight.getAlpha() - sLeft.getAlpha()) / distance;
  m_density = (sRight.getDensity() - sLeft.getDensity()) / distance;
  m_pressure = (sRight.getPressure() - sLeft.getPressure()) / distance;
}

//***************************************************************************

void PhaseMultiP::setToZero()
{
  m_alpha = 0.; m_density = 0.; m_pressure = 0.;
}

//***************************************************************************

void PhaseMultiP::extrapolate(const Phase &slope, const double &distance)
{
  m_alpha += slope.getAlpha() * distance;
  m_density += slope.getDensity() * distance;
  m_pressure += slope.getPressure() * distance;
}

//***************************************************************************

void PhaseMultiP::limitSlopes(const Phase &slopeGauche, const Phase &slopeDroite, Limiter &globalLimiter, Limiter &volumeFractionLimiter)
{
  m_alpha = volumeFractionLimiter.limiteSlope(slopeGauche.getAlpha(), slopeDroite.getAlpha());
  m_density = globalLimiter.limiteSlope(slopeGauche.getDensity(), slopeDroite.getDensity());
  m_pressure = globalLimiter.limiteSlope(slopeGauche.getPressure(), slopeDroite.getPressure());
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int PhaseMultiP::numberOfTransmittedSlopes() const
{
	return 3;
}

//***************************************************************************

void PhaseMultiP::fillBufferSlopes(double *buffer, int &counter) const
{
	buffer[++counter] = m_alpha;
	buffer[++counter] = m_density;
	buffer[++counter] = m_pressure;
}

//***************************************************************************

void PhaseMultiP::getBufferSlopes(double *buffer, int &counter)
{
	m_alpha = buffer[++counter];
	m_density = buffer[++counter];
	m_pressure = buffer[++counter];
}

//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseMultiP::verifyPhase(const string &message) const
{
  if (m_alpha <= 1e-10) errors.push_back(Errors(message + "too small alpha in verifyPhase"));
  if (m_alpha >= 1. - 1e-10) errors.push_back(Errors(message + "too big alpha in verifyPhase"));
  if (m_density <= 1.e-10) errors.push_back(Errors(message + "too small density in verifyPhase"));
  m_eos->verifyPressure(m_pressure, message);
}

//***************************************************************************

void PhaseMultiP::verifyAndCorrectPhase()
{
  if (m_alpha < 1e-10) m_alpha = 1e-9;
  if (m_alpha > 1. - 1e-10) m_alpha = 1. - 1e-9;
  if (m_density < 1.e-10) m_density = 1.e-9;
  m_eos->verifyAndModifyPressure(m_pressure);
}

//****************************************************************************
//**************************** DATA ACCESSORS ********************************
//****************************************************************************

double PhaseMultiP::getAlpha() const { return m_alpha; }

//***************************************************************************

double PhaseMultiP::getDensity() const { return m_density; }

//***************************************************************************

double PhaseMultiP::getPressure() const { return m_pressure; }

//***************************************************************************

double PhaseMultiP::getY() const { return m_Y; }

//***************************************************************************

Eos* PhaseMultiP::getEos() const { return m_eos; }

//***************************************************************************

double PhaseMultiP::getEnergy() const { return m_energie; }

//***************************************************************************

double PhaseMultiP::getSoundSpeed() const { return m_soundSpeed; }

//***************************************************************************

double PhaseMultiP::getTotalEnergy() const { return m_totalEnergy; }

//***************************************************************************

double PhaseMultiP::getTemperature() const { return m_eos->computeTemperature(m_density, m_pressure); }

//***************************************************************************

void PhaseMultiP::setAlpha(double alpha) { m_alpha = alpha; }

//***************************************************************************

void PhaseMultiP::setDensity(double density) { m_density = density; }

//***************************************************************************

void PhaseMultiP::setPressure(double pressure) { m_pressure = pressure; }

//***************************************************************************

void PhaseMultiP::setEos(Eos *eos) { m_eos = eos; }

//***************************************************************************

void PhaseMultiP::setEnergy(double energie) { m_energie = energie; }

//***************************************************************************

void PhaseMultiP::setSoundSpeed(double soundSpeed) { m_soundSpeed = soundSpeed; }

//***************************************************************************

void PhaseMultiP::setTotalEnergy(double totalEnergy) { m_totalEnergy = totalEnergy; }

//****************************************************************************
//****************************** OPERATORS ***********************************
//****************************************************************************

void PhaseMultiP::changeSign()
{
  m_alpha = -m_alpha;
  m_density = -m_density;
  m_pressure = -m_pressure;
}

//***************************************************************************

void PhaseMultiP::multiplyAndAdd(const Phase &slopesPhasesTemp, const double &coeff)
{
  m_alpha += slopesPhasesTemp.getAlpha()*coeff;
  m_density += slopesPhasesTemp.getDensity()*coeff;
  m_pressure += slopesPhasesTemp.getPressure()*coeff;
}

//***************************************************************************

void PhaseMultiP::divide(const double &coeff)
{
  m_alpha /= coeff;
  m_density /= coeff;
  m_pressure /= coeff;
}
