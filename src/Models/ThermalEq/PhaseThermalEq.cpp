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

//! \file      PhaseThermalEq.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "PhaseThermalEq.h"
#include "../../Eos/Eos.h"

using namespace tinyxml2;

//***************************************************************************

PhaseThermalEq::PhaseThermalEq() :m_alpha(1.0), m_density(0.), m_pressure(0.), m_eos(0), m_energie(0.), m_totalEnergy(0.), m_soundSpeed(0.) {}

//***************************************************************************

PhaseThermalEq::PhaseThermalEq(XMLElement *material, Eos *eos, std::string fileName) : m_eos(eos), m_energie(0.), m_totalEnergy(0.), m_soundSpeed(0.)
{
  XMLElement *sousElement(material->FirstChildElement("dataFluid"));
  if (sousElement == NULL) throw ErrorXMLElement("dataFluid", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  //alpha
  error = sousElement->QueryDoubleAttribute("alpha", &m_alpha);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("alpha", fileName, __FILE__, __LINE__);
}

//***************************************************************************

PhaseThermalEq::~PhaseThermalEq(){}

//***************************************************************************

void PhaseThermalEq::allocateAndCopyPhase(Phase **vecPhase)
{
  *vecPhase = new PhaseThermalEq(*this);
}

//***************************************************************************

void PhaseThermalEq::copyPhase(Phase &phase)
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

void PhaseThermalEq::extendedCalculusPhase(const Coord &velocity)
{
  m_energie = m_eos->computeEnergy(m_density, m_pressure);
  m_soundSpeed = m_eos->computeSoundSpeed(m_density, m_pressure);
  m_totalEnergy = m_energie + 0.5*velocity.squaredNorm();
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhaseThermalEq::returnScalar(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_alpha; break;
  case 2:
    return m_density; break;
  default:
    return 0.; break;
  }
}

//***************************************************************************

std::string PhaseThermalEq::returnNameScalar(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Alpha"; break;
  case 2:
    return "Density"; break;
  default:
    return "NoName"; break;
  }
}

//****************************************************************************
//************************* READING FROM FILE ********************************
//****************************************************************************

void PhaseThermalEq::setScalar(const int &numVar, const double &value)
{
  switch (numVar)
  {
  case 1:
    m_alpha = value; break;
  case 2:
    m_density = value; break;
  default:
    Errors::errorMessage("numVar not found in Phase::setScalar"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseThermalEq::numberOfTransmittedVariables() const
{
  //1 variable + number EOS
  return 2;
}

//***************************************************************************

void PhaseThermalEq::fillBuffer(double *buffer, int &counter) const
{
  buffer[++counter] = m_alpha;
  buffer[++counter] = static_cast<double>(m_eos->getNumber());
}

//***************************************************************************

void PhaseThermalEq::fillBuffer(std::vector<double> &dataToSend) const
{
  dataToSend.push_back(m_alpha);
  dataToSend.push_back(static_cast<double>(m_eos->getNumber()));
}

//***************************************************************************

void PhaseThermalEq::getBuffer(double *buffer, int &counter, Eos **eos)
{
  m_alpha = buffer[++counter];
  m_eos = eos[static_cast<int>(buffer[++counter])];
}

//***************************************************************************

void PhaseThermalEq::getBuffer(std::vector<double> &dataToReceive, int &counter, Eos **eos)
{
  m_alpha = dataToReceive[counter++];
  m_eos = eos[static_cast<int>(dataToReceive[counter++])];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void PhaseThermalEq::computeSlopesPhase(const Phase &sLeft, const Phase &sRight, const double &distance)
{
  m_alpha = (sRight.getAlpha() - sLeft.getAlpha()) / distance;
}

//***************************************************************************

void PhaseThermalEq::setToZero()
{
  m_alpha = 0.;
}

//***************************************************************************

void PhaseThermalEq::extrapolate(const Phase &slope, const double &distance)
{
  m_alpha += slope.getAlpha() * distance;
}

//***************************************************************************

void PhaseThermalEq::limitSlopes(const Phase &slopeGauche, const Phase &slopeDroite, Limiter &globalLimiter, Limiter &volumeFractionLimiter)
{
  m_alpha = volumeFractionLimiter.limiteSlope(slopeGauche.getAlpha(), slopeDroite.getAlpha());
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int PhaseThermalEq::numberOfTransmittedSlopes() const
{
	return 1;
}

//***************************************************************************

void PhaseThermalEq::fillBufferSlopes(double *buffer, int &counter) const
{
	buffer[++counter] = m_alpha;
}

//***************************************************************************

void PhaseThermalEq::getBufferSlopes(double *buffer, int &counter)
{
	m_alpha = buffer[++counter];
}

//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseThermalEq::verifyPhase(const std::string &message) const
{
  if (m_alpha <= 1e-10) errors.push_back(Errors(message + "too small alpha in verifyPhase"));
  if (m_density <= 1.e-10) errors.push_back(Errors(message + "too small density in verifyPhase"));
  m_eos->verifyPressure(m_pressure);
}

//***************************************************************************

void PhaseThermalEq::verifyAndCorrectPhase()
{
  if (m_alpha < 1e-10) this->setAlpha(1e-10);
  if (m_alpha > 1. - 1e-10) this->setAlpha(1. - 1e-10);
  if (m_density < 1.e-10) this->setDensity(1.e-10);
  m_eos->verifyAndModifyPressure(m_pressure);
}

//****************************************************************************
//**************************** DATA ACCESSORS ********************************
//****************************************************************************

void PhaseThermalEq::setAlpha(double alpha) { m_alpha = alpha; }

//***************************************************************************

void PhaseThermalEq::setDensity(double density) { m_density = density; }

//***************************************************************************

void PhaseThermalEq::setPressure(double pressure) { m_pressure = pressure; }

//***************************************************************************

void PhaseThermalEq::setEos(Eos *eos) { m_eos = eos; }

//***************************************************************************

void PhaseThermalEq::setEnergy(double energie) { m_energie = energie; }

//***************************************************************************

void PhaseThermalEq::setSoundSpeed(double soundSpeed) { m_soundSpeed = soundSpeed; }

//***************************************************************************

void PhaseThermalEq::setTotalEnergy(double totalEnergy) { m_totalEnergy = totalEnergy; }

//****************************************************************************
//****************************** OPERATORS ***********************************
//****************************************************************************

void PhaseThermalEq::changeSign()
{
  m_alpha = -m_alpha;
}

//***************************************************************************

void PhaseThermalEq::multiplyAndAdd(const Phase &slopesPhasesTemp, const double &coeff)
{
  m_alpha += slopesPhasesTemp.getAlpha()*coeff;
}

//***************************************************************************

void PhaseThermalEq::divide(const double &coeff)
{
  m_alpha /= coeff;
}
