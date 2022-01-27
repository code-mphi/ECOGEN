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

#include "PhaseUEqTotEnergy.h"
#include "../../Eos/Eos.h"

using namespace tinyxml2;

//***************************************************************************

PhaseUEqTotEnergy::PhaseUEqTotEnergy() : m_alpha(1.0), m_density(0.), m_pressure(0.), m_eos(0), m_totEnergy(0.), m_soundSpeed(0.) {}

//***************************************************************************

PhaseUEqTotEnergy::PhaseUEqTotEnergy(XMLElement* material, Eos* eos, std::string fileName) : m_alpha(1.0), m_density(0.), m_pressure(0.), m_eos(eos), m_totEnergy(0.), m_soundSpeed(0.)
{
  XMLElement* sousElement(material->FirstChildElement("dataFluid"));
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
  if (presenceDensity&&presencePressure) m_temperature = m_eos->computeTemperature(m_density,m_pressure);
  if (presenceDensity&&presenceTemperature) throw ErrorXMLAttribut("impossible to initialize UEq phase with density and temperature", fileName, __FILE__, __LINE__);
}

//***************************************************************************

PhaseUEqTotEnergy::~PhaseUEqTotEnergy(){}

//***************************************************************************

void PhaseUEqTotEnergy::allocateAndCopyPhase(Phase** vecPhase)
{
  *vecPhase = new PhaseUEqTotEnergy(*this);
}

//***************************************************************************

void PhaseUEqTotEnergy::copyPhase(Phase &phase)
{
  m_alpha = phase.getAlpha();
  m_density = phase.getDensity();
  m_pressure = phase.getPressure();
  m_eos = phase.getEos();
  m_totEnergy = phase.getTotalEnergy();
  m_soundSpeed = phase.getSoundSpeed();
}

//***************************************************************************

void PhaseUEqTotEnergy::extendedCalculusPhase(const Coord& vel)
{
  m_temperature = m_eos->computeTemperature(m_density, m_pressure);
  m_totEnergy = m_eos->computeEnergy(m_density, m_pressure) + 0.5 * vel.squaredNorm();
  m_soundSpeed = m_eos->computeSoundSpeed(m_density, m_pressure);
}

//****************************************************************************

void PhaseUEqTotEnergy::computeMassFraction(const double& density)
{
  m_Y = m_alpha*m_density / std::max(density, epsilonAlphaNull);
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhaseUEqTotEnergy::returnScalar(const int& numVar) const
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
  case 5:
    return m_pressure;
    break;
  default:
    return 0.; break;
  }
}

//***************************************************************************

std::string PhaseUEqTotEnergy::returnNameScalar(const int& numVar) const
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

void PhaseUEqTotEnergy::setScalar(const int& numVar, const double& value)
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
    Errors::errorMessage("numVar not found in PhaseUEqTotEnergy::setScalar"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseUEqTotEnergy::numberOfTransmittedVariables() const
{
  //3 variables + number EOS
  return 4;
}

//***************************************************************************

void PhaseUEqTotEnergy::fillBuffer(double* buffer, int& counter) const
{
  buffer[++counter] = m_alpha;
  buffer[++counter] = m_density;
  buffer[++counter] = m_pressure;
  buffer[++counter] = static_cast<double>(m_eos->getNumber());
}

//***************************************************************************

void PhaseUEqTotEnergy::fillBuffer(std::vector<double>& dataToSend) const
{
  dataToSend.push_back(m_alpha);
  dataToSend.push_back(m_density);
  dataToSend.push_back(m_pressure);
  dataToSend.push_back(static_cast<double>(m_eos->getNumber()));
}

//***************************************************************************

void PhaseUEqTotEnergy::getBuffer(double* buffer, int& counter, Eos** eos)
{
  m_alpha = buffer[++counter];
  m_density = buffer[++counter];
  m_pressure = buffer[++counter];
  m_eos = eos[static_cast<int>(buffer[++counter])];
}

//***************************************************************************

void PhaseUEqTotEnergy::getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** eos)
{
  m_alpha = dataToReceive[counter++];
  m_density = dataToReceive[counter++];
  m_pressure = dataToReceive[counter++];
  m_eos = eos[static_cast<int>(std::round(dataToReceive[counter++]))];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void PhaseUEqTotEnergy::computeSlopesPhase(const Phase &sLeft, const Phase &sRight, const double& distance)
{
  m_alpha = (sRight.getAlpha() - sLeft.getAlpha()) / distance;
  m_density = (sRight.getDensity() - sLeft.getDensity()) / distance;
  m_pressure = (sRight.getPressure() - sLeft.getPressure()) / distance;
}

//***************************************************************************

void PhaseUEqTotEnergy::setToZero()
{
  m_alpha = 0.; m_density = 0.; m_pressure = 0.;
}

//***************************************************************************

void PhaseUEqTotEnergy::extrapolate(const Phase &slope, const double& distance)
{
  m_alpha += slope.getAlpha() * distance;
  m_density += slope.getDensity() * distance;
  m_pressure += slope.getPressure() * distance;
}

//***************************************************************************

void PhaseUEqTotEnergy::limitSlopes(const Phase &slopeGauche, const Phase &slopeDroite, Limiter& globalLimiter, Limiter& volumeFractionLimiter)
{
  m_alpha = volumeFractionLimiter.limiteSlope(slopeGauche.getAlpha(), slopeDroite.getAlpha());
  m_density = globalLimiter.limiteSlope(slopeGauche.getDensity(), slopeDroite.getDensity());
  m_pressure = globalLimiter.limiteSlope(slopeGauche.getPressure(), slopeDroite.getPressure());
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int PhaseUEqTotEnergy::numberOfTransmittedSlopes() const
{
	return 3;
}

//***************************************************************************

void PhaseUEqTotEnergy::fillBufferSlopes(double* buffer, int& counter) const
{
	buffer[++counter] = m_alpha;
	buffer[++counter] = m_density;
	buffer[++counter] = m_pressure;
}

//***************************************************************************

void PhaseUEqTotEnergy::getBufferSlopes(double* buffer, int& counter)
{
	m_alpha = buffer[++counter];
	m_density = buffer[++counter];
	m_pressure = buffer[++counter];
}

//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseUEqTotEnergy::verifyPhase(const std::string& message) const
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

void PhaseUEqTotEnergy::verifyAndCorrectPhase()
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

//***************************************************************************

void PhaseUEqTotEnergy::verifyAndCorrectDensityMax(const double& mass)
{
  m_eos->verifyAndCorrectDensityMax(mass, m_alpha, m_density);
}

//***************************************************************************

void PhaseUEqTotEnergy::verifyAndCorrectDensityMax()
{
  m_eos->verifyAndCorrectDensityMax(m_density);
}

//****************************************************************************
//**************************** DATA ACCESSORS ********************************
//****************************************************************************

void PhaseUEqTotEnergy::setAlpha(double alpha) { m_alpha = alpha; }

//***************************************************************************

void PhaseUEqTotEnergy::setDensity(double density) { m_density = density; }

//***************************************************************************

void PhaseUEqTotEnergy::setPressure(double pressure) { m_pressure = pressure; }

//***************************************************************************

void PhaseUEqTotEnergy::setEos(Eos* eos) { m_eos = eos; }

//***************************************************************************

void PhaseUEqTotEnergy::setTotalEnergy(double totalEnergy) { m_totEnergy = totalEnergy; }

//***************************************************************************

void PhaseUEqTotEnergy::setTotalEnergy(const double &energy, const Coord& vel)
{ 
  m_totEnergy = energy + 0.5 * vel.squaredNorm();
}

//***************************************************************************

void PhaseUEqTotEnergy::setSoundSpeed(double soundSpeed) { m_soundSpeed = soundSpeed; }

//****************************************************************************
//****************************** OPERATORS ***********************************
//****************************************************************************

void PhaseUEqTotEnergy::changeSign()
{
  m_alpha = -m_alpha;
  m_density = -m_density;
  m_pressure = -m_pressure;
}

//***************************************************************************

void PhaseUEqTotEnergy::multiplyAndAdd(const Phase &slopesPhasesTemp, const double& coeff)
{
  m_alpha += slopesPhasesTemp.getAlpha()*coeff;
  m_density += slopesPhasesTemp.getDensity()*coeff;
  m_pressure += slopesPhasesTemp.getPressure()*coeff;
}

//***************************************************************************

void PhaseUEqTotEnergy::divide(const double& coeff)
{
  m_alpha /= coeff;
  m_density /= coeff;
  m_pressure /= coeff;
}
