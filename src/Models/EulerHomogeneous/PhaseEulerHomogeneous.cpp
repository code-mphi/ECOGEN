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

#include "PhaseEulerHomogeneous.h"
#include "../../Eos/Eos.h"

using namespace tinyxml2;

//***************************************************************************

PhaseEulerHomogeneous::PhaseEulerHomogeneous() :m_alpha(1.0), m_density(0.), m_pressure(0.), m_Y(0.), m_eos(0), m_energy(0.), m_totalEnergy(0.), m_soundSpeed(0.) {}

//***************************************************************************

PhaseEulerHomogeneous::PhaseEulerHomogeneous(XMLElement* material, Eos* eos, std::string fileName) : m_eos(eos), m_energy(0.), m_totalEnergy(0.), m_soundSpeed(0.)
{
  XMLElement* sousElement(material->FirstChildElement("dataFluid"));
  if (sousElement == NULL) throw ErrorXMLElement("dataFluid", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  //alpha
  error = sousElement->QueryDoubleAttribute("alpha", &m_alpha);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("alpha", fileName, __FILE__, __LINE__);
}

//***************************************************************************

PhaseEulerHomogeneous::~PhaseEulerHomogeneous() {}

//***************************************************************************

void PhaseEulerHomogeneous::allocateAndCopyPhase(Phase** vecPhase)
{
  *vecPhase = new PhaseEulerHomogeneous(*this);
}

//***************************************************************************

void PhaseEulerHomogeneous::copyPhase(Phase &phase)
{
  m_alpha = phase.getAlpha();
  m_density = phase.getDensity();
  m_pressure = phase.getPressure();
  m_eos = phase.getEos();
  m_energy = phase.getEnergy();
  m_soundSpeed = phase.getSoundSpeed();
  m_totalEnergy = phase.getTotalEnergy();
}

//***************************************************************************

void PhaseEulerHomogeneous::extendedCalculusPhase(const Coord& velocity)
{
  m_energy = m_eos->computeEnergy(m_density, m_pressure);
  m_soundSpeed = m_eos->computeSoundSpeed(m_density, m_pressure);
  m_totalEnergy = m_energy + 0.5*velocity.squaredNorm();
}

//****************************************************************************

void PhaseEulerHomogeneous::computeMassFraction(const double& density)
{
  m_Y = m_alpha*m_density / density;
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhaseEulerHomogeneous::returnScalar(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_alpha; break;
  case 2:
    return m_density; break;
  case 3:
    return m_Y; break;
  default:
    return 0.; break;
  }
}

//***************************************************************************

std::string PhaseEulerHomogeneous::returnNameScalar(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Alpha"; break;
  case 2:
    return "Density"; break;
  case 3:
    return "Mass Fraction"; break;
  default:
    return "NoName"; break;
  }
}

//****************************************************************************
//************************* READING FROM FILE ********************************
//****************************************************************************

void PhaseEulerHomogeneous::setScalar(const int& numVar, const double& value)
{
  switch (numVar)
  {
  case 1:
    m_alpha = value; break;
  case 2:
    m_density = value; break;
  case 3:
    m_Y = value; break;
  default:
    Errors::errorMessage("numVar not found in Phase::setScalar"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseEulerHomogeneous::numberOfTransmittedVariables() const
{
  //3 variables + number EOS
  return 4;
}

//***************************************************************************

void PhaseEulerHomogeneous::fillBuffer(double* buffer, int& counter) const
{
  buffer[++counter] = m_alpha;
  buffer[++counter] = m_density;
  buffer[++counter] = m_pressure;
  buffer[++counter] = static_cast<double>(m_eos->getNumber());
}

//***************************************************************************

void PhaseEulerHomogeneous::getBuffer(double* buffer, int& counter, Eos** eos)
{
  m_alpha = buffer[++counter];
  m_density = buffer[++counter];
  m_pressure = buffer[++counter];
  m_eos = eos[static_cast<int>(buffer[++counter])];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void PhaseEulerHomogeneous::computeSlopesPhase(const Phase &sLeft, const Phase &sRight, const double& distance)
{
  m_alpha = (sRight.getAlpha() - sLeft.getAlpha()) / distance;
}

//***************************************************************************

void PhaseEulerHomogeneous::setToZero()
{
  m_alpha = 0.;
}

//***************************************************************************

void PhaseEulerHomogeneous::extrapolate(const Phase &slope, const double& distance)
{
  m_alpha += slope.getAlpha() * distance;
}

//***************************************************************************

void PhaseEulerHomogeneous::limitSlopes(const Phase& slopeGauche, const Phase& slopeDroite, Limiter& /*globalLimiter*/, Limiter& volumeFractionLimiter)
{
  m_alpha = volumeFractionLimiter.limiteSlope(slopeGauche.getAlpha(), slopeDroite.getAlpha());
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int PhaseEulerHomogeneous::numberOfTransmittedSlopes() const
{
	return 1;
}

//***************************************************************************

void PhaseEulerHomogeneous::fillBufferSlopes(double* buffer, int& counter) const
{
	buffer[++counter] = m_alpha;
}

//***************************************************************************

void PhaseEulerHomogeneous::getBufferSlopes(double* buffer, int& counter)
{
	m_alpha = buffer[++counter];
}

//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseEulerHomogeneous::verifyPhase(const std::string& message) const
{
  if (m_alpha <= 1e-10) errors.push_back(Errors(message + "too small alpha in verifyPhase"));
  if (m_density <= 1.e-10) errors.push_back(Errors(message + "too small density in verifyPhase"));
  m_eos->verifyPressure(m_pressure);
}

//***************************************************************************

void PhaseEulerHomogeneous::verifyAndCorrectPhase()
{
  if (m_alpha < 1e-10) this->setAlpha(1e-10);
  if (m_alpha > 1. - 1e-10) this->setAlpha(1. - 1e-10);
  if (m_density < 1.e-10) this->setDensity(1.e-10);
  m_eos->verifyAndModifyPressure(m_pressure);
}

//***************************************************************************

void PhaseEulerHomogeneous::verifyAndCorrectDensityMax(const double& mass)
{
  m_eos->verifyAndCorrectDensityMax(mass, m_alpha, m_density);
}

//***************************************************************************

void PhaseEulerHomogeneous::verifyAndCorrectDensityMax()
{
  m_eos->verifyAndCorrectDensityMax(m_density);
}

//****************************************************************************
//**************************** DATA ACCESSORS ********************************
//****************************************************************************

void PhaseEulerHomogeneous::setAlpha(double alpha) { m_alpha = alpha; }

//***************************************************************************

void PhaseEulerHomogeneous::setDensity(double density) { m_density = density; }

//***************************************************************************

void PhaseEulerHomogeneous::setPressure(double pressure) { m_pressure = pressure; }

//***************************************************************************

void PhaseEulerHomogeneous::setEos(Eos* eos) { m_eos = eos; }

//***************************************************************************

void PhaseEulerHomogeneous::setEnergy(double energy) { m_energy = energy; }

//***************************************************************************

void PhaseEulerHomogeneous::setSoundSpeed(double soundSpeed) { m_soundSpeed = soundSpeed; }

//***************************************************************************

void PhaseEulerHomogeneous::setTotalEnergy(double totalEnergy) { m_totalEnergy = totalEnergy; }

//****************************************************************************
//***************************** OPERATORS ************************************
//****************************************************************************

void PhaseEulerHomogeneous::changeSign()
{
  m_alpha = -m_alpha;
}

//***************************************************************************

void PhaseEulerHomogeneous::multiplyAndAdd(const Phase &slopesPhasesTemp, const double& coeff)
{
  m_alpha += slopesPhasesTemp.getAlpha()*coeff;
}

//***************************************************************************

void PhaseEulerHomogeneous::divide(const double& coeff)
{
  m_alpha /= coeff;
}
