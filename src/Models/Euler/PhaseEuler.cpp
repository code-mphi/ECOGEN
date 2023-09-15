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

#include "PhaseEuler.h"

using namespace tinyxml2;

//***************************************************************************

PhaseEuler::PhaseEuler() : m_density(0.), m_pressure(0.), m_eos(0), m_energy(0.), m_totalEnergy(0.), m_soundSpeed(0.)
{
  m_velocity.setXYZ(0., 0., 0.);
}

//***************************************************************************

PhaseEuler::PhaseEuler(XMLElement* material, Eos* eos, std::string fileName) : m_eos(eos), m_energy(0.), m_totalEnergy(0.), m_soundSpeed(0.)
{
  XMLElement* sousElement(material->FirstChildElement("dataFluid"));
  if (sousElement == NULL) throw ErrorXMLElement("dataFluid", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;

  //Thermodynamic data reading
  int presenceDensity(0), presencePressure(0), presenceTemperature(0);
  if (sousElement->QueryDoubleAttribute("density", &m_density) == XML_NO_ERROR) presenceDensity = 1;
  if (sousElement->QueryDoubleAttribute("pressure", &m_pressure) == XML_NO_ERROR) presencePressure = 1;
  if (sousElement->QueryDoubleAttribute("temperature", &m_temperature) == XML_NO_ERROR) presenceTemperature = 1;

  //Attribute error gestion
  if (presenceDensity + presencePressure + presenceTemperature != 2) throw ErrorXMLAttribut("only two of following is required : density, pressure, temperature", fileName, __FILE__, __LINE__);

  //Thermodynamic reconstruction if needed
  if (presencePressure&&presenceTemperature) m_density = m_eos->computeDensity(m_pressure, m_temperature);
  if (presencePressure&&presenceDensity) m_temperature = m_eos->computeTemperature(m_density,m_pressure);
  if (presenceDensity&&presenceTemperature) throw ErrorXMLAttribut("impossible to initialize Euler phase with density and temperature", fileName, __FILE__, __LINE__);

  //velocity
  XMLElement* velocity(sousElement->FirstChildElement("velocity"));
  if (velocity == NULL) throw ErrorXMLElement("velocity", fileName, __FILE__, __LINE__);
  double velocityX(0.), velocityY(0.), velocityZ(0.);
  error = velocity->QueryDoubleAttribute("x", &velocityX);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName, __FILE__, __LINE__);
  error = velocity->QueryDoubleAttribute("y", &velocityY);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName, __FILE__, __LINE__);
  error = velocity->QueryDoubleAttribute("z", &velocityZ);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName, __FILE__, __LINE__);
  m_velocity.setXYZ(velocityX, velocityY, velocityZ);
}

//***************************************************************************

PhaseEuler::~PhaseEuler(){}

//***************************************************************************

void PhaseEuler::allocateAndCopyPhase(Phase** vecPhase)
{
  *vecPhase = new PhaseEuler(*this);
}

//***************************************************************************

void PhaseEuler::copyPhase(Phase &phase)
{
  m_density = phase.getDensity();
  m_pressure = phase.getPressure();
  m_velocity = phase.getVelocity();
  m_eos = phase.getEos();
  m_energy = phase.getEnergy();
  m_soundSpeed = phase.getSoundSpeed();
  m_totalEnergy = phase.getTotalEnergy();
}

//***************************************************************************

void PhaseEuler::extendedCalculusPhase(const Coord& /*velocity*/)
{
  m_energy = m_eos->computeEnergy(m_density, m_pressure);
  m_soundSpeed = m_eos->computeSoundSpeed(m_density, m_pressure);
  m_totalEnergy = m_energy + 0.5*m_velocity.squaredNorm();
  m_temperature = m_eos->computeTemperature(m_density, m_pressure);
}

//***************************************************************************

void PhaseEuler::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void PhaseEuler::reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhaseEuler::returnScalar(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_density; break;
  case 2:
    return m_pressure; break;
  case 3:
    return m_temperature; break;
  default:
    return 0.; break;
  }
}

//***************************************************************************

Coord PhaseEuler::returnVector(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_velocity; break;
  default:
    return 0; break;
  }
}

//***************************************************************************

std::string PhaseEuler::returnNameScalar(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Density"; break;
  case 2:
    return "Pressure"; break;
  case 3:
    return "Temperature"; break;
  default:
    return "NoName"; break;
  }
}

//***************************************************************************

std::string PhaseEuler::returnNameVector(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Velocity"; break;
  default:
    return "NoName"; break;
  }
}

//****************************************************************************
//************************* READING FROM FILE ********************************
//****************************************************************************

void PhaseEuler::setScalar(const int& numVar, const double& value)
{
  switch (numVar)
  {
  case 1:
    m_density = value; break;
  case 2:
    m_pressure = value; break;
  case 3:
    m_temperature = value; break;
  default:
    Errors::errorMessage("numVar not found in PhaseEuler::setScalar"); break;
  }
}

//****************************************************************************

void PhaseEuler::setVector(const int& numVar, const Coord& value)
{
  switch (numVar)
  {
  case 1:
    m_velocity = value; break;
  default:
    Errors::errorMessage("numVar not found in PhaseEuler::setVector"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseEuler::numberOfTransmittedVariables() const
{
  //5 variables + number EOS
  return 6;
}

//***************************************************************************

void PhaseEuler::fillBuffer(double* buffer, int& counter) const
{
  buffer[++counter] = m_density;
  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();
  buffer[++counter] = m_pressure;
  buffer[++counter] = static_cast<double>(m_eos->getNumber());
}

//***************************************************************************

void PhaseEuler::fillBuffer(std::vector<double>& dataToSend) const
{
  dataToSend.push_back(m_density);
  dataToSend.push_back(m_velocity.getX());
  dataToSend.push_back(m_velocity.getY());
  dataToSend.push_back(m_velocity.getZ());
  dataToSend.push_back(m_pressure);
  dataToSend.push_back(static_cast<double>(m_eos->getNumber()));
}

//***************************************************************************

void PhaseEuler::getBuffer(double* buffer, int& counter, Eos** eos)
{
  m_density = buffer[++counter];
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);
  m_pressure = buffer[++counter];
  m_eos = eos[static_cast<int>(buffer[++counter])];
}

//***************************************************************************

void PhaseEuler::getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** eos)
{
  m_density = dataToReceive[counter++];
  m_velocity.setX(dataToReceive[counter++]);
  m_velocity.setY(dataToReceive[counter++]);
  m_velocity.setZ(dataToReceive[counter++]);
  m_pressure = dataToReceive[counter++];
  m_eos = eos[static_cast<int>(dataToReceive[counter++])];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void PhaseEuler::computeSlopesPhase(const Phase &sLeft, const Phase &sRight, const double& distance)
{
  m_density = (sRight.getDensity() - sLeft.getDensity()) / distance;
  m_pressure = (sRight.getPressure() - sLeft.getPressure()) / distance;
  m_velocity.setX((sRight.getVelocity().getX() - sLeft.getVelocity().getX()) / distance);
  m_velocity.setY((sRight.getVelocity().getY() - sLeft.getVelocity().getY()) / distance);
  m_velocity.setZ((sRight.getVelocity().getZ() - sLeft.getVelocity().getZ()) / distance);
}

//***************************************************************************

void PhaseEuler::setToZero()
{
  m_density = 0.; m_pressure = 0.;
  m_velocity.setX(0.); m_velocity.setY(0.); m_velocity.setZ(0.);
}

//***************************************************************************

void PhaseEuler::setToMax()
{
  m_density = 1.e15; m_pressure = 1.e15;
  m_velocity.setX(1.e15); m_velocity.setY(1.e15); m_velocity.setZ(1.e15);
}

//***************************************************************************

void PhaseEuler::extrapolate(const Phase &slope, const double& distance)
{
  m_density += slope.getDensity() * distance;
  m_pressure += slope.getPressure() * distance;
  m_velocity.setX(m_velocity.getX() + slope.getVelocity().getX() * distance);
  m_velocity.setY(m_velocity.getY() + slope.getVelocity().getY() * distance);
  m_velocity.setZ(m_velocity.getZ() + slope.getVelocity().getZ() * distance);
}

//***************************************************************************

void PhaseEuler::limitSlopes(const Phase& slopeGauche, const Phase& slopeDroite, Limiter& globalLimiter, Limiter& /*volumeFractionLimiter*/)
{
  m_density = globalLimiter.limiteSlope(slopeGauche.getDensity(), slopeDroite.getDensity());
  m_pressure = globalLimiter.limiteSlope(slopeGauche.getPressure(), slopeDroite.getPressure());
  m_velocity.setX(globalLimiter.limiteSlope(slopeGauche.getVelocity().getX(), slopeDroite.getVelocity().getX()));
  m_velocity.setY(globalLimiter.limiteSlope(slopeGauche.getVelocity().getY(), slopeDroite.getVelocity().getY()));
  m_velocity.setZ(globalLimiter.limiteSlope(slopeGauche.getVelocity().getZ(), slopeDroite.getVelocity().getZ()));
}

//****************************************************************************

void PhaseEuler::setMin(const Phase& phase1, const Phase& phase2)
{
  m_density = std::min(phase1.getDensity(), phase2.getDensity());
  m_pressure = std::min(phase1.getPressure(), phase2.getPressure());
  
  m_velocity.setX(std::min(phase1.getVelocity().getX(), phase2.getVelocity().getX()));
  m_velocity.setY(std::min(phase1.getVelocity().getY(), phase2.getVelocity().getY()));
  m_velocity.setZ(std::min(phase1.getVelocity().getZ(), phase2.getVelocity().getZ()));
}

//****************************************************************************

void PhaseEuler::setMax(const Phase& phase1, const Phase& phase2)
{
  m_density = std::max(phase1.getDensity(), phase2.getDensity());
  m_pressure = std::max(phase1.getPressure(), phase2.getPressure());
  
  m_velocity.setX(std::max(phase1.getVelocity().getX(), phase2.getVelocity().getX()));
  m_velocity.setY(std::max(phase1.getVelocity().getY(), phase2.getVelocity().getY()));
  m_velocity.setZ(std::max(phase1.getVelocity().getZ(), phase2.getVelocity().getZ()));
}

//****************************************************************************

void PhaseEuler::computeGradientLimiter(const Limiter& globalLimiter, const Phase& phase, const Phase& phaseMin, const Phase& phaseMax, const Phase& slope)
{
  m_density = std::min(m_density, globalLimiter.computeGradientLimiter(phase.getDensity(), phaseMin.getDensity(), phaseMax.getDensity(), slope.getDensity()));
  m_pressure = std::min(m_pressure, globalLimiter.computeGradientLimiter(phase.getPressure(), phaseMin.getPressure(), phaseMax.getPressure(), slope.getPressure()));

  m_velocity.setX(std::min(m_velocity.getX(), globalLimiter.computeGradientLimiter(phase.getVelocity().getX(), phaseMin.getVelocity().getX(), phaseMax.getVelocity().getX(), slope.getVelocity().getX())));
  m_velocity.setY(std::min(m_velocity.getY(), globalLimiter.computeGradientLimiter(phase.getVelocity().getY(), phaseMin.getVelocity().getY(), phaseMax.getVelocity().getY(), slope.getVelocity().getY())));
  m_velocity.setZ(std::min(m_velocity.getZ(), globalLimiter.computeGradientLimiter(phase.getVelocity().getZ(), phaseMin.getVelocity().getZ(), phaseMax.getVelocity().getZ(), slope.getVelocity().getZ())));
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int PhaseEuler::numberOfTransmittedSlopes() const
{
	return 5;
}

//***************************************************************************

void PhaseEuler::fillBufferSlopes(double* buffer, int& counter) const
{
	buffer[++counter] = m_density;
	buffer[++counter] = m_velocity.getX();
	buffer[++counter] = m_velocity.getY();
	buffer[++counter] = m_velocity.getZ();
	buffer[++counter] = m_pressure;
}

//***************************************************************************

void PhaseEuler::getBufferSlopes(double* buffer, int& counter)
{
	m_density = buffer[++counter];
	m_velocity.setX(buffer[++counter]);
	m_velocity.setY(buffer[++counter]);
	m_velocity.setZ(buffer[++counter]);
	m_pressure = buffer[++counter];
}


//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseEuler::verifyPhase(const std::string& message) const
{
  if (m_density <= 1.e-10) errors.push_back(Errors(message + "too small density in verifyPhase"));
  m_eos->verifyPressure(m_pressure);
}

//***************************************************************************

void PhaseEuler::verifyAndCorrectPhase()
{
  if (m_density < 1.e-10) this->setDensity(1.e-10);
  m_eos->verifyAndModifyPressure(m_pressure);
}

//***************************************************************************

void PhaseEuler::verifyAndCorrectDensityMax(const double& /*mass*/)
{
  m_eos->verifyAndCorrectDensityMax(m_density);
}

//***************************************************************************

void PhaseEuler::verifyAndCorrectDensityMax()
{
  m_eos->verifyAndCorrectDensityMax(m_density);
}

//****************************************************************************
//**************************** DATA ACCESSORS ********************************
//****************************************************************************

void PhaseEuler::setDensity(double density) { m_density = density; }

//***************************************************************************

void PhaseEuler::setPressure(double pressure) { m_pressure = pressure; }

//***************************************************************************

void PhaseEuler::setVelocity(const double& u, const double& v, const double& w) { m_velocity.setXYZ(u, v, w); }

//***************************************************************************

void PhaseEuler::setVelocity(const Coord& vit) { m_velocity = vit; }

//***************************************************************************

void PhaseEuler::setU(const double& u) { m_velocity.setX(u); }

//***************************************************************************

void PhaseEuler::setV(const double& v) { m_velocity.setY(v); }

//***************************************************************************

void PhaseEuler::setW(const double& w) { m_velocity.setZ(w); }

//***************************************************************************

void PhaseEuler::setEos(Eos* eos) { m_eos = eos; }

//***************************************************************************

void PhaseEuler::setEnergy(double energy) { m_energy = energy; }

//***************************************************************************

void PhaseEuler::setSoundSpeed(double soundSpeed) { m_soundSpeed = soundSpeed; }

//***************************************************************************

void PhaseEuler::setTotalEnergy(double totalEnergy) { m_totalEnergy = totalEnergy; }

//***************************************************************************

void PhaseEuler::setTemperature(double temperature) { m_temperature = temperature; }

//****************************************************************************
//***************************** OPERATORS ************************************
//****************************************************************************

void PhaseEuler::changeSign()
{
  m_density = -m_density;
  m_pressure = -m_pressure;
  m_velocity = m_velocity*-1.;
}

//***************************************************************************

void PhaseEuler::multiplyAndAdd(const Phase &slopesPhasesTemp, const double& coeff)
{
  m_density += slopesPhasesTemp.getDensity()*coeff;
  m_pressure += slopesPhasesTemp.getPressure()*coeff;
  m_velocity += slopesPhasesTemp.getVelocity()*coeff;
}

//***************************************************************************

void PhaseEuler::divide(const double& coeff)
{
  m_density /= coeff;
  m_pressure /= coeff;
  m_velocity /= coeff;
}
