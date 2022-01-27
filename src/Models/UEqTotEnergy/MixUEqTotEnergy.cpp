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

#include <cmath>
#include "MixUEqTotEnergy.h"

using namespace tinyxml2;

//***************************************************************************

MixUEqTotEnergy::MixUEqTotEnergy() : m_density(0.), m_pressure(0.), m_velocity(0), m_frozenSoundSpeed(0.), m_woodSoundSpeed(0.) {}

//***************************************************************************

MixUEqTotEnergy::MixUEqTotEnergy(XMLElement* state, std::string fileName) :
  m_density(0.), m_pressure(0.), m_velocity(0), m_frozenSoundSpeed(0.), m_woodSoundSpeed(0.)
{
  XMLElement* sousElement(state->FirstChildElement("mixture"));
  if (sousElement == NULL) throw ErrorXMLElement("mixture", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
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

MixUEqTotEnergy::~MixUEqTotEnergy(){}

//***************************************************************************

void MixUEqTotEnergy::allocateAndCopyMixture(Mixture** mixture)
{
  *mixture = new MixUEqTotEnergy(*this);
}

//***************************************************************************

void MixUEqTotEnergy::copyMixture(Mixture &mixture)
{
  m_density = mixture.getDensity();
  m_pressure = mixture.getPressure();
  m_velocity = mixture.getVelocity();
  m_frozenSoundSpeed = mixture.getFrozenSoundSpeed();
  m_woodSoundSpeed = mixture.getWoodSoundSpeed();
}

//***************************************************************************

double MixUEqTotEnergy::computeDensity(const double* alphak, const double* rhok)
{
  double rho(0.);
  for(int k=0;k<numberPhases;k++)
  {
      rho += alphak[k]*rhok[k];
  }
  return rho;
}

//***************************************************************************

double MixUEqTotEnergy::computePressure(const double* alphak, const double* pk)
{
  double p(0.);
  for(int k=0;k<numberPhases;k++)
  {
      p += alphak[k]*pk[k];
  }
  return p;
}

//***************************************************************************

double MixUEqTotEnergy::computeInternalEnergy(const double* Yk, const double* ek)
{
  double e(0.);
  for(int k=0;k<numberPhases;k++)
  {
      e += Yk[k]*ek[k];
  }
  return e;
}

//***************************************************************************

double MixUEqTotEnergy::computeFrozenSoundSpeed(const double* Yk, const double* ck)
{
  double cF(0.);
  for(int k=0;k<numberPhases;k++)
  {
      cF += Yk[k]*ck[k]*ck[k];
  }
  return sqrt(cF);
}

//***************************************************************************

void MixUEqTotEnergy::computeMixtureVariables(Phase** vecPhase)
{
  //mixture density and pressure
  m_density = 0.;
  m_pressure = 0.;
  for (int k = 0; k < numberPhases; k++) {
    m_density += vecPhase[k]->getAlpha()*vecPhase[k]->getDensity();
    m_pressure += vecPhase[k]->getAlpha()*vecPhase[k]->getPressure();
  }
  //Mass fraction
  for (int k = 0; k < numberPhases; k++) {
    vecPhase[k]->computeMassFraction(m_density);
  }
  //Speed of sounds
  m_frozenSoundSpeed = 0.;
  m_woodSoundSpeed = 0.;
  for (int k = 0; k < numberPhases; k++) {
    m_frozenSoundSpeed += vecPhase[k]->getY() * vecPhase[k]->getSoundSpeed()*vecPhase[k]->getSoundSpeed();
    m_woodSoundSpeed += vecPhase[k]->getAlpha() / std::max((vecPhase[k]->getDensity()*vecPhase[k]->getSoundSpeed()*vecPhase[k]->getSoundSpeed()), epsilonAlphaNull);
  }
  m_frozenSoundSpeed = sqrt(m_frozenSoundSpeed);
  m_woodSoundSpeed = 1. / sqrt(m_density*m_woodSoundSpeed);
}

//***************************************************************************

void MixUEqTotEnergy::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void MixUEqTotEnergy::reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
//**************************** DATA PRINTING *********************************
//****************************************************************************

double MixUEqTotEnergy::returnScalar(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_density; break;
  case 2:
    return m_pressure; break;
  default:
    return 0.; break;
  }
}

//***************************************************************************

Coord MixUEqTotEnergy::returnVector(const int& numVar) const
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

std::string MixUEqTotEnergy::returnNameScalar(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Density_Mixture"; break;
  case 2:
    return "Pressure_Mixture"; break;
  default:
    return "NoName"; break;
  }
}

//***************************************************************************

std::string MixUEqTotEnergy::returnNameVector(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Velocity_Mixture"; break;
  default:
    return "NoName"; break;
  }
}

//****************************************************************************
//**************************** DATA READING **********************************
//****************************************************************************

void MixUEqTotEnergy::setScalar(const int& numVar, const double& value)
{
  switch (numVar)
  {
  case 1:
    m_density = value; break;
  case 2:
    m_pressure = value; break;
  default:
    Errors::errorMessage("numVar not found in MixUEqTotEnergy::setScalar"); break;
  }
}

//***************************************************************************

void MixUEqTotEnergy::setVector(const int& numVar, const Coord& value)
{
  switch (numVar)
  {
  case 1:
    m_velocity = value; break;
  default:
    Errors::errorMessage("numVar not found in MixUEqTotEnergy::setVector"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int MixUEqTotEnergy::numberOfTransmittedVariables() const
{
  //1 vector : 3 variables
  return 3;
}

//***************************************************************************

void MixUEqTotEnergy::fillBuffer(double* buffer, int& counter) const
{
  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();
}

//***************************************************************************

void MixUEqTotEnergy::fillBuffer(std::vector<double>& dataToSend) const
{
  dataToSend.push_back(m_velocity.getX());
  dataToSend.push_back(m_velocity.getY());
  dataToSend.push_back(m_velocity.getZ());
}

//***************************************************************************

void MixUEqTotEnergy::getBuffer(double* buffer, int& counter)
{
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);
}

//***************************************************************************

void MixUEqTotEnergy::getBuffer(std::vector<double>& dataToReceive, int& counter)
{
  m_velocity.setX(dataToReceive[counter++]);
  m_velocity.setY(dataToReceive[counter++]);
  m_velocity.setZ(dataToReceive[counter++]);
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void MixUEqTotEnergy::computeSlopesMixture(const Mixture &sLeft, const Mixture &sRight, const double& distance)
{
  m_velocity.setX((sRight.getVelocity().getX() - sLeft.getVelocity().getX()) / distance);
  m_velocity.setY((sRight.getVelocity().getY() - sLeft.getVelocity().getY()) / distance);
  m_velocity.setZ((sRight.getVelocity().getZ() - sLeft.getVelocity().getZ()) / distance);
}

//***************************************************************************

void MixUEqTotEnergy::setToZero()
{
  m_velocity.setX(0.); m_velocity.setY(0.); m_velocity.setZ(0.);
}

//***************************************************************************

void MixUEqTotEnergy::extrapolate(const Mixture &slope, const double& distance)
{
  m_velocity.setX(m_velocity.getX() + slope.getVelocity().getX() * distance);
  m_velocity.setY(m_velocity.getY() + slope.getVelocity().getY() * distance);
  m_velocity.setZ(m_velocity.getZ() + slope.getVelocity().getZ() * distance);
}

//***************************************************************************

void MixUEqTotEnergy::limitSlopes(const Mixture &slopeGauche, const Mixture &slopeDroite, Limiter& globalLimiter)
{
  m_velocity.setX(globalLimiter.limiteSlope(slopeGauche.getVelocity().getX(), slopeDroite.getVelocity().getX()));
  m_velocity.setY(globalLimiter.limiteSlope(slopeGauche.getVelocity().getY(), slopeDroite.getVelocity().getY()));
  m_velocity.setZ(globalLimiter.limiteSlope(slopeGauche.getVelocity().getZ(), slopeDroite.getVelocity().getZ()));
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int MixUEqTotEnergy::numberOfTransmittedSlopes() const
{
	return 3;
}

//***************************************************************************

void MixUEqTotEnergy::fillBufferSlopes(double* buffer, int& counter) const
{
	buffer[++counter] = m_velocity.getX();
	buffer[++counter] = m_velocity.getY();
	buffer[++counter] = m_velocity.getZ();
}

//***************************************************************************

void MixUEqTotEnergy::getBufferSlopes(double* buffer, int& counter)
{
	m_velocity.setX(buffer[++counter]);
	m_velocity.setY(buffer[++counter]);
	m_velocity.setZ(buffer[++counter]);
}

//****************************************************************************
//******************************* ACCESSORS **********************************
//****************************************************************************

void MixUEqTotEnergy::setPressure(const double& p) { m_pressure = p; }

//***************************************************************************

void MixUEqTotEnergy::setVelocity(const double& u, const double& v, const double& w) { m_velocity.setXYZ(u, v, w); }

//***************************************************************************

void MixUEqTotEnergy::setVelocity(const Coord& vit) { m_velocity = vit; }

//***************************************************************************

void MixUEqTotEnergy::setU(const double& u) { m_velocity.setX(u); }

//***************************************************************************

void MixUEqTotEnergy::setV(const double& v) { m_velocity.setY(v); }

//***************************************************************************

void MixUEqTotEnergy::setW(const double& w) { m_velocity.setZ(w); }

//****************************************************************************
//***************************** OPERATORS ************************************
//****************************************************************************

void MixUEqTotEnergy::changeSign()
{
  m_density = -m_density;
  m_pressure = -m_pressure;
  m_velocity = m_velocity*-1.;
}

//***************************************************************************

void MixUEqTotEnergy::multiplyAndAdd(const Mixture &slopesMixtureTemp, const double& coeff)
{
  m_density += slopesMixtureTemp.getDensity()*coeff;
  m_pressure += slopesMixtureTemp.getPressure()*coeff;
  m_velocity += slopesMixtureTemp.getVelocity()*coeff;
}

//***************************************************************************

void MixUEqTotEnergy::divide(const double& coeff)
{
  m_density /= coeff;
  m_pressure /= coeff;
  m_velocity /= coeff;
}
