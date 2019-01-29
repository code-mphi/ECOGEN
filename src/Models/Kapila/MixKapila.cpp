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

//! \file      MixKapila.cpp
//! \author    K. Schmidmayer, F. Petitpas
//! \version   1.0
//! \date      December 19 2017

#include <cmath>
#include "MixKapila.h"

using namespace std;
using namespace tinyxml2;

//***************************************************************************

MixKapila::MixKapila() :m_density(0.), m_pressure(0.), m_velocity(0), m_energie(0.), m_totalEnergy(0.), m_frozenSoundSpeed(0.), m_woodSoundSpeed(0.) {}

//***************************************************************************

MixKapila::MixKapila(XMLElement *state, string fileName) :
  m_density(0.), m_pressure(0.), m_energie(0.), m_totalEnergy(0.), m_frozenSoundSpeed(0.), m_woodSoundSpeed(0.)
{
  XMLElement *sousElement(state->FirstChildElement("mixture"));
  if (sousElement == NULL) throw ErrorXMLElement("mixture", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  XMLElement *dataMix(sousElement->FirstChildElement("dataMix"));
  if (dataMix == NULL) throw ErrorXMLElement("dataMix", fileName, __FILE__, __LINE__);
  //pressure
  error = dataMix->QueryDoubleAttribute("pressure", &m_pressure);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("pressure", fileName, __FILE__, __LINE__);
  //velocity
  XMLElement *velocity(sousElement->FirstChildElement("velocity"));
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

MixKapila::~MixKapila(){}

//***************************************************************************

void MixKapila::allocateAndCopyMixture(Mixture **mixture)
{
  *mixture = new MixKapila(*this);
}

//***************************************************************************

void MixKapila::copyMixture(Mixture &mixture)
{
  m_density = mixture.getDensity();
  m_pressure = mixture.getPressure();
  m_velocity = mixture.getVelocity();
  m_energie = mixture.getEnergy();
  m_totalEnergy = mixture.getTotalEnergy();
  m_frozenSoundSpeed = mixture.getFrozenSoundSpeed();
  m_woodSoundSpeed = mixture.getWoodSoundSpeed();
}

//***************************************************************************

double MixKapila::computeDensity(const double *alphak, const double *rhok, const int &numberPhases)
{
  double rho(0.);
  for(int k=0;k<numberPhases;k++)
  {
      rho += alphak[k]*rhok[k];
  }
  return rho;
}

//***************************************************************************

double MixKapila::computePressure(const double *alphak, const double *pk, const int &numberPhases)
{
  double p(0.);
  for(int k=0;k<numberPhases;k++)
  {
      p += alphak[k]*pk[k];
  }
  return p;
}

//***************************************************************************

double MixKapila::computeInternalEnergy(const double *Yk, const double *ek, const int &numberPhases)
{
  double e(0.);
  for(int k=0;k<numberPhases;k++)
  {
      e += Yk[k]*ek[k];
  }
  return e;
}

//***************************************************************************

double MixKapila::computeFrozenSoundSpeed(const double *Yk, const double *ck, const int &numberPhases)
{
  double cF(0.);
  for(int k=0;k<numberPhases;k++)
  {
      cF += Yk[k]*ck[k]*ck[k];
  }
  return sqrt(cF);
}

//***************************************************************************

void MixKapila::computeMixtureVariables(Phase **vecPhase, const int &numberPhases)
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
  //Specific internal energy and speed of sounds
  m_energie = 0.;
  m_frozenSoundSpeed = 0.;
  m_woodSoundSpeed = 0.;
  for (int k = 0; k < numberPhases; k++) {
    m_energie += vecPhase[k]->getY() * vecPhase[k]->getEnergy();
    m_frozenSoundSpeed += vecPhase[k]->getY() * vecPhase[k]->getSoundSpeed()*vecPhase[k]->getSoundSpeed();
    m_woodSoundSpeed += vecPhase[k]->getAlpha() / max((vecPhase[k]->getDensity()*vecPhase[k]->getSoundSpeed()*vecPhase[k]->getSoundSpeed()), epsilon);
  }
  m_frozenSoundSpeed = sqrt(m_frozenSoundSpeed);
  m_woodSoundSpeed = 1. / sqrt(m_density*m_woodSoundSpeed);
  //m_totalEnergy cannot be computed here because depending on extra additional energies
}

//***************************************************************************

void MixKapila::internalEnergyToTotalEnergy(vector<QuantitiesAddPhys*> &vecGPA)
{
  m_totalEnergy = m_energie + 0.5*m_velocity.squaredNorm();
  for (unsigned int pa = 0; pa < vecGPA.size(); pa++) {
    m_totalEnergy += vecGPA[pa]->computeEnergyAddPhys()/m_density; //Caution /m_density important
  }
}

//***************************************************************************

void MixKapila::totalEnergyToInternalEnergy(vector<QuantitiesAddPhys*> &vecGPA)
{
  m_energie = m_totalEnergy - 0.5*m_velocity.squaredNorm();
  for (unsigned int pa = 0; pa < vecGPA.size(); pa++) {
    m_energie -= vecGPA[pa]->computeEnergyAddPhys()/m_density; //Caution /m_density important
  }
}

//***************************************************************************

void MixKapila::localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void MixKapila::reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
//**************************** DATA PRINTING *********************************
//****************************************************************************

double MixKapila::returnScalar(const int &numVar) const
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

Coord MixKapila::returnVector(const int &numVar) const
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

string MixKapila::returnNameScalar(const int &numVar) const
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

string MixKapila::returnNameVector(const int &numVar) const
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

void MixKapila::setScalar(const int &numVar, const double &value)
{
  switch (numVar)
  {
  case 1:
    m_density = value; break;
  case 2:
    m_pressure = value; break;
  default:
    Errors::errorMessage("numVar not found in MixKapila::setScalar"); break;
  }
}

//***************************************************************************

void MixKapila::setVector(const int &numVar, const Coord &value)
{
  switch (numVar)
  {
  case 1:
    m_velocity = value; break;
  default:
    Errors::errorMessage("numVar not found in MixKapila::setVector"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int MixKapila::numberOfTransmittedVariables() const
{
  //1 scalar + 1 vector : 4 variables
  return 4;
}

//***************************************************************************

void MixKapila::fillBuffer(double *buffer, int &counter) const
{
  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();
  buffer[++counter] = m_totalEnergy;
}

//***************************************************************************

void MixKapila::getBuffer(double *buffer, int &counter)
{
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);
  m_totalEnergy = buffer[++counter];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void MixKapila::computeSlopesMixture(const Mixture &sLeft, const Mixture &sRight, const double &distance)
{
  m_velocity.setX((sRight.getVelocity().getX() - sLeft.getVelocity().getX()) / distance);
  m_velocity.setY((sRight.getVelocity().getY() - sLeft.getVelocity().getY()) / distance);
  m_velocity.setZ((sRight.getVelocity().getZ() - sLeft.getVelocity().getZ()) / distance);
}

//***************************************************************************

void MixKapila::setToZero()
{
  m_velocity.setX(0.); m_velocity.setY(0.); m_velocity.setZ(0.);
}

//***************************************************************************

void MixKapila::extrapolate(const Mixture &slope, const double &distance)
{
  m_velocity.setX(m_velocity.getX() + slope.getVelocity().getX() * distance);
  m_velocity.setY(m_velocity.getY() + slope.getVelocity().getY() * distance);
  m_velocity.setZ(m_velocity.getZ() + slope.getVelocity().getZ() * distance);
}

//***************************************************************************

void MixKapila::limitSlopes(const Mixture &slopeGauche, const Mixture &slopeDroite, Limiter &globalLimiter)
{
  m_velocity.setX(globalLimiter.limiteSlope(slopeGauche.getVelocity().getX(), slopeDroite.getVelocity().getX()));
  m_velocity.setY(globalLimiter.limiteSlope(slopeGauche.getVelocity().getY(), slopeDroite.getVelocity().getY()));
  m_velocity.setZ(globalLimiter.limiteSlope(slopeGauche.getVelocity().getZ(), slopeDroite.getVelocity().getZ()));
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int MixKapila::numberOfTransmittedSlopes() const
{
	return 3;
}

//***************************************************************************

void MixKapila::fillBufferSlopes(double *buffer, int &counter) const
{
	buffer[++counter] = m_velocity.getX();
	buffer[++counter] = m_velocity.getY();
	buffer[++counter] = m_velocity.getZ();
}

//***************************************************************************

void MixKapila::getBufferSlopes(double *buffer, int &counter)
{
	m_velocity.setX(buffer[++counter]);
	m_velocity.setY(buffer[++counter]);
	m_velocity.setZ(buffer[++counter]);
}

//****************************************************************************
//******************************* ACCESSORS **********************************
//****************************************************************************

double MixKapila::getDensity() const
{
  return m_density;
}

//***************************************************************************

double MixKapila::getPressure() const
{
  return m_pressure;
}

//***************************************************************************

double MixKapila::getU() const { return m_velocity.getX(); }
double MixKapila::getV() const { return m_velocity.getY(); }
double MixKapila::getW() const { return m_velocity.getZ(); }

//***************************************************************************

Coord MixKapila::getVelocity() const
{
  return m_velocity;
}

//***************************************************************************

double MixKapila::getEnergy() const
{
  return m_energie;
}

//***************************************************************************

double MixKapila::getTotalEnergy() const
{
  return m_totalEnergy;
}

//***************************************************************************

double MixKapila::getFrozenSoundSpeed() const
{
  return m_frozenSoundSpeed;
}

//***************************************************************************

double MixKapila::getWoodSoundSpeed() const
{
  return m_woodSoundSpeed;
}

//***************************************************************************

void MixKapila::setPressure(const double &p) { m_pressure = p; }

//***************************************************************************

void MixKapila::setVelocity(const double &u, const double &v, const double &w) { m_velocity.setXYZ(u, v, w); }

//***************************************************************************

void MixKapila::setVelocity(const Coord &vit) { m_velocity = vit; }

//***************************************************************************

void MixKapila::setU(const double &u) { m_velocity.setX(u); }

//***************************************************************************

void MixKapila::setV(const double &v) { m_velocity.setY(v); }

//***************************************************************************

void MixKapila::setW(const double &w) { m_velocity.setZ(w); }

//***************************************************************************

void MixKapila::setTotalEnergy(double &totalEnergy)
{
  m_totalEnergy = totalEnergy;
}

//****************************************************************************
//***************************** OPERATORS ************************************
//****************************************************************************

void MixKapila::changeSign()
{
  m_density = -m_density;
  m_pressure = -m_pressure;
  m_velocity = m_velocity*-1.;
}

//***************************************************************************

void MixKapila::multiplyAndAdd(const Mixture &slopesMixtureTemp, const double &coeff)
{
  m_density += slopesMixtureTemp.getDensity()*coeff;
  m_pressure += slopesMixtureTemp.getPressure()*coeff;
  m_velocity += slopesMixtureTemp.getVelocity()*coeff;
}

//***************************************************************************

void MixKapila::divide(const double &coeff)
{
  m_density /= coeff;
  m_pressure /= coeff;
  m_velocity /= coeff;
}
