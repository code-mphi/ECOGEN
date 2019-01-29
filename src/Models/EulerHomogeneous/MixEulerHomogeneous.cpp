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

//! \file      MixEulerHomogeneous.cpp
//! \author    K. Schmidmayer, F. Petitpas
//! \version   1.0
//! \date      December 19 2017

#include <cmath>
#include "MixEulerHomogeneous.h"

using namespace std;
using namespace tinyxml2;

//***************************************************************************

MixEulerHomogeneous::MixEulerHomogeneous() :m_density(0.), m_pressure(0.), m_velocity(0), m_energie(0.), m_totalEnergy(0.), m_EqSoundSpeed(0.) {}

//***************************************************************************

MixEulerHomogeneous::MixEulerHomogeneous(XMLElement *state, string fileName) :
  m_density(0.), m_pressure(0.), m_energie(0.), m_totalEnergy(0.), m_EqSoundSpeed(0.)
{
  XMLElement *sousElement(state->FirstChildElement("mixture"));
  if (sousElement == NULL) throw ErrorXMLElement("mixture", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;

  XMLElement *dataMix(sousElement->FirstChildElement("dataMix"));
  if (dataMix == NULL) throw ErrorXMLElement("dataMix", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  //pressure
  error = dataMix->QueryDoubleAttribute("pressure", &m_pressure);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("pressure", fileName, __FILE__, __LINE__);
  ////temperature
  //error = dataMix->QueryDoubleAttribute("temperature", &m_temperature);
  //if (error != XML_NO_ERROR) throw ErrorXMLAttribut("temperature", fileName, __FILE__, __LINE__);

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

MixEulerHomogeneous::~MixEulerHomogeneous() {}

//***************************************************************************

void MixEulerHomogeneous::allocateAndCopyMixture(Mixture **mixture)
{
  *mixture = new MixEulerHomogeneous(*this);
}

//***************************************************************************

void MixEulerHomogeneous::copyMixture(Mixture &mixture)
{
  m_density = mixture.getDensity();
  m_pressure = mixture.getPressure();
  m_velocity = mixture.getVelocity();
  m_energie = mixture.getEnergy();
  m_totalEnergy = mixture.getTotalEnergy();
  m_EqSoundSpeed = mixture.getMixSoundSpeed();
}

//***************************************************************************

double MixEulerHomogeneous::computeDensity(const double *alphak, const double *rhok, const int &numberPhases)
{
  double rho(0.);
  for (int k = 0; k<numberPhases; k++)
  {
    rho += alphak[k] * rhok[k];
  }
  return rho;
}

//***************************************************************************

double MixEulerHomogeneous::computePressure(const double *alphak, const double *pk, const int &numberPhases)
{
  double p(0.);
  for (int k = 0; k<numberPhases; k++)
  {
    p += alphak[k] * pk[k];
  }
  return p;
}

//***************************************************************************

double MixEulerHomogeneous::computePressure(double masse, const double &internalEnergy, Phase **phases, Mixture *mixture, const int &numberPhases, const int &liq, const int &vap)
{
  double pressure(0.);

  //Restrictions
  if (numberPhases > 2) Errors::errorMessage("more than two phases not permitted in thermodynamical equilibrium model : MixEulerHomogeneous::computePressure");
  if (phases[vap]->getEos()->getType() != "IG" && phases[liq]->getEos()->getType() != "SG") { Errors::errorMessage("Only IG for vapor and SG for liquid permitted in thermodyanmical equilibrium model : MixEulerHomogeneous::computePressure"); }

  //Iterative process for pressure determination based on energy conservation (e=Sum(Yk*ek))
  int iteration(0);
  pressure = phases[liq]->getPressure(); //pressure estimate
  double f(0.), df(1.);
  double alphaVap, dalphaVap, rhoLiq, drhoLiq, rhoVap, drhoVap, rhoeLiq, drhoeLiq, rhoeVap, drhoeVap, Tsat, dTsat;
  do {
    pressure -= f / df; iteration++;
    if (iteration > 50) {
      errors.push_back(Errors("not converged in MixEulerHomogeneous::computePressure", __FILE__, __LINE__));
      break;
    }
    Tsat = mixture->computeTsat(phases[liq]->getEos(), phases[vap]->getEos(), pressure, &dTsat);
    rhoVap = phases[vap]->getEos()->computeDensitySaturation(pressure, Tsat, dTsat, &drhoVap);
    rhoLiq = phases[liq]->getEos()->computeDensitySaturation(pressure, Tsat, dTsat, &drhoLiq);
    alphaVap = (masse - rhoLiq) / (rhoVap - rhoLiq);
    dalphaVap = (-drhoLiq*(rhoVap - rhoLiq) - (masse - rhoLiq)*(drhoVap - drhoLiq)) / ((rhoVap - rhoLiq)*(rhoVap - rhoLiq));
    rhoeVap = phases[vap]->getEos()->computeDensityEnergySaturation(pressure, rhoVap, drhoVap, &drhoeVap);
    rhoeLiq = phases[liq]->getEos()->computeDensityEnergySaturation(pressure, rhoLiq, drhoLiq, &drhoeLiq);
    f = masse*internalEnergy - alphaVap*(rhoeVap - rhoeLiq) - rhoeLiq;
    df = -dalphaVap*(rhoeVap - rhoeLiq) - alphaVap*(drhoeVap - drhoeLiq) - drhoeLiq;
  } while (abs(f / (masse*internalEnergy)) > 1e-10);

  return pressure;
}

//***************************************************************************

double MixEulerHomogeneous::computeInternalEnergy(const double *Yk, const double *ek, const int &numberPhases)
{
  double e(0.);
  for (int k = 0; k<numberPhases; k++)
  {
    e += Yk[k] * ek[k];
  }
  return e;
}

//***************************************************************************

double MixEulerHomogeneous::computeFrozenSoundSpeed(const double *Yk, const double *ck, const int &numberPhases)
{
  double cF(0.);
  for (int k = 0; k<numberPhases; k++)
  {
    cF += Yk[k] * ck[k] * ck[k];
  }
  return sqrt(cF);
}

//***************************************************************************

void MixEulerHomogeneous::computeMixtureVariables(Phase **vecPhase, const int &numberPhases)
{
  //mixture density and pressure
  m_density = 0.;
  for (int k = 0; k < numberPhases; k++) {
    m_density += vecPhase[k]->getAlpha()*vecPhase[k]->getDensity();
  }
  //Mass fraction
  for (int k = 0; k < numberPhases; k++) {
    vecPhase[k]->computeMassFraction(m_density);
  }
  //Specific internal energy, speed of sounds and total specific energy
  m_energie = 0.;
  m_EqSoundSpeed = 0.;
  for (int k = 0; k < numberPhases; k++) {
    m_energie += vecPhase[k]->getY() *vecPhase[k]->getEnergy();
    m_EqSoundSpeed += vecPhase[k]->getY() * vecPhase[k]->getSoundSpeed()*vecPhase[k]->getSoundSpeed();
  }
  m_EqSoundSpeed = sqrt(m_EqSoundSpeed);
  //m_totalEnergy cannot be computed here because depending on extra additional energies
}

//***************************************************************************

void MixEulerHomogeneous::internalEnergyToTotalEnergy(vector<QuantitiesAddPhys*> &vecGPA)
{
  m_totalEnergy = m_energie + 0.5*m_velocity.squaredNorm();
  for (unsigned int pa = 0; pa < vecGPA.size(); pa++) {
    m_totalEnergy += vecGPA[pa]->computeEnergyAddPhys() / m_density; //Caution /m_density important
  }
}

//***************************************************************************

void MixEulerHomogeneous::totalEnergyToInternalEnergy(vector<QuantitiesAddPhys*> &vecGPA)
{
  m_energie = m_totalEnergy - 0.5*m_velocity.squaredNorm();
  for (unsigned int pa = 0; pa < vecGPA.size(); pa++) {
    m_energie -= vecGPA[pa]->computeEnergyAddPhys() / m_density; //Caution /m_density important
  }
}

//***************************************************************************

void MixEulerHomogeneous::localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void MixEulerHomogeneous::reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
//**************************** DATA PRINTING *********************************
//****************************************************************************

double MixEulerHomogeneous::returnScalar(const int &numVar) const
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

Coord MixEulerHomogeneous::returnVector(const int &numVar) const
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

string MixEulerHomogeneous::returnNameScalar(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Density_Mixture"; break;
  case 2:
    return "Pressure_Mixture"; break;
  case 3:
    return "Temperature_Mixture"; break;
  default:
    return "NoName"; break;
  }
}

//***************************************************************************

string MixEulerHomogeneous::returnNameVector(const int &numVar) const
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

void MixEulerHomogeneous::setScalar(const int &numVar, const double &value)
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
    Errors::errorMessage("numVar not found in MixEulerHomogeneous::setScalar"); break;
  }
}

//***************************************************************************

void MixEulerHomogeneous::setVector(const int &numVar, const Coord &value)
{
  switch (numVar)
  {
  case 1:
    m_velocity = value; break;
  default:
    Errors::errorMessage("numVar not found in MixEulerHomogeneous::setVector"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int MixEulerHomogeneous::numberOfTransmittedVariables() const
{
  //2 scalar + 1 vector : 5 variables
  return 5;
}

//***************************************************************************

void MixEulerHomogeneous::fillBuffer(double *buffer, int &counter) const
{
  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();
  buffer[++counter] = m_pressure;
  buffer[++counter] = m_totalEnergy;
}

//***************************************************************************

void MixEulerHomogeneous::getBuffer(double *buffer, int &counter)
{
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);
  m_pressure = buffer[++counter];
  m_totalEnergy = buffer[++counter];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void MixEulerHomogeneous::computeSlopesMixture(const Mixture &sLeft, const Mixture &sRight, const double &distance)
{
  m_pressure = (sRight.getPressure() - sLeft.getPressure()) / distance;
  m_velocity.setX((sRight.getVelocity().getX() - sLeft.getVelocity().getX()) / distance);
  m_velocity.setY((sRight.getVelocity().getY() - sLeft.getVelocity().getY()) / distance);
  m_velocity.setZ((sRight.getVelocity().getZ() - sLeft.getVelocity().getZ()) / distance);
}

//***************************************************************************

void MixEulerHomogeneous::setToZero()
{
  m_pressure = 0.;
  m_velocity.setX(0.); m_velocity.setY(0.); m_velocity.setZ(0.);
}

//***************************************************************************

void MixEulerHomogeneous::extrapolate(const Mixture &slope, const double &distance)
{
  m_pressure += slope.getPressure()*distance;
  m_velocity.setX(m_velocity.getX() + slope.getVelocity().getX() * distance);
  m_velocity.setY(m_velocity.getY() + slope.getVelocity().getY() * distance);
  m_velocity.setZ(m_velocity.getZ() + slope.getVelocity().getZ() * distance);
}

//***************************************************************************

void MixEulerHomogeneous::limitSlopes(const Mixture &slopeGauche, const Mixture &slopeDroite, Limiter &globalLimiter)
{
  m_pressure = globalLimiter.limiteSlope(slopeGauche.getPressure(), slopeDroite.getPressure());
  m_velocity.setX(globalLimiter.limiteSlope(slopeGauche.getVelocity().getX(), slopeDroite.getVelocity().getX()));
  m_velocity.setY(globalLimiter.limiteSlope(slopeGauche.getVelocity().getY(), slopeDroite.getVelocity().getY()));
  m_velocity.setZ(globalLimiter.limiteSlope(slopeGauche.getVelocity().getZ(), slopeDroite.getVelocity().getZ()));
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int MixEulerHomogeneous::numberOfTransmittedSlopes() const
{
	return 4;
}

//***************************************************************************

void MixEulerHomogeneous::fillBufferSlopes(double *buffer, int &counter) const
{
  buffer[++counter] = m_pressure;
	buffer[++counter] = m_velocity.getX();
	buffer[++counter] = m_velocity.getY();
	buffer[++counter] = m_velocity.getZ();
}

//***************************************************************************

void MixEulerHomogeneous::getBufferSlopes(double *buffer, int &counter)
{
  m_pressure = buffer[++counter];
	m_velocity.setX(buffer[++counter]);
	m_velocity.setY(buffer[++counter]);
	m_velocity.setZ(buffer[++counter]);
}

//****************************************************************************
//****************************** ACCESSORS  **********************************
//****************************************************************************

double MixEulerHomogeneous::getDensity() const
{
  return m_density;
}

//***************************************************************************

double MixEulerHomogeneous::getPressure() const
{
  return m_pressure;
}

//***************************************************************************

double MixEulerHomogeneous::getU() const { return m_velocity.getX(); }
double MixEulerHomogeneous::getV() const { return m_velocity.getY(); }
double MixEulerHomogeneous::getW() const { return m_velocity.getZ(); }

//***************************************************************************

Coord MixEulerHomogeneous::getVelocity() const
{
  return m_velocity;
}

//***************************************************************************

double MixEulerHomogeneous::getEnergy() const
{
  return m_energie;
}

//***************************************************************************

double MixEulerHomogeneous::getTotalEnergy() const
{
  return m_totalEnergy;
}

//***************************************************************************

double MixEulerHomogeneous::getMixSoundSpeed() const
{
  return m_EqSoundSpeed;
}

//***************************************************************************

void MixEulerHomogeneous::setPressure(const double &p) { m_pressure = p; }

//***************************************************************************

void MixEulerHomogeneous::setTemperature(const double &T) { m_temperature = T; }

//***************************************************************************

void MixEulerHomogeneous::setVelocity(const double &u, const double &v, const double &w) { m_velocity.setXYZ(u, v, w); }

//***************************************************************************

void MixEulerHomogeneous::setVelocity(const Coord &vit) { m_velocity = vit; }

//***************************************************************************

void MixEulerHomogeneous::setU(const double &u) { m_velocity.setX(u); }

//***************************************************************************

void MixEulerHomogeneous::setV(const double &v) { m_velocity.setY(v); }

//***************************************************************************

void MixEulerHomogeneous::setW(const double &w) { m_velocity.setZ(w); }

//***************************************************************************

void MixEulerHomogeneous::setTotalEnergy(double &totalEnergy)
{
  m_totalEnergy = totalEnergy;
}

//****************************************************************************
//****************************** OPERATORS ***********************************
//****************************************************************************

void MixEulerHomogeneous::changeSign()
{
  m_pressure = -m_pressure;
  m_velocity = m_velocity*-1.;
}

//***************************************************************************

void MixEulerHomogeneous::multiplyAndAdd(const Mixture &slopesMixtureTemp, const double &coeff)
{
  m_pressure += slopesMixtureTemp.getPressure()*coeff;
  m_velocity += slopesMixtureTemp.getVelocity()*coeff;
}

//***************************************************************************

void MixEulerHomogeneous::divide(const double &coeff)
{
  m_pressure /= coeff;
  m_velocity /= coeff;
}
