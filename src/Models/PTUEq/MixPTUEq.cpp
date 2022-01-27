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

#include "MixPTUEq.h"

using namespace tinyxml2;

//***************************************************************************

MixPTUEq::MixPTUEq() :m_density(0.), m_pressure(0.), m_velocity(0), m_energy(0.), m_totalEnergy(0.), m_PTUEqSoundSpeed(0.) {}

//***************************************************************************

MixPTUEq::MixPTUEq(XMLElement* state, std::string fileName) :
  m_density(0.), m_pressure(0.), m_energy(0.), m_totalEnergy(0.), m_PTUEqSoundSpeed(0.)
{
  XMLElement* sousElement(state->FirstChildElement("mixture"));
  if (sousElement == NULL) throw ErrorXMLElement("mixture", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  XMLElement* dataMix(sousElement->FirstChildElement("dataMix"));
  if (dataMix == NULL) throw ErrorXMLElement("dataMix", fileName, __FILE__, __LINE__);
  //pressure
  error = dataMix->QueryDoubleAttribute("pressure", &m_pressure);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("pressure", fileName, __FILE__, __LINE__);
  //temperature
  error = dataMix->QueryDoubleAttribute("temperature", &m_temperature);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("temperature", fileName, __FILE__, __LINE__);

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

MixPTUEq::~MixPTUEq(){}

//***************************************************************************

void MixPTUEq::allocateAndCopyMixture(Mixture** mixture)
{
  *mixture = new MixPTUEq(*this);
}

//***************************************************************************

void MixPTUEq::copyMixture(Mixture &mixture)
{
  m_density = mixture.getDensity();
  m_pressure = mixture.getPressure();
  m_temperature = mixture.getTemperature();
  m_velocity = mixture.getVelocity();
  m_energy = mixture.getEnergy();
  m_totalEnergy = mixture.getTotalEnergy();
  m_PTUEqSoundSpeed = mixture.getMixSoundSpeed();
}

//***************************************************************************

double MixPTUEq::computeDensity(const double* alphak, const double* rhok)
{
  double rho(0.);
  for(int k=0;k<numberPhases;k++)
  {
      rho += alphak[k]*rhok[k];
  }
  return rho;
}

//***************************************************************************

double MixPTUEq::computePressure(const double* alphak, const double* pk)
{
  double p(0.);
  for(int k=0;k<numberPhases;k++)
  {
      p += alphak[k]*pk[k];
  }
  return p;
}

//***************************************************************************

double MixPTUEq::computePressure(double* masses, const double& mixInternalEnerg, Phase** phases)
{
  //Restrictions //FP//TODO// to improve
  if (numberPhases > 2) Errors::errorMessage("more than two phases not permitted in thermal equilibrium model : MixPTUEq::computePressure");
  for (int k = 0; k < numberPhases; k++) {
    if (phases[k]->getEos()->getType() != TypeEOS::IG && phases[k]->getEos()->getType() != TypeEOS::SG) { Errors::errorMessage("Only IG and SG permitted in thermal equilibrium model: MixPTUEq::computePressure"); }
  }

  double rhoMel(0.);
  for (int k = 0; k < numberPhases; k++) {
    rhoMel += masses[k];
  }

  //Formulae of pressure for 2 phases goverened by SG EOS (Le Martelot, 2013, thesis)
  double gamma1 = phases[0]->getEos()->getGamma();
  double pInf1 = phases[0]->getEos()->getPInf();
  double cv1 = phases[0]->getEos()->getCv();
  double e01 = phases[0]->getEos()->getERef();
  double Y1 = masses[0] / rhoMel;

  double gamma2 = phases[1]->getEos()->getGamma();
  double pInf2 = phases[1]->getEos()->getPInf();
  double cv2 = phases[1]->getEos()->getCv();
  double e02 = phases[1]->getEos()->getERef();
  double Y2 = masses[1] / rhoMel;

  double q = Y1*e01 + Y2*e02;
  double cvMel = Y1*cv1 + Y2*cv2;
  
  double A1 = Y1*(gamma1 - 1.)*cv1 / cvMel*(rhoMel*(mixInternalEnerg - q) - pInf1);
  double A2 = Y2*(gamma2 - 1.)*cv2 / cvMel*(rhoMel*(mixInternalEnerg - q) - pInf2);

  m_pressure = 0.5*(A1 + A2 - (pInf1 + pInf2)) + sqrt(0.25*(A2 - A1 - (pInf2 - pInf1))*(A2 - A1 - (pInf2 - pInf1)) + A1*A2);
 
  return m_pressure;
}

//***************************************************************************

double MixPTUEq::computeTemperature(double* masses, const double& pressure, Phase** phases)
{
  //Restrictions //FP//TODO// to improve
  for (int k = 0; k < numberPhases; k++) {
    if (phases[k]->getEos()->getType() != TypeEOS::IG && phases[k]->getEos()->getType() != TypeEOS::SG) { Errors::errorMessage("Only IG and SG permitted in thermal equilibrium model: MixPTUEq::computePressure"); }
  }

  double rhoMel(0.);
  for (int k = 0; k < numberPhases; k++) {
    rhoMel += masses[k];
  }

  //Formulae for phases goverened by SG EOS (Le Martelot, 2013, phd thesis)
  double gammak, pInfk, cvk, Yk;
  m_temperature = 0.;
  for (int k = 0; k < numberPhases; k++) {
    gammak = phases[k]->getEos()->getGamma();
    pInfk = phases[k]->getEos()->getPInf();
    cvk = phases[k]->getEos()->getCv();
    Yk = masses[k] / rhoMel;
    m_temperature += Yk*(gammak - 1.)*cvk / (pressure + pInfk);
  }
  m_temperature = 1. / (m_temperature*rhoMel);

  return m_temperature;
}

//***************************************************************************

double MixPTUEq::computeInternalEnergy(const double* Yk, const double* ek)
{
  double e(0.);
  for(int k=0;k<numberPhases;k++)
  {
      e += Yk[k]*ek[k];
  }
  return e;
}

//***************************************************************************

double MixPTUEq::computeFrozenSoundSpeed(const double* Yk, const double* ck)
{
  double cF(0.);
  for(int k=0;k<numberPhases;k++)
  {
      cF += Yk[k]*ck[k]*ck[k];
  }
  return sqrt(cF);
}

//***************************************************************************

double MixPTUEq::computeTemperatureIsentrope(const double* Yk, const double& p0, const double& T0, const double& p, double* dTdp)
{
  //Restrictions //FP//TODO// to improve
  if (numberPhases > 2) Errors::errorMessage("more than two phases not permitted in thermal equilibrium model : MixPTUEq::computeTemperatureIsentrope");
  for (int k = 0; k < numberPhases; k++) {
    if (TB->eos[k]->getType() != TypeEOS::IG && TB->eos[k]->getType() != TypeEOS::SG) { Errors::errorMessage("Only IG and SG permitted in thermal equilibrium model: MixPTUEq::computeTemperatureIsentrope"); }
  }

  //Formulae for phases goverened by SG EOS
  double T(T0), cM(0.); 
  if (dTdp != NULL) *dTdp = 0.;
  for (int k = 0; k < numberPhases; k++) {
    cM += Yk[k] * TB->eos[k]->getGamma()* TB->eos[k]->getCv();
  }
  double puissance(0.), fk;
  for (int k = 0; k < numberPhases; k++) {
    puissance = (TB->eos[k]->getGamma() - 1.)*Yk[k] * TB->eos[k]->getCv()/cM;
    fk = std::pow((p + TB->eos[k]->getPInf()) / ((p0 + TB->eos[k]->getPInf())), puissance);
    T *= fk;
    if (dTdp != NULL) *dTdp += puissance/ (p + TB->eos[k]->getPInf());
  }
  if (dTdp != NULL) *dTdp *= T;

  return T;
}

//***************************************************************************

double MixPTUEq::computeEnthalpyIsentrope(const double* Yk, const double& p0, const double& T0, const double& p, double* dhdp)
{
  //Restrictions //FP//TODO// to improve
  for (int k = 0; k < numberPhases; k++) {
    if (TB->eos[k]->getType() != TypeEOS::IG && TB->eos[k]->getType() != TypeEOS::SG) { Errors::errorMessage("Only IG and SG permitted in thermal equilibrium model: MixPTUEq::computeEnthalpyIsentrope"); }
  }

  double dTdp(0.);
  double T = this->computeTemperatureIsentrope(Yk, p0, T0, p, &dTdp);
  //Formulae for phases goverened by SG EOS
  double h(0.);
  if (dhdp != NULL) *dhdp = 0.;
  for (int k = 0; k < numberPhases; k++) {
    h += Yk[k] * (TB->eos[k]->getGamma()*TB->eos[k]->getCv()*T + TB->eos[k]->getERef());
    if (dhdp != NULL) *dhdp += Yk[k] * TB->eos[k]->getGamma()*TB->eos[k]->getCv()*dTdp;
  }
  
  return h;
}

//***************************************************************************

double MixPTUEq::computeVolumeIsentrope(const double* Yk, const double& p0, const double& T0, const double& p, double* dvdp)
{
  //Restrictions //FP//TODO// to improve
  for (int k = 0; k < numberPhases; k++) {
    if (TB->eos[k]->getType() != TypeEOS::IG && TB->eos[k]->getType() != TypeEOS::SG) { Errors::errorMessage("Only IG and SG permitted in thermal equilibrium model: MixPTUEq::computeVolumeIsentrope"); }
  }

  double dTdp(0.);
  double T = this->computeTemperatureIsentrope(Yk, p0, T0, p, &dTdp);
  //Formulae for phases goverened by SG EOS
  double v(0.), vk(0.), dvk(0.); 
  if (dvdp != NULL) *dvdp = 0.;
  for (int k = 0; k < numberPhases; k++) {
    vk = ((TB->eos[k]->getGamma() - 1.)*TB->eos[k]->getCv()*T)/ (p + TB->eos[k]->getPInf());
    dvk = ((TB->eos[k]->getGamma() - 1.)*TB->eos[k]->getCv()*dTdp - vk) / (p + TB->eos[k]->getPInf());
    v += Yk[k] * vk;
    if (dvdp != NULL) *dvdp += Yk[k] * dvk;
  }

  return v;
}

//***************************************************************************

void MixPTUEq::computeMixtureVariables(Phase** vecPhase)
{
  //mixture density and pressure
  m_density = 0.;
  for (int k = 0; k < numberPhases; k++) {
    m_density += vecPhase[k]->getAlpha()*vecPhase[k]->getDensity();
  }
  //Mass fraction
  for (int k = 0; k < numberPhases; k++) {
    TB->Yk[k] = vecPhase[k]->getAlpha()*vecPhase[k]->getDensity() / m_density;
  }
  //Specific internal energy, speed of sound (frozen for now //FP//TODO//CHANGE sound speed)
  m_energy = 0.;
  m_PTUEqSoundSpeed = 0.;
  for (int k = 0; k < numberPhases; k++) {
    m_energy += TB->Yk[k] * vecPhase[k]->getEnergy();
    m_PTUEqSoundSpeed += TB->Yk[k] * vecPhase[k]->getSoundSpeed()*vecPhase[k]->getSoundSpeed();
  }
  m_PTUEqSoundSpeed = sqrt(m_PTUEqSoundSpeed);  
  //m_totalEnergy cannot be computed here because depending on extra additional energies
}

//***************************************************************************

void MixPTUEq::internalEnergyToTotalEnergy(std::vector<QuantitiesAddPhys*>& vecGPA)
{
  m_totalEnergy = m_energy + 0.5*m_velocity.squaredNorm();
  for (unsigned int pa = 0; pa < vecGPA.size(); pa++) {
    m_totalEnergy += vecGPA[pa]->computeEnergyAddPhys()/m_density; //Caution /m_density important
  }
}

//***************************************************************************

void MixPTUEq::totalEnergyToInternalEnergy(std::vector<QuantitiesAddPhys*>& vecGPA)
{
  m_energy = m_totalEnergy - 0.5*m_velocity.squaredNorm();
  for (unsigned int pa = 0; pa < vecGPA.size(); pa++) {
    m_energy -= vecGPA[pa]->computeEnergyAddPhys()/m_density; //Caution /m_density important
  }
}

//***************************************************************************

void MixPTUEq::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void MixPTUEq::reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
//**************************** DATA PRINTING *********************************
//****************************************************************************

double MixPTUEq::returnScalar(const int& numVar) const
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

Coord MixPTUEq::returnVector(const int& numVar) const
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

std::string MixPTUEq::returnNameScalar(const int& numVar) const
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

std::string MixPTUEq::returnNameVector(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Velocity_Mixture_norm"; break;
  default:
    return "NoName"; break;
  }
}

//****************************************************************************
//**************************** DATA READING **********************************
//****************************************************************************

void MixPTUEq::setScalar(const int& numVar, const double& value)
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
    Errors::errorMessage("numVar not found in MixPTUEq::setScalar"); break;
  }
}

//***************************************************************************

void MixPTUEq::setVector(const int& numVar, const Coord& value)
{
  switch (numVar)
  {
  case 1:
    m_velocity = value; break;
  default:
    Errors::errorMessage("numVar not found in MixPTUEq::setVector"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int MixPTUEq::numberOfTransmittedVariables() const
{
  //3 scalar + 1 vector : 6 variables
  return 6;
}

//***************************************************************************

void MixPTUEq::fillBuffer(double* buffer, int& counter) const
{
  buffer[++counter] = m_pressure;
  buffer[++counter] = m_temperature;
  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();
  buffer[++counter] = m_totalEnergy;
}

//***************************************************************************

void MixPTUEq::getBuffer(double* buffer, int& counter)
{
  m_pressure = buffer[++counter];
  m_temperature = buffer[++counter];
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);
  m_totalEnergy = buffer[++counter];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void MixPTUEq::computeSlopesMixture(const Mixture &sLeft, const Mixture &sRight, const double& distance)
{
  m_pressure = (sRight.getPressure() - sLeft.getPressure()) / distance;
  m_temperature = (sRight.getTemperature() - sLeft.getTemperature()) / distance;
  m_velocity.setX((sRight.getVelocity().getX() - sLeft.getVelocity().getX()) / distance);
  m_velocity.setY((sRight.getVelocity().getY() - sLeft.getVelocity().getY()) / distance);
  m_velocity.setZ((sRight.getVelocity().getZ() - sLeft.getVelocity().getZ()) / distance);
}

//***************************************************************************

void MixPTUEq::setToZero()
{
  m_pressure = 0.; m_temperature = 0.;
  m_velocity.setX(0.); m_velocity.setY(0.); m_velocity.setZ(0.);
}

//***************************************************************************

void MixPTUEq::extrapolate(const Mixture &slope, const double& distance)
{
  m_pressure += slope.getPressure()*distance;
  m_temperature += slope.getTemperature()*distance;
  m_velocity.setX(m_velocity.getX() + slope.getVelocity().getX() * distance);
  m_velocity.setY(m_velocity.getY() + slope.getVelocity().getY() * distance);
  m_velocity.setZ(m_velocity.getZ() + slope.getVelocity().getZ() * distance);
}

//***************************************************************************

void MixPTUEq::limitSlopes(const Mixture &slopeGauche, const Mixture &slopeDroite, Limiter& globalLimiter)
{
  m_pressure = globalLimiter.limiteSlope(slopeGauche.getPressure(), slopeDroite.getPressure());
  m_temperature = globalLimiter.limiteSlope(slopeGauche.getTemperature(), slopeDroite.getTemperature());
  m_velocity.setX(globalLimiter.limiteSlope(slopeGauche.getVelocity().getX(), slopeDroite.getVelocity().getX()));
  m_velocity.setY(globalLimiter.limiteSlope(slopeGauche.getVelocity().getY(), slopeDroite.getVelocity().getY()));
  m_velocity.setZ(globalLimiter.limiteSlope(slopeGauche.getVelocity().getZ(), slopeDroite.getVelocity().getZ()));
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int MixPTUEq::numberOfTransmittedSlopes() const
{
	return 5;
}

//***************************************************************************

void MixPTUEq::fillBufferSlopes(double* buffer, int& counter) const
{
  buffer[++counter] = m_pressure;
  buffer[++counter] = m_temperature;
	buffer[++counter] = m_velocity.getX();
	buffer[++counter] = m_velocity.getY();
	buffer[++counter] = m_velocity.getZ();
}

//***************************************************************************

void MixPTUEq::getBufferSlopes(double* buffer, int& counter)
{
  m_pressure = buffer[++counter];
  m_temperature = buffer[++counter];
	m_velocity.setX(buffer[++counter]);
	m_velocity.setY(buffer[++counter]);
	m_velocity.setZ(buffer[++counter]);
}

//****************************************************************************
//******************************* ACCESSORS **********************************
//****************************************************************************

void MixPTUEq::setPressure(const double& p) { m_pressure = p; }

//***************************************************************************

void MixPTUEq::setVelocity(const double& u, const double& v, const double& w) { m_velocity.setXYZ(u, v, w); }

//***************************************************************************

void MixPTUEq::setVelocity(const Coord& vit) { m_velocity = vit; }

//***************************************************************************

void MixPTUEq::setU(const double& u) { m_velocity.setX(u); }

//***************************************************************************

void MixPTUEq::setV(const double& v) { m_velocity.setY(v); }

//***************************************************************************

void MixPTUEq::setW(const double& w) { m_velocity.setZ(w); }

//***************************************************************************

void MixPTUEq::setTotalEnergy(double& totalEnergy)
{
  m_totalEnergy = totalEnergy;
}

//****************************************************************************
//***************************** OPERATORS ************************************
//****************************************************************************

void MixPTUEq::changeSign()
{
  m_pressure = -m_pressure;
  m_temperature = -m_temperature;
  m_velocity = m_velocity*-1.;
}

//***************************************************************************

void MixPTUEq::multiplyAndAdd(const Mixture &slopesMixtureTemp, const double& coeff)
{
  m_pressure += slopesMixtureTemp.getPressure()*coeff;
  m_temperature += slopesMixtureTemp.getTemperature()*coeff;
  m_velocity += slopesMixtureTemp.getVelocity()*coeff;
}

//***************************************************************************

void MixPTUEq::divide(const double& coeff)
{
  m_pressure /= coeff;
  m_temperature /= coeff;
  m_velocity /= coeff;
}
