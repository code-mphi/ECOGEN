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

#include "PhaseEulerKorteweg.h"
#include "../../Eos/Eos.h"

using namespace tinyxml2;

//***************************************************************************

PhaseEulerKorteweg::PhaseEulerKorteweg() : m_density(0.), m_omega(0.), m_eta(0.), m_pressure(0.), m_eos(0)
{
  m_velocity.setXYZ(0., 0., 0.);
  m_vectorP.setXYZ(0., 0., 0.);
}

//***************************************************************************

PhaseEulerKorteweg::PhaseEulerKorteweg(XMLElement* material, Eos* eos, std::string fileName) : m_density(0.), m_omega(0.), m_eta(0.), m_pressure(0.), m_velocity(0), m_vectorP(0), m_eos(eos)
{
  XMLElement* sousElement(material->FirstChildElement("dataFluid"));
  if (sousElement == NULL) throw ErrorXMLElement("dataFluid", fileName, __FILE__, __LINE__);

  //Attributes reading
  //------------------
  XMLError error;
  //Density
  error = sousElement->QueryDoubleAttribute("density", &m_density);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("density", fileName, __FILE__, __LINE__);
  //Velocity
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

PhaseEulerKorteweg::~PhaseEulerKorteweg(){}

//***************************************************************************

void PhaseEulerKorteweg::allocateAndCopyPhase(Phase** vecPhase)
{
  *vecPhase = new PhaseEulerKorteweg(*this);
}

//***************************************************************************

void PhaseEulerKorteweg::copyPhase(Phase &phase)
{
  m_density = phase.getDensity();
  m_omega = phase.getOmega();
  m_eta = phase.getEta();
  m_velocity = phase.getVelocity();
  m_vectorP = phase.getVectorP();
  m_eos = phase.getEos();
}

//***************************************************************************

void PhaseEulerKorteweg::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
  m_vectorP.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void PhaseEulerKorteweg::reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
  m_vectorP.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhaseEulerKorteweg::returnScalar(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_density; break;
  case 2:
    return m_omega; break;
  case 3:
    return m_eta; break;
  case 4:
    return m_pressure; break;
  default:
    return 0.; break;
  }
}

//***************************************************************************

Coord PhaseEulerKorteweg::returnVector(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_velocity; break;
  case 2:
    return m_vectorP; break;
  default:
    return 0; break;
  }
}

//***************************************************************************

std::string PhaseEulerKorteweg::returnNameScalar(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Density"; break;
  case 2:
    return "Omega"; break;
  case 3:
    return "Eta"; break;
  case 4:
    return "Pressure"; break;
  default:
    return "NoName"; break;
  }
}

//***************************************************************************

std::string PhaseEulerKorteweg::returnNameVector(const int& numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Velocity"; break;
  case 2:
    return "VectorP"; break;
  default:
    return "NoName"; break;
  }
}

//****************************************************************************
//************************* READING FROM FILE ********************************
//****************************************************************************

void PhaseEulerKorteweg::setScalar(const int& numVar, const double& value)
{
  switch (numVar)
  {
  case 1:
    m_density = value; break;
  case 2:
    m_omega = value; break;
  case 3:
    m_eta = value; break;
  case 4:
    m_pressure = value; break;
  default:
    Errors::errorMessage("numVar not found in PhaseEulerKorteweg::setScalar"); break;
  }
}

//****************************************************************************

void PhaseEulerKorteweg::setVector(const int& numVar, const Coord& value)
{
  switch (numVar)
  {
  case 1:
    m_velocity = value; break;
  case 2:
    m_vectorP = value; break;
  default:
    Errors::errorMessage("numVar not found in PhaseEulerKorteweg::setVector"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseEulerKorteweg::numberOfTransmittedVariables() const
{
  //10 variables (3 scalar + 2*3 vector + number EOS)
  return 10;
}

//***************************************************************************

void PhaseEulerKorteweg::fillBuffer(double* buffer, int& counter) const
{
  buffer[++counter] = m_density;
  buffer[++counter] = m_omega;
  buffer[++counter] = m_eta;

  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();

  buffer[++counter] = m_vectorP.getX();
  buffer[++counter] = m_vectorP.getY();
  buffer[++counter] = m_vectorP.getZ();

  buffer[++counter] = static_cast<double>(m_eos->getNumber());
}

//***************************************************************************

void PhaseEulerKorteweg::fillBuffer(std::vector<double>& dataToSend) const
{
  dataToSend.push_back(m_density);
  dataToSend.push_back(m_omega);
  dataToSend.push_back(m_eta);
  
  dataToSend.push_back(m_velocity.getX());
  dataToSend.push_back(m_velocity.getY());
  dataToSend.push_back(m_velocity.getZ());
    
  dataToSend.push_back(m_vectorP.getX());
  dataToSend.push_back(m_vectorP.getY());
  dataToSend.push_back(m_vectorP.getZ());

  dataToSend.push_back(static_cast<double>(m_eos->getNumber()));
}

//***************************************************************************

void PhaseEulerKorteweg::getBuffer(double* buffer, int& counter, Eos** eos)
{
  m_density = buffer[++counter];
  m_omega = buffer[++counter];
  m_eta = buffer[++counter];
  
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);

  m_vectorP.setX(buffer[++counter]);
  m_vectorP.setY(buffer[++counter]);
  m_vectorP.setZ(buffer[++counter]);

  m_eos = eos[static_cast<int>(buffer[++counter])];
}

//***************************************************************************

void PhaseEulerKorteweg::getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** eos)
{
  m_density = dataToReceive[counter++];
  m_omega = dataToReceive[counter++];
  m_eta = dataToReceive[counter++];

  m_velocity.setX(dataToReceive[counter++]);
  m_velocity.setY(dataToReceive[counter++]);
  m_velocity.setZ(dataToReceive[counter++]);

  m_vectorP.setX(dataToReceive[counter++]);
  m_vectorP.setY(dataToReceive[counter++]);
  m_vectorP.setZ(dataToReceive[counter++]);

  m_eos = eos[static_cast<int>(dataToReceive[counter++])];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void PhaseEulerKorteweg::computeSlopesPhase(const Phase &sLeft, const Phase &sRight, const double& distance)
{
  m_density = (sRight.getDensity() - sLeft.getDensity()) / distance;
  m_omega = (sRight.getOmega() - sLeft.getOmega()) / distance;
  m_eta =   (sRight.getEta() -   sLeft.getEta()) / distance;

  m_velocity.setX((sRight.getVelocity().getX() - sLeft.getVelocity().getX()) / distance);
  m_velocity.setY((sRight.getVelocity().getY() - sLeft.getVelocity().getY()) / distance);
  m_velocity.setZ((sRight.getVelocity().getZ() - sLeft.getVelocity().getZ()) / distance);

  m_vectorP.setX((sRight.getVectorP().getX() - sLeft.getVectorP().getX()) / distance);
  m_vectorP.setY((sRight.getVectorP().getY() - sLeft.getVectorP().getY()) / distance);
  m_vectorP.setZ((sRight.getVectorP().getZ() - sLeft.getVectorP().getZ()) / distance);
}

//***************************************************************************

void PhaseEulerKorteweg::setToZero()
{
  m_density = 0.;  m_omega = 0.; m_eta = 0.;
  m_velocity.setX(0.); m_velocity.setY(0.); m_velocity.setZ(0.);
  m_vectorP.setX(0.); m_vectorP.setY(0.); m_vectorP.setZ(0.);
}

//***************************************************************************

void PhaseEulerKorteweg::extrapolate(const Phase &slope, const double& distance)
{
  m_density += slope.getDensity() * distance;
  m_omega += slope.getOmega() * distance;
  m_eta += slope.getEta() * distance;
  
  m_velocity.setX(m_velocity.getX() + slope.getVelocity().getX() * distance);
  m_velocity.setY(m_velocity.getY() + slope.getVelocity().getY() * distance);
  m_velocity.setZ(m_velocity.getZ() + slope.getVelocity().getZ() * distance);

  m_vectorP.setX(m_vectorP.getX() + slope.getVectorP().getX() * distance);
  m_vectorP.setY(m_vectorP.getY() + slope.getVectorP().getY() * distance);
  m_vectorP.setZ(m_vectorP.getZ() + slope.getVectorP().getZ() * distance);
}

//***************************************************************************

void PhaseEulerKorteweg::limitSlopes(const Phase& slopeGauche, const Phase& slopeDroite, Limiter& globalLimiter, Limiter& /*volumeFractionLimiter*/)
{
  m_density = globalLimiter.limiteSlope(slopeGauche.getDensity(), slopeDroite.getDensity());
  m_omega = globalLimiter.limiteSlope(slopeGauche.getOmega(), slopeDroite.getOmega());
  m_eta = globalLimiter.limiteSlope(slopeGauche.getEta(), slopeDroite.getEta());
  
  m_velocity.setX(globalLimiter.limiteSlope(slopeGauche.getVelocity().getX(), slopeDroite.getVelocity().getX()));
  m_velocity.setY(globalLimiter.limiteSlope(slopeGauche.getVelocity().getY(), slopeDroite.getVelocity().getY()));
  m_velocity.setZ(globalLimiter.limiteSlope(slopeGauche.getVelocity().getZ(), slopeDroite.getVelocity().getZ()));

  m_vectorP.setX(globalLimiter.limiteSlope(slopeGauche.getVectorP().getX(), slopeDroite.getVectorP().getX()));
  m_vectorP.setY(globalLimiter.limiteSlope(slopeGauche.getVectorP().getY(), slopeDroite.getVectorP().getY()));
  m_vectorP.setZ(globalLimiter.limiteSlope(slopeGauche.getVectorP().getZ(), slopeDroite.getVectorP().getZ()));
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int PhaseEulerKorteweg::numberOfTransmittedSlopes() const
{
  return 9;
}

//***************************************************************************

void PhaseEulerKorteweg::fillBufferSlopes(double* buffer, int& counter) const
{
	buffer[++counter] = m_density;
	buffer[++counter] = m_omega;
  buffer[++counter] = m_eta;
  
  buffer[++counter] = m_velocity.getX();
	buffer[++counter] = m_velocity.getY();
	buffer[++counter] = m_velocity.getZ();
	
  buffer[++counter] = m_vectorP.getX();
	buffer[++counter] = m_vectorP.getY();
	buffer[++counter] = m_vectorP.getZ();
}

//***************************************************************************

void PhaseEulerKorteweg::getBufferSlopes(double* buffer, int& counter)
{
	m_density = buffer[++counter];
  m_omega = buffer[++counter];
  m_eta = buffer[++counter];
  
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);
  
  m_vectorP.setX(buffer[++counter]);
  m_vectorP.setY(buffer[++counter]);
  m_vectorP.setZ(buffer[++counter]);
}

//****************************************************************************
//**************************** VERIFICATION **********************************
//****************************************************************************

void PhaseEulerKorteweg::verifyPhase(const std::string& message) const
{
  if (m_density <= 1.e-10) errors.push_back(Errors(message + "too small density in verifyPhase"));
}

//***************************************************************************

void PhaseEulerKorteweg::verifyAndCorrectPhase()
{
  if (m_density < 1.e-10) m_density = 1.e-10;
}

//***************************************************************************

void PhaseEulerKorteweg::verifyAndCorrectDensityMax()
{
  m_eos->verifyAndCorrectDensityMax(m_density);
}

//****************************************************************************
//**************************** DATA ACCESSORS ********************************
//****************************************************************************

void PhaseEulerKorteweg::setDensity(double density) { m_density = density; }

//***************************************************************************

void PhaseEulerKorteweg::setOmega(const double& omega) { m_omega = omega; }

//***************************************************************************

void PhaseEulerKorteweg::setEta(const double& eta) { m_eta = eta; }

//***************************************************************************

void PhaseEulerKorteweg::setPressure(double pressure) { m_pressure = pressure; }

//***************************************************************************

void PhaseEulerKorteweg::setVelocity(const double& u, const double& v, const double& w) { m_velocity.setXYZ(u, v, w); }

//***************************************************************************

void PhaseEulerKorteweg::setVelocity(const Coord& vit) { m_velocity = vit; }

//***************************************************************************

void PhaseEulerKorteweg::setU(const double& u) { m_velocity.setX(u); }

//***************************************************************************

void PhaseEulerKorteweg::setV(const double& v) { m_velocity.setY(v); }

//***************************************************************************

void PhaseEulerKorteweg::setW(const double& w) { m_velocity.setZ(w); }

//***************************************************************************

void PhaseEulerKorteweg::setVectorP(const double& Px, const double& Py, const double& Pz) { m_vectorP.setXYZ(Px, Py, Pz); }

//***************************************************************************

void PhaseEulerKorteweg::setVectorP(const Coord& vecP) { m_vectorP = vecP; }

//***************************************************************************

void PhaseEulerKorteweg::setVectorPX(const double& Px) { m_vectorP.setX(Px); }

//***************************************************************************

void PhaseEulerKorteweg::setVectorPY(const double& Py) { m_vectorP.setY(Py); }

//***************************************************************************

void PhaseEulerKorteweg::setVectorPZ(const double& Pz) { m_vectorP.setZ(Pz); }

//***************************************************************************

void PhaseEulerKorteweg::setEos(Eos* eos) { m_eos = eos; }

//****************************************************************************
//***************************** OPERATORS ************************************
//****************************************************************************

void PhaseEulerKorteweg::changeSign()
{
  m_density = -m_density;
  m_omega = -m_omega;
  m_eta = -m_eta;
  m_velocity = m_velocity*-1.;
  m_vectorP = m_vectorP*-1;
}

//***************************************************************************

void PhaseEulerKorteweg::multiplyAndAdd(const Phase &slopesPhasesTemp, const double& coeff)
{
  m_density += slopesPhasesTemp.getDensity()*coeff;
  m_omega += slopesPhasesTemp.getOmega()*coeff;
  m_eta +=   slopesPhasesTemp.getEta()*coeff;
  m_velocity += slopesPhasesTemp.getVelocity()*coeff;
  m_vectorP +=  slopesPhasesTemp.getVectorP()*coeff;
}

//***************************************************************************

void PhaseEulerKorteweg::divide(const double& coeff)
{
  m_density /= coeff;
  m_omega /= coeff;
  m_eta /= coeff;
  m_velocity /= coeff;
  m_vectorP /= coeff;
}

//***************************************************************************