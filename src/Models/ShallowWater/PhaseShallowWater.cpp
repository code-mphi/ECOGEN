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

#include "PhaseShallowWater.h"

using namespace tinyxml2;

//***************************************************************************

PhaseShallowWater::PhaseShallowWater() : m_height(0.), m_pressure(0.), m_soundSpeed(0.), m_eos(0) { m_velocity.setXYZ(0., 0., 0.); }

//***************************************************************************

PhaseShallowWater::PhaseShallowWater(XMLElement* material, Eos* eos, std::string fileName) : m_soundSpeed(0.), m_eos(eos)
{
  XMLElement* sousElement(material->FirstChildElement("dataFluid"));
  if (sousElement == NULL) throw ErrorXMLElement("dataFluid", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;

  //height: return an error if keyword height (and value) is not found in input file. Assign readed value otherwise.
  if (sousElement->QueryDoubleAttribute("height", &m_height) != XML_NO_ERROR) throw ErrorXMLElement("height", fileName, __FILE__, __LINE__);

  //velocity
  XMLElement* velocity(sousElement->FirstChildElement("velocity"));
  if (velocity == NULL) throw ErrorXMLElement("velocity", fileName, __FILE__, __LINE__);
  double velocityX(0.), velocityY(0.), velocityZ(0.);
  // Read x-component of velocity and trhow and error if not found
  error = velocity->QueryDoubleAttribute("x", &velocityX);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName, __FILE__, __LINE__);
  // Read y-component of velocity and trhow and error if not found
  error = velocity->QueryDoubleAttribute("y", &velocityY);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName, __FILE__, __LINE__);
  // z-component of velocity is unexpected in this model: throw an error if it is found
  error = velocity->QueryDoubleAttribute("z", &velocityZ);
  if (error == XML_NO_ERROR) throw ErrorXMLAttribut("Input z-component of velocity is unexpected for this model", fileName, __FILE__, __LINE__);
}

//***************************************************************************

PhaseShallowWater::~PhaseShallowWater() {}

//***************************************************************************

void PhaseShallowWater::allocateAndCopyPhase(Phase** vecPhase) { *vecPhase = new PhaseShallowWater(*this); }

//***************************************************************************

void PhaseShallowWater::copyPhase(Phase& phase)
{
  m_height     = phase.getHeight();
  m_velocity   = phase.getVelocity();
  m_pressure   = phase.getPressure();
  m_eos        = phase.getEos();
  m_soundSpeed = phase.getSoundSpeed();
}

//***************************************************************************

void PhaseShallowWater::extendedCalculusPhase(const Coord& /*velocity*/)
{
  m_pressure   = m_eos->computePressure(m_height);
  m_soundSpeed = m_eos->computeSoundSpeed(m_height);
}

//***************************************************************************

void PhaseShallowWater::localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void PhaseShallowWater::reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
//****************************** DATA PRINTING *******************************
//****************************************************************************

double PhaseShallowWater::returnScalar(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return m_height;
    break;
  case 2:
    return m_pressure;
    break;
  default:
    return 0.;
    break;
  }
}

//***************************************************************************

Coord PhaseShallowWater::returnVector(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return m_velocity;
    break;
  default:
    return 0;
    break;
  }
}

//***************************************************************************

std::string PhaseShallowWater::returnNameScalar(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return "Height";
    break;
  case 2:
    return "Pressure";
    break;
  default:
    return "NoName";
    break;
  }
}

//***************************************************************************

std::string PhaseShallowWater::returnNameVector(const int& numVar) const
{
  switch (numVar) {
  case 1:
    return "Velocity";
    break;
  default:
    return "NoName";
    break;
  }
}

//****************************************************************************
//************************* READING FROM FILE ********************************
//****************************************************************************

void PhaseShallowWater::setScalar(const int& numVar, const double& value)
{
  switch (numVar) {
  case 1:
    m_height = value;
    break;
  default:
    Errors::errorMessage("numVar not found in PhaseShallowWater::setScalar");
    break;
  }
}

//****************************************************************************

void PhaseShallowWater::setVector(const int& numVar, const Coord& value)
{
  switch (numVar) {
  case 1:
    m_velocity = value;
    break;
  default:
    Errors::errorMessage("numVar not found in PhaseShallowWater::setVector");
    break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int PhaseShallowWater::numberOfTransmittedVariables() const
{
  //3 variables (height and two components of velocity) + eos
  return 4;
}

//***************************************************************************

void PhaseShallowWater::fillBuffer(double* buffer, int& counter) const
{
  buffer[++counter] = m_height;
  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = static_cast<double>(m_eos->getNumber());
}

//***************************************************************************

void PhaseShallowWater::fillBuffer(std::vector<double>& dataToSend) const
{
  dataToSend.push_back(m_height);
  dataToSend.push_back(m_velocity.getX());
  dataToSend.push_back(m_velocity.getY());
  dataToSend.push_back(static_cast<double>(m_eos->getNumber()));
}

//***************************************************************************

void PhaseShallowWater::getBuffer(double* buffer, int& counter, Eos** eos)
{
  m_height = buffer[++counter];
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_eos = eos[static_cast<int>(buffer[++counter])];
}

//***************************************************************************

void PhaseShallowWater::getBuffer(std::vector<double>& dataToReceive, int& counter, Eos** eos)
{
  m_height = dataToReceive[counter++];
  m_velocity.setX(dataToReceive[counter++]);
  m_velocity.setY(dataToReceive[counter++]);
  m_eos = eos[static_cast<int>(dataToReceive[counter++])];
}

//****************************************************************************
//**************************** DATA ACCESSORS ********************************
//****************************************************************************

void PhaseShallowWater::setHeight(const double& height) { m_height = height; }

//***************************************************************************

void PhaseShallowWater::setPressure(double pressure) { m_pressure = pressure; }

//***************************************************************************

void PhaseShallowWater::setVelocity(const double& u, const double& v, const double& w) { m_velocity.setXYZ(u, v, w); }

//***************************************************************************

void PhaseShallowWater::setU(const double& u) { m_velocity.setX(u); }

//***************************************************************************

void PhaseShallowWater::setV(const double& v) { m_velocity.setY(v); }

//***************************************************************************

void PhaseShallowWater::setEos(Eos* eos) { m_eos = eos; }

//***************************************************************************

void PhaseShallowWater::setSoundSpeed(double soundSpeed) { m_soundSpeed = soundSpeed; }
