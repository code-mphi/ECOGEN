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

#include "BoundCondMovingWall.h"

using namespace tinyxml2;

//****************************************************************************

BoundCondMovingWall::BoundCondMovingWall(const BoundCondMovingWall& Source, const int& lvl) : BoundCondWall(Source, lvl)
{
  m_amplitude       = Source.m_amplitude;
  m_frequency       = Source.m_frequency;
  m_frequencyOffset = Source.m_frequencyOffset;
}

//****************************************************************************

BoundCondMovingWall::BoundCondMovingWall(int numPhysique, XMLElement* element, std::string fileName) :
  BoundCondWall(numPhysique), m_amplitude(0.), m_frequency(0.), m_frequencyOffset(0.)
{
  XMLElement* subElement(element->FirstChildElement("dataMovingWall"));
  if (subElement != NULL) {
    XMLError error;
    XMLElement* elementPhysWall;

    //Moving wall amplitude and frequency
    error = subElement->QueryDoubleAttribute("amplitude", &m_amplitude);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("amplitude", fileName, __FILE__, __LINE__);
    error = subElement->QueryDoubleAttribute("frequency", &m_frequency);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("frequency", fileName, __FILE__, __LINE__);
    error = subElement->QueryDoubleAttribute("frequencyOffset", &m_frequencyOffset);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("frequencyOffset", fileName, __FILE__, __LINE__);

    // Heat transfer
    elementPhysWall = subElement->FirstChildElement("dataWallHeatTransfer");
    if (elementPhysWall != NULL) {
      std::string heatCondition(elementPhysWall->Attribute("heatCondition"));
      Tools::uppercase(heatCondition);
      // One could use wall with imposed temperature, imposed flux density or adiabatic (default)
      // This option requires the conductivity additionnal physic
      if (heatCondition == "TEMPERATURE") {
        m_heatCondition = TypeBCHeat::IMPOSEDTEMP;
        error           = elementPhysWall->QueryDoubleAttribute("temperature", &m_imposedHeatQuantity);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("temperature", fileName, __FILE__, __LINE__);
      }
      else if (heatCondition == "FLUX") {
        m_heatCondition = TypeBCHeat::IMPOSEDFLUX;
        error           = elementPhysWall->QueryDoubleAttribute("flux", &m_imposedHeatQuantity);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("flux", fileName, __FILE__, __LINE__);
      }
    }

    // MRF wall
    elementPhysWall = subElement->FirstChildElement("dataWallMRF");
    if (elementPhysWall != NULL) {
      double omegaX(0.), omegaY(0.), omegaZ(0.);
      error = elementPhysWall->QueryDoubleAttribute("omegaX", &omegaX);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("omegaX", fileName, __FILE__, __LINE__);
      error = elementPhysWall->QueryDoubleAttribute("omegaY", &omegaY);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("omegaY", fileName, __FILE__, __LINE__);
      error = elementPhysWall->QueryDoubleAttribute("omegaZ", &omegaZ);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("omegaZ", fileName, __FILE__, __LINE__);
      m_omegaWall.setXYZ(omegaX, omegaY, omegaZ);
      m_isMRFWall = true;
    }
  }
}

//****************************************************************************

BoundCondMovingWall::BoundCondMovingWall(int numPhysique) : BoundCondWall(numPhysique), m_amplitude(0.), m_frequency(0.), m_frequencyOffset(0.) {}

//****************************************************************************

BoundCondMovingWall::~BoundCondMovingWall() {}

//****************************************************************************

void BoundCondMovingWall::createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondMovingWall(*(this)));
}

//****************************************************************************

void BoundCondMovingWall::solveRiemannBoundary(Cell& cellLeft, const double& dxLeft, double& dtMax)
{
  //Compute moving wall velocity in the wall normal direction
  //Careful, if moving wall on the left, bottom or back, the sign of the velocity is the opposite
  //if (TB->physicalTime > 1.e-6) m_amplitude = -1.;
  //if (TB->physicalTime > 2.e-6) m_amplitude = 1.;
  double wallVelocity = m_amplitude * std::sin(2. * M_PI * m_frequency * TB->physicalTime + m_frequencyOffset);
  model->solveRiemannMovingWall(cellLeft, dxLeft, dtMax, wallVelocity, m_boundData);
}

//****************************************************************************

void BoundCondMovingWall::solveRiemannTransportBoundary(Cell& cellLeft) const
{
  //Compute moving wall velocity in the wall normal direction
  //Careful, if moving wall on the left, bottom or back, the sign of the velocity is the opposite
  //if (TB->physicalTime > 1.e-6) m_amplitude = -1.;
  //if (TB->physicalTime > 2.e-6) m_amplitude = 1.;
  double wallVelocity = m_amplitude * std::sin(2. * M_PI * m_frequency * TB->physicalTime + m_frequencyOffset);
  model->solveRiemannTransportMovingWall(cellLeft, wallVelocity);
}

//****************************************************************************
//*******************************Accessors************************************
//****************************************************************************

double BoundCondMovingWall::getWallVelocity()
{
  //Compute moving wall velocity in the wall normal direction
  //Careful, if moving wall on the left, bottom or back, the sign of the velocity is the opposite
  //if (TB->physicalTime > 1.e-6) m_amplitude = -1.;
  //if (TB->physicalTime > 2.e-6) m_amplitude = 1.;
  return m_amplitude * std::sin(2. * M_PI * m_frequency * TB->physicalTime + m_frequencyOffset);
}

//****************************************************************************
//*******************************AMR Method***********************************
//****************************************************************************

void BoundCondMovingWall::creerCellInterfaceChild() { m_cellInterfacesChildren.push_back(new BoundCondMovingWall(*this, m_lvl + 1)); }

//****************************************************************************
