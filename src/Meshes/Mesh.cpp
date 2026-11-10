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

#include "Mesh.h"

//***********************************************************************

Mesh::Mesh() : m_numFile(0), m_problemDimension(0) {}

//***********************************************************************

Mesh::~Mesh() {}

//***********************************************************************

void Mesh::writeResultsGnuplot(std::vector<Cell*>* cellsLvl, std::ofstream& fileStream, GeometricObject* objet, bool recordPsat) const
{
  for (unsigned int c = 0; c < cellsLvl[0].size(); c++) {
    if (cellsLvl[0][c]->printGnuplotAMR(fileStream, m_problemDimension, objet, recordPsat)) break;
  }
}

//****************************************************************************
//****************************** Parallele ***********************************
//****************************************************************************

void Mesh::initializePersistentCommunications(const TypeMeshContainer<Cell*>& cells, SimulationOrder schemeOrder)
{
  int numberVariablesPhaseToSend(0);
  for (int k = 0; k < numberPhases; k++) {
    numberVariablesPhaseToSend += cells[0]->getPhase(k)->numberOfTransmittedVariables();
  }
  int numberVariablesMixtureToSend = cells[0]->getMixture()->numberOfTransmittedVariables();
  int m_numberPrimitiveVariables   = numberVariablesPhaseToSend + numberVariablesMixtureToSend + numberTransports;
  int m_numberSlopeVariables(0);
  if (schemeOrder == SimulationOrder::second) {

    int numberSlopesPhaseToSend(0);
    int numberSlopesMixtureToSend(0);
    int numberSlopesTransportToSend(0);

    if (this->getType() != TypeM::UNS) {
      for (int k = 0; k < numberPhases; k++) {
        numberSlopesPhaseToSend += cells[0]->getPhase(k)->numberOfTransmittedSlopes();
      }
      numberSlopesMixtureToSend   = cells[0]->getMixture()->numberOfTransmittedSlopes();
      m_numberSlopeVariables      = 1 + 1; //+1 for the interface detection + 1 for slope index
      numberSlopesTransportToSend = numberTransports;
    }
    else {
      for (int k = 0; k < numberPhases; k++) {
        numberSlopesPhaseToSend += cells[0]->getGradPhase(k)->numberOfTransmittedGradients();
      }
      numberSlopesMixtureToSend = cells[0]->getGradMixture()->numberOfTransmittedGradients();
      for (int t = 0; t < numberTransports; t++) {
        numberSlopesTransportToSend += cells[0]->getGradTransport(t)->numberOfTransmittedGradients();
      }
    }
    m_numberSlopeVariables += numberSlopesPhaseToSend + numberSlopesMixtureToSend + numberSlopesTransportToSend;
  }
  parallel.initializePersistentCommunications(m_numberPrimitiveVariables, m_numberSlopeVariables, numberTransports, m_problemDimension);
}

//***********************************************************************

void Mesh::finalizeParallele(const int& lvlMax) { parallel.finalize(lvlMax); }

//***********************************************************************

template <typename T> void Mesh::getData_body(Cell* cell, std::vector<T>& dataset, const int var, const int phase, const bool isAMR) const
{
  constexpr T eps20{1.e-20};

  if (var > 0) { //We want to get the scalar data
    if (phase >= 0) {
      //phases data
      dataset.push_back(static_cast<T>(cell->getPhase(phase)->returnScalar(var)));
    }
    else if (phase == -1) {
      //mixture data
      dataset.push_back(static_cast<T>(cell->getMixture()->returnScalar(var)));
    }
    else if (phase == -2) {
      T transport = static_cast<T>(cell->getTransport(var - 1).getValue());
      if (transport < eps20) {
        transport = 0.;
      }
      dataset.push_back(static_cast<T>(transport));
    }
    else if (phase == -3) {
      dataset.push_back(static_cast<T>(cell->getXi()));
    }
    else if (phase == -4) {
      dataset.push_back(static_cast<T>(cell->getDensityGradient()));
    }
    else if (phase == -5 && isAMR) {
      dataset.push_back(static_cast<T>(rankCpu));
    }
    else if (phase == -6 && isAMR) {
      dataset.push_back(static_cast<T>(cell->getElement()->getKey().getIndex()));
    }
    else if (phase == -7) {
      // Saturation pressure
      dataset.push_back(static_cast<T>(cell->getPsat()));
    }
    else if (phase == -8) {
      // Turbulent viscosity (nu_T)
      dataset.push_back(static_cast<T>(cell->getQPA(turbulenceNumQPA)->getTurbulentViscosity()));
    }
    else {
      Errors::errorMessage("getData: unknown number of phase: ", phase);
    }
  }
  else {
    //We want to get the vector data
    if (phase >= 0) {
      //Phases data
      dataset.push_back(static_cast<T>(cell->getPhase(phase)->returnVector(-var).getX()));
      dataset.push_back(static_cast<T>(cell->getPhase(phase)->returnVector(-var).getY()));
      dataset.push_back(static_cast<T>(cell->getPhase(phase)->returnVector(-var).getZ()));
    }
    else if (phase == -1) {
      //Mixture data
      dataset.push_back(static_cast<T>(cell->getMixture()->returnVector(-var).getX()));
      dataset.push_back(static_cast<T>(cell->getMixture()->returnVector(-var).getY()));
      dataset.push_back(static_cast<T>(cell->getMixture()->returnVector(-var).getZ()));
    }
    else {
      Errors::errorMessage("getData: unknown number of phase: ", phase);
    }
  } //End vector
}

//***********************************************************************

//! \brief Instanciation of previous function with double type to force compilation of the function body and avoid link error
template void Mesh::getData_body<double>(Cell*, std::vector<double>&, const int, const int, const bool) const;

//***********************************************************************

//! \brief Instanciation of previous function with float type to force compilation of the function body and avoid link error
template void Mesh::getData_body<float>(Cell*, std::vector<float>&, const int, const int, const bool) const;

//***********************************************************************
