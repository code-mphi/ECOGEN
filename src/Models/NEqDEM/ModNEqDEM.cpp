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

#include "ModNEqDEM.h"

const std::string ModNEqDEM::NAME = "NONEQ";

//***********************************************************************

ModNEqDEM::ModNEqDEM(const int& numbTransports, const int& numbPhases, const bool& HLLCRiemannSolver, TypeContactAreas typeContactAreas) :
  ModNEq(NAME, numbTransports)
{
  fluxBuff    = new FluxNEqDEM(numbPhases);
  fluxBuffMRF = new FluxNEqDEM(numbPhases);
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxNEqDEM(numbPhases));
  }
  if (HLLCRiemannSolver) {
    m_RiemannSolver = new RiemannSolverHLLC();
  }
  else {
    m_RiemannSolver = new RiemannSolverExact();
  }
  //Type of contact areas: 1 = symmetric, 2 = nonsymmetric
  if (typeContactAreas == TypeContactAreas::symmetric) {
    contactAreas    = new SymmetricContactAreas(numbPhases);
    contactAreas_02 = new SymmetricContactAreas(numbPhases);
  }
  else if (typeContactAreas == TypeContactAreas::nonsymmetric) {
    contactAreas    = new NonSymmetricContactAreas(numbPhases);
    contactAreas_02 = new NonSymmetricContactAreas(numbPhases);
  }
  m_phaseSlope.resize(numbPhases);
  m_alphaLeft.resize(numbPhases);
  m_alphaRight.resize(numbPhases);
}

//***********************************************************************

ModNEqDEM::ModNEqDEM(const std::string& name, const int& numbTransports) : ModNEq(name, numbTransports) {}

//***********************************************************************

ModNEqDEM::~ModNEqDEM()
{
  // Delete operators are called in ModNEq destructor
  delete contactAreas;
  delete contactAreas_02;
  delete m_RiemannSolver;
}

//***********************************************************************

void ModNEqDEM::allocateCons(Flux** cons) { *cons = new FluxNEqDEM(numberPhases); }

//***********************************************************************

void ModNEqDEM::setPhaseSlope(Phase& phaseSlope, int k) { m_phaseSlope[k].copyPhase(phaseSlope); }

//***********************************************************************

void ModNEqDEM::setSlopeDistance(const double& slopeDistance) { m_slopeDistance = slopeDistance; }

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModNEqDEM::solveRiemannIntern(
  Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double>& /*boundData*/) const
{
  //Initialise alpha left and right
  //-------------------------------
  for (int i = 0; i < numberPhases; i++) {
    m_alphaLeft[i]  = cellLeft.getPhase(i)->getAlpha();
    m_alphaRight[i] = cellRight.getPhase(i)->getAlpha();
  }

  //Computation of contact areas
  //----------------------------
  //Classical two-phase computation
  // contactAreas->m_contacts[0][0] = std::min(cellLeft.getPhase(0)->getAlpha(), cellRight.getPhase(0)->getAlpha());
  // contactAreas->m_contacts[0][1] = std::max(0., cellLeft.getPhase(0)->getAlpha() - cellRight.getPhase(0)->getAlpha());
  // contactAreas->m_contacts[1][0] = std::max(0., cellRight.getPhase(0)->getAlpha() - cellLeft.getPhase(0)->getAlpha());
  // contactAreas->m_contacts[1][1] = std::min(cellLeft.getPhase(1)->getAlpha(), cellRight.getPhase(1)->getAlpha());
  // contactAreas->recalibrateContactAreas();

  double epsilon(1.e-12);

  // Reduce oscillations problem for more than 2 phases and 2nd order:
  if (numberPhases > 2 && simulationOrder == SimulationOrder::second) {
    double one_L(0.), one_R(0.);
    for (int i = 0; i < numberPhases; i++) {
      if (m_alphaLeft[i] < epsilon) m_alphaLeft[i] = 0.;
      one_L += m_alphaLeft[i];
      if (m_alphaRight[i] < epsilon) m_alphaRight[i] = 0;
      ;
      one_R += m_alphaRight[i];
    }
    for (int i = 0; i < numberPhases; i++) {
      m_alphaLeft[i]  = m_alphaLeft[i] / one_L;
      m_alphaRight[i] = m_alphaRight[i] / one_R;
    }
  }

  //N-phase computation
  contactAreas->computeContactAreas(m_alphaLeft, m_alphaRight);

  //Fluxes for each contact areas
  //------------------------------

  //Initialize the flux vector to zero (for summation)
  fluxBuff->setToZero();

  for (int i = 0; i < numberPhases; i++) {
    for (int j = 0; j < numberPhases; j++) {

      if (contactAreas->m_contacts[i][j] < epsilon) continue;

      //Left and right states for Riemann problem
      //-----------------------------------------
      double pStar(0.);

      this->m_RiemannSolver->solver(*(cellLeft.getPhase(i)), *(cellRight.getPhase(j)), dxLeft, dxRight, dtMax, pStar, m_fluxEuler);

      static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j] = m_fluxEuler.m_sM;
      static_cast<FluxNEqDEM*>(fluxBuff)->m_pI[i][j] = pStar;

      //Convective Flux
      if (i == j || (i != j && m_fluxEuler.m_sM >= 0.)) { //Otherwise it will be taken into account for the couple j-i
        static_cast<FluxNEqDEM*>(fluxBuff)->m_mass[i]     += contactAreas->m_contacts[i][j] * m_fluxEuler.m_mass;
        static_cast<FluxNEqDEM*>(fluxBuff)->m_momentum[i] += contactAreas->m_contacts[i][j] * m_fluxEuler.m_momentum;
        static_cast<FluxNEqDEM*>(fluxBuff)->m_energ[i]    += contactAreas->m_contacts[i][j] * m_fluxEuler.m_energ;
      }
      else if (i != j && m_fluxEuler.m_sM < 0.) {
        static_cast<FluxNEqDEM*>(fluxBuff)->m_mass[j]     += contactAreas->m_contacts[i][j] * m_fluxEuler.m_mass;
        static_cast<FluxNEqDEM*>(fluxBuff)->m_momentum[j] += contactAreas->m_contacts[i][j] * m_fluxEuler.m_momentum;
        static_cast<FluxNEqDEM*>(fluxBuff)->m_energ[j]    += contactAreas->m_contacts[i][j] * m_fluxEuler.m_energ;
      }
    }
  }

  //Second-order computation
  //------------------------
  if (simulationOrder == SimulationOrder::second) {
    double one(0.);
    // Compute left-extrapolated volume fractions of the cell left, stored in the "right" alpha (since the original is not used anymore)
    for (int i = 0; i < numberPhases; i++) {
      m_alphaRight[i]  = m_alphaLeft[i] - m_phaseSlope[i].getAlpha() * m_slopeDistance;
      one             += m_alphaRight[i];
    }
    for (int i = 0; i < numberPhases; i++) {
      m_alphaRight[i] /= one;
    }

    contactAreas_02->computeContactAreas(m_alphaRight, m_alphaLeft);

    // Riemann problem for the non-conservative flux internal to the cell,
    // only for 2nd order and if the left cell is of type CellO2 (because it contains the slopes needed for the 2nd order flux)
    double pStar(0.);
    for (int i = 0; i < numberPhases; i++) { //Reconstruct the extrapolated left states from the left cell
      m_phaseSlope[i].setDensity(cellLeft.getPhase(i)->getDensity() - m_phaseSlope[i].getDensity() * m_slopeDistance);
      m_phaseSlope[i].setPressure(cellLeft.getPhase(i)->getPressure() - m_phaseSlope[i].getPressure() * m_slopeDistance);
      m_phaseSlope[i].setU(cellLeft.getPhase(i)->getU() - m_phaseSlope[i].getU() * m_slopeDistance);
      m_phaseSlope[i].setV(cellLeft.getPhase(i)->getV() - m_phaseSlope[i].getV() * m_slopeDistance);
      m_phaseSlope[i].setW(cellLeft.getPhase(i)->getW() - m_phaseSlope[i].getW() * m_slopeDistance);
    }
    // for (int i = 0; i < numberPhases; i++) { //Reconstruct the cell-centered states from the left cell
    //   m_phaseSlope[i].setDensity(cellLeft.getPhase(i)->getDensity() - m_phaseSlope[i].getDensity() * dxLeft * 0.5);
    //   m_phaseSlope[i].setPressure(cellLeft.getPhase(i)->getPressure() - m_phaseSlope[i].getPressure() * dxLeft * 0.5);
    //   m_phaseSlope[i].setVelocity(cellLeft.getPhase(i)->getVelocity() - m_phaseSlope[i].getVelocity() * dxLeft * 0.5);
    // }
    for (int i = 0; i < numberPhases; i++) {
      for (int j = 0; j < numberPhases; j++) {
        if (i != j && contactAreas_02->m_contacts[i][j] > epsilon) {
          // this->m_RiemannSolver->solver(*(cellLeft.getPhase(i)), *(cellLeft.getPhase(j)), dxLeft, dxLeft, dtMax, pStar, m_fluxEuler); //Riemann problem between the two right extrapolated states within the left cell
          this->m_RiemannSolver->solver(m_phaseSlope[i],
                                        *(cellLeft.getPhase(j)),
                                        dxLeft,
                                        dxLeft,
                                        dtMax,
                                        pStar,
                                        m_fluxEuler); //Riemann problem between both extrapolated states within the left cell
          // this->m_RiemannSolver->solver(m_phaseSlope[i], m_phaseSlope[j], dxLeft, dxLeft, dtMax, pStar, m_fluxEuler); //Riemann problem between both cell-centered states within the left cell
          static_cast<FluxNEqDEM*>(fluxBuff)->m_uI_02[i][j] = m_fluxEuler.m_sM;
          static_cast<FluxNEqDEM*>(fluxBuff)->m_pI_02[i][j] = pStar;
        }
      }
    }
  }
}
//***********************************************************************

void ModNEqDEM::solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double>& /*boundData*/) const
{
  //Initialise alpha left and right
  //-------------------------------
  for (int i = 0; i < numberPhases; i++) {
    m_alphaLeft[i] = cellLeft.getPhase(i)->getAlpha();
  }

  //Computation of contact areas
  //----------------------------
  contactAreas->computeContactAreas(m_alphaLeft, m_alphaLeft);

  //Fluxes for each contact length
  //------------------------------

  //Initialize the flux vector to zero (for summation)
  fluxBuff->setToZero();

  double epsilon(1.e-12);
  for (int i = 0; i < numberPhases; i++) {
    for (int j = 0; j < numberPhases; j++) {

      if (contactAreas->m_contacts[i][j] < epsilon) continue;

      //Left and right states for Riemann problem
      //-----------------------------------------
      double pStar{0.};
      this->m_RiemannSolver->solverWall(*(cellLeft.getPhase(i)), dxLeft, dtMax, pStar, m_fluxEuler);

      static_cast<FluxNEqDEM*>(fluxBuff)->m_uI[i][j] = 0.;
      static_cast<FluxNEqDEM*>(fluxBuff)->m_pI[i][j] = pStar;

      //Convective Flux, only for i=j couples
      static_cast<FluxNEqDEM*>(fluxBuff)->m_mass[i]     += contactAreas->m_contacts[i][j] * m_fluxEuler.m_mass;
      static_cast<FluxNEqDEM*>(fluxBuff)->m_momentum[i] += contactAreas->m_contacts[i][j] * m_fluxEuler.m_momentum;
      static_cast<FluxNEqDEM*>(fluxBuff)->m_energ[i]    += contactAreas->m_contacts[i][j] * m_fluxEuler.m_energ;
    }
  }
}

//****************************************************************************
//*** Half Riemann solver for MRF interface between static/rotating region ***
//****************************************************************************

//FP//TODO

//****************************************************************************
//********************** Transport Riemann solvers ***************************
//****************************************************************************

//FP//TODO

//****************************************************************************
