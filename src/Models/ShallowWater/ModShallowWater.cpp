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

#include "ModShallowWater.h"

const std::string ModShallowWater::NAME = "SHALLOWWATER";

//****************************************************************************

ModShallowWater::ModShallowWater(const int& numbTransports) : Model(NAME, numbTransports)
{
  fluxBuff = new FluxShallowWater();
  for (int i = 0; i < 3; i++) {
    sourceCons.push_back(new FluxShallowWater());
  }
}

//****************************************************************************

ModShallowWater::~ModShallowWater()
{
  delete fluxBuff;
  for (int i = 0; i < 3; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//****************************************************************************

void ModShallowWater::allocateCons(Flux** cons) { *cons = new FluxShallowWater; }

//***********************************************************************

void ModShallowWater::allocatePhase(Phase** phase) { *phase = new PhaseShallowWater; }

//***********************************************************************

void ModShallowWater::allocateMixture(Mixture** mixture) { *mixture = new MixShallowWater; }

//***********************************************************************

void ModShallowWater::fulfillState(Phase** phases, Mixture* /*mixture*/)
{
  // Compute pressure and sound speed from EoS (velocity is not used but required for sake of genericity).
  phases[0]->extendedCalculusPhase(phases[0]->getVelocity());
}

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModShallowWater::solveRiemannIntern(
  Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double>& /*boundData*/) const
{
  double cL, cR, sL, sR, mL, mR, sM;
  double uL, uR, vL, vR, hL, hR, pL, pR;

  Phase* phaseLeft{cellLeft.getPhase(0)};   // Get phase '0' (our only phase) of left cell
  Phase* phaseRight{cellRight.getPhase(0)}; // Get phase '0' (our only phase) of right cell

  // Remark: phase variables have already been projected in the normal plane
  // Get left state
  uL = phaseLeft->getU();
  vL = phaseLeft->getV();
  hL = phaseLeft->getHeight();
  pL = phaseLeft->getPressure();
  cL = phaseLeft->getSoundSpeed();

  // Get right state
  uR = phaseRight->getU();
  vR = phaseRight->getV();
  hR = phaseRight->getHeight();
  pR = phaseRight->getPressure();
  cR = phaseRight->getSoundSpeed();

  // TO FILL
  assert(0 && "Compute fastest left and right characteristic speed.");
  // Compute the fastest left characteristic speed
  //sL =

  // Compute the fastest right characteristic speed
  //sR =

  // Uncomment computation of maximal time step
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  if (std::fabs(sR) > 1.e-3) dtMax = std::min(dtMax, dxRight / std::fabs(sR));

  // TO FILL

  // Compute left and right mass flow rates
  assert(0 && "Compute fastest left and right mass flow rates.");
  //mL =
  //mR =

  // Compute the velocity of middle wave
  assert(0 && "Compute the velocity of middle wave.");
  //sM =

  if (sL > 0.) {
    // TO FILL
    assert(0 && "Compute left flux.");

    // Use left flux
    // static_cast<FluxShallowWater*>(fluxBuff)->m_mass =
    // static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX( ? );
    // static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY( ? );
  }
  else if (sR < 0.) {
    // TO FILL
    assert(0 && "Compute right flux.");

    // Use right flux
    // static_cast<FluxShallowWater*>(fluxBuff)->m_mass =
    // static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX( ? );
    // static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY( ? );
  }
  // HLLC
  else if (sM > 0.) {
    // TO FILL
    assert(0 && "Compute p* and h* from left state.");

    // Compute pStar and hStar from left state
    // double pStar =
    // double hStar =

    // Compute flux from pStar, hStar and left state
    // static_cast<FluxShallowWater*>(fluxBuff)->m_mass =
    // static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX(?);
    // static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY(?);
  }
  else {
    // TO FILL
    assert(0 && "Compute p* and h* from right state.");

    // Compute pStar and hStar from right state
    // double pStar =
    // double hStar =

    // Compute flux from pStar, hStar and right state
    // static_cast<FluxShallowWater*>(fluxBuff)->m_mass =
    // static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX(?);
    // static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY(?);
  }

  // Contact discontinuity velocity
  static_cast<FluxShallowWater*>(fluxBuff)->m_sM = sM;
  // Height of left and right states
  static_cast<FluxShallowWater*>(fluxBuff)->m_hL = hL;
  static_cast<FluxShallowWater*>(fluxBuff)->m_hR = hR;
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions ****************
//****************************************************************************

void ModShallowWater::solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double>& /*boundData*/) const
{
  double cL, sL;
  double uL, pL, hL;
  double pStar{0.};

  // By convention, at bondary cells, the state inside the cell is considered as the left state.
  // Get left state (so internal state)
  Phase* phaseLeft{cellLeft.getPhase(0)};
  uL = phaseLeft->getU();
  pL = phaseLeft->getPressure();
  hL = phaseLeft->getHeight();
  cL = phaseLeft->getSoundSpeed();

  // TO FILL
  assert(0 && "Compute the left characteristic speed.");
  // Compute the left characteristic speed
  //sL =

  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));

  // TO FILL
  assert(0 && "Compute p*.");
  // Compute p*
  // pStar =

  // Compute Fluxes
  assert(0 && "Compute Fluxes.");
  // static_cast<FluxShallowWater*>(fluxBuff)->m_mass =
  // static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX(?);
  // static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY(?);

  // Middle wave velocity (null at wall)
  // static_cast<FluxShallowWater*>(fluxBuff)->m_sM =

  // Water height for left state
  static_cast<FluxShallowWater*>(fluxBuff)->m_hL = hL;
}

//****************************************************************************
//******************************* Accessors **********************************
//****************************************************************************

double ModShallowWater::selectScalar(Phase** phases, Mixture* /*mixture*/, Transport* /*transports*/, Variable nameVariable, int /*num*/) const
{
  switch (nameVariable) {
  case Variable::pressure:
    return phases[0]->getPressure();
    break;
  case Variable::height:
    return phases[0]->getHeight();
    break;
  case Variable::velocityU:
    return phases[0]->getU();
    break;
  case Variable::velocityV:
    return phases[0]->getV();
    break;
  case Variable::velocityMag:
    return phases[0]->getVelocity().norm();
    break;
  default:
    Errors::errorMessage("nameVariable unknown in selectScalar.");
    return 0;
    break;
  }
}

//****************************************************************************
//***************************** Others methods *******************************
//****************************************************************************

void ModShallowWater::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
