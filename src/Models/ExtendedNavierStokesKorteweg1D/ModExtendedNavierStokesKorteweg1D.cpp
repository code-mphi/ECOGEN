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

#include "ModExtendedNavierStokesKorteweg1D.h"
#include "PhaseExtendedNavierStokesKorteweg1D.h"

const std::string ModExtendedNavierStokesKorteweg1D::NAME = "EXTENDEDNAVIERSTOKESKORTEWEG1D";

//****************************************************************************

ModExtendedNavierStokesKorteweg1D::ModExtendedNavierStokesKorteweg1D(const int& numbTransports) : Model(NAME, numbTransports)
{
  fluxBuff = new FluxExtendedNavierStokesKorteweg1D();
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxExtendedNavierStokesKorteweg1D());
  }
}

//***********************************************************************

ModExtendedNavierStokesKorteweg1D::ModExtendedNavierStokesKorteweg1D(
  const int& numbTransports, const double& lambda, const double& beta, const double& tau, const double& eps) :
  Model(NAME, numbTransports)
{
  fluxBuff = new FluxExtendedNavierStokesKorteweg1D();
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxExtendedNavierStokesKorteweg1D());
  }
  lambdaEK = lambda;
  betaEK   = beta;
  tauEK    = tau;
  epsEK    = eps;
}

//***********************************************************************

ModExtendedNavierStokesKorteweg1D::ModExtendedNavierStokesKorteweg1D(
  const std::string& name, const int& numbTransports, const double& lambda, const double& beta, const double& tau, const double& eps) :
  Model(name, numbTransports)
{
  lambdaEK = lambda;
  betaEK   = beta;
  tauEK    = tau;
  epsEK    = eps;
}

//****************************************************************************

ModExtendedNavierStokesKorteweg1D::~ModExtendedNavierStokesKorteweg1D()
{
  delete fluxBuff;
  for (int i = 0; i < 4; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//****************************************************************************

void ModExtendedNavierStokesKorteweg1D::allocateCons(Flux** cons) { *cons = new FluxExtendedNavierStokesKorteweg1D; }

//***********************************************************************

void ModExtendedNavierStokesKorteweg1D::allocatePhase(Phase** phase) { *phase = new PhaseExtendedNavierStokesKorteweg1D; }

//***********************************************************************

void ModExtendedNavierStokesKorteweg1D::allocateMixture(Mixture** mixture) { *mixture = new MixExtendedNavierStokesKorteweg1D; }

//***********************************************************************

void ModExtendedNavierStokesKorteweg1D::fulfillState(Phase** phases, Mixture* /*mixture*/)
{
  phases[0]->extendedCalculusPhase(phases[0]->getVelocity());
}

//***********************************************************************

void ModExtendedNavierStokesKorteweg1D::initializeAugmentedVariables(Cell* cell)
{
  Phase* phase(cell->getPhase(0));

  //Density sur Eta = 1 at t = 0
  phase->setDensitySurEta(1.);
  //Omega
  cell->computeGradients(gradVelocityX, variableVelocity, numeratorDefault);
  double omega(0.);
  omega = -phase->getDensity() * gradVelocityX[0].getX(); //QD//TODO replace U_x by divergence of U = U_x + U_y + U_z or Trace of (gradient U)
  phase->setOmega(omega);
  //VectorP
  cell->computeGradients(gradRho, variableDensity, numeratorDefault);
  phase->setVectorP(gradRho[0]);
  //TotalEnergy
  double varTotEnergy(0.);
  varTotEnergy = phase->getEnergy() + 0.5 * phase->getVelocity().squaredNorm() + 0.5 * betaEK * omega * omega +
                 0.5 * lambdaEK * (1. - 1. / phase->getDensitySurEta()) * (1. - 1. / phase->getDensitySurEta()) +
                 0.5 * epsEK * epsEK * phase->getEos()->kappa(phase->getDensity()) / phase->getDensity() * phase->getVectorP().squaredNorm();
  phase->setTotalEnergy(varTotEnergy);
}

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModExtendedNavierStokesKorteweg1D::solveRiemannIntern(
  Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double>& /*boundData*/) const
{
  //Get variables
  double rhoL, rhoR, omegaL, omegaR, rhoSurEtaL, rhoSurEtaR;
  double uL, uR;         // vL, vR;
  double vecPxL, vecPxR; // vecPyL, vecPyR;
  double pressureL, pressureR;
  double EL, ER;
  double kappaL, kappaR;
  double kappaPrimeL, kappaPrimeR;
  double kappaSecondL, kappaSecondR;
  double c2L, c2R;
  double sL(0.), sR(0.);

  Phase *phaseLeft(0), *phaseRight(0);
  phaseLeft  = cellLeft.getPhase(0);
  phaseRight = cellRight.getPhase(0);

  rhoL       = phaseLeft->getDensity();
  omegaL     = phaseLeft->getOmega();
  rhoSurEtaL = phaseLeft->getDensitySurEta();
  uL         = phaseLeft->getU();
  // vL = phaseLeft->getV();
  // wL = phaseLeft->getW();
  vecPxL = phaseLeft->getVectorPX();
  // vecPyL = phaseLeft->getVectorPY();
  // vecPzL = phaseLeft->getVectorPZ();
  pressureL    = phaseLeft->getPressure();
  kappaL       = phaseLeft->getEos()->kappa(rhoL);
  kappaPrimeL  = phaseLeft->getEos()->kappaPrime(rhoL);
  kappaSecondL = phaseLeft->getEos()->kappaSecond(rhoL);
  c2L          = phaseLeft->getSquareSoundSpeed();
  EL           = phaseLeft->getTotalEnergy();

  rhoR       = phaseRight->getDensity();
  omegaR     = phaseRight->getOmega();
  rhoSurEtaR = phaseRight->getDensitySurEta();
  uR         = phaseRight->getU();
  // vR = phaseRight->getV();
  // wR = phaseRight->getW();
  vecPxR = phaseRight->getVectorPX();
  // vecPyR = phaseRight->getVectorPY();
  // vecPzR = phaseRight->getVectorPZ();
  pressureR    = phaseRight->getPressure();
  kappaR       = phaseRight->getEos()->kappa(rhoR);
  kappaPrimeR  = phaseRight->getEos()->kappaPrime(rhoR);
  kappaSecondR = phaseRight->getEos()->kappaSecond(rhoR);
  c2R          = phaseRight->getSquareSoundSpeed();
  ER           = phaseRight->getTotalEnergy();

  //HLLC approximation solver
  double alphaEK = epsEK * epsEK;

  //Left and right tensor P in momentum equation
  double piL = pressureL + 0.5 * alphaEK * (rhoL * kappaPrimeL - kappaL) * phaseLeft->getVectorP().squaredNorm() +
               lambdaEK * rhoL / rhoSurEtaL * (1. - 1. / rhoSurEtaL);
  double tensorPxxL = piL + alphaEK * kappaL * vecPxL * vecPxL;
  //double tensorPyxL= kappaL * vecPxL * vecPyL;
  //double tensorPzxL= kappaL * vecPxL * vecPzL;

  double piR = pressureR + 0.5 * alphaEK * (rhoR * kappaPrimeR - kappaR) * phaseRight->getVectorP().squaredNorm() +
               lambdaEK * rhoR / rhoSurEtaR * (1. - 1. / rhoSurEtaR);
  double tensorPxxR = piR + alphaEK * kappaR * vecPxR * vecPxR;
  //double tensorPyxR= kappaR * vecPxR * vecPyR;
  //double tensorPzxR= kappaR * vecPxR * vecPzR;

  //Compute maximal wave speed using Davis approximation and eigenvalues of the hyperbolic equations
  computeMaxWaveSpeed(
    rhoL, rhoR, uL, uR, rhoSurEtaL, rhoSurEtaR, vecPxL, vecPxR, kappaL, kappaR, kappaPrimeL, kappaPrimeR, kappaSecondL, kappaSecondR, c2L, c2R, sL, sR);
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  if (std::fabs(sR) > 1.e-3) dtMax = std::min(dtMax, dxRight / std::fabs(sR));

  double mL = rhoL * (uL - sL), mR = rhoR * (uR - sR);
  double sM = (tensorPxxR - tensorPxxL + uR * mR - uL * mL) / (mR - mL);
  if (std::fabs(sM) < 1.e-8) sM = 0.;

  double tensorPxxStar  = (mR * tensorPxxL - mL * tensorPxxR + mR * mL * (uL - uR)) / (mR - mL);
  double omegaStar      = (omegaR * mR - omegaL * mL + alphaEK / betaEK * (kappaL * vecPxL - kappaR * vecPxR)) / (mR - mL);
  double KappaVecPxStar = (mR * kappaL * vecPxL - mL * kappaR * vecPxR + betaEK / alphaEK * mR * mL * (omegaR - omegaL)) / (mR - mL);

  if (sL > 0.) {
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_mass = rhoL * uL;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_momentum.setX(rhoL * uL * uL + tensorPxxL);
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqDensitySurEta = rhoSurEtaL * uL;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqOmega         = rhoL * uL * omegaL - alphaEK / betaEK * kappaL * vecPxL;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqVectorP.setX(vecPxL * uL - omegaL);
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_energ = (rhoL * EL + tensorPxxL) * uL - alphaEK * vecPxL * omegaL * kappaL;
  }
  else if (sR < 0.) {
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_mass = rhoR * uR;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_momentum.setX(rhoR * uR * uR + tensorPxxR);
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqDensitySurEta = rhoSurEtaR * uR;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqOmega         = rhoR * uR * omegaR - alphaEK / betaEK * kappaR * vecPxR;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqVectorP.setX(vecPxR * uR - omegaR);
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_energ = (rhoR * ER + tensorPxxR) * uR - alphaEK * vecPxR * omegaR * kappaR;
  }
  else if (sM >= 0.) {
    double rhoStar       = mL / (sM - sL);
    double rhoSurEtaStar = rhoStar / rhoL * rhoSurEtaL;
    double vecPxStar     = (vecPxL * (uL - sL) + omegaStar - omegaL) / (sM - sL);
    double EStar         = EL + (tensorPxxL * uL - tensorPxxStar * sM + alphaEK * (KappaVecPxStar * omegaStar - kappaL * vecPxL * omegaL)) / (mL);

    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_mass = rhoStar * sM;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_momentum.setX(rhoStar * sM * sM + tensorPxxStar);
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqDensitySurEta = rhoSurEtaStar * sM;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqOmega         = rhoStar * sM * omegaStar - alphaEK / betaEK * KappaVecPxStar;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqVectorP.setX(vecPxStar * sM - omegaStar);
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_energ = (rhoStar * EStar + tensorPxxStar) * sM -
                                                                          alphaEK * KappaVecPxStar * omegaStar;
  }
  else {
    double rhoStar       = mR / (sM - sR);
    double rhoSurEtaStar = rhoStar / rhoR * rhoSurEtaR;
    double vecPxStar     = (vecPxR * (uR - sR) + omegaStar - omegaR) / (sM - sR);
    double EStar         = ER + (tensorPxxR * uR - tensorPxxStar * sM + alphaEK * (KappaVecPxStar * omegaStar - kappaR * vecPxR * omegaR)) / (mR);

    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_mass = rhoStar * sM;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_momentum.setX(rhoStar * sM * sM + tensorPxxStar);
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqDensitySurEta = rhoSurEtaStar * sM;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqOmega         = rhoStar * sM * omegaStar - alphaEK / betaEK * KappaVecPxStar;
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqVectorP.setX(vecPxStar * sM - omegaStar);
    static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_energ = (rhoStar * EStar + tensorPxxStar) * sM -
                                                                          alphaEK * KappaVecPxStar * omegaStar;
  }

  //Contact discontinuity velocity
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_sM = sM;
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions** **************
//****************************************************************************

void ModExtendedNavierStokesKorteweg1D::solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double>& /*boundData*/) const
{
  //Get variables
  double rhoL, omegaL, rhoSurEtaL;
  double uL;     // vL, vR;
  double vecPxL; // vecPyL, vecPyR;
  double pressureL;
  double kappaL;
  double kappaPrimeL;
  double kappaSecondL;
  double c2L;
  double sL(0.), sR(0.);

  Phase* phaseLeft(0);
  phaseLeft = cellLeft.getPhase(0);

  rhoL       = phaseLeft->getDensity();
  omegaL     = phaseLeft->getOmega();
  rhoSurEtaL = phaseLeft->getDensitySurEta();
  uL         = phaseLeft->getU();
  // vL = phaseLeft->getV();
  // wL = phaseLeft->getW();
  vecPxL = phaseLeft->getVectorPX();
  // vecPyL = phaseLeft->getVectorPY();
  // vecPzL = phaseLeft->getVectorPZ();
  pressureL    = phaseLeft->getPressure();
  kappaL       = phaseLeft->getEos()->kappa(rhoL);
  kappaPrimeL  = phaseLeft->getEos()->kappaPrime(rhoL);
  kappaSecondL = phaseLeft->getEos()->kappaSecond(rhoL);
  c2L          = phaseLeft->getSquareSoundSpeed();

  //HLLC approximation solver
  double alphaEK = epsEK * epsEK;

  //Left and right tensor P in momentum equation
  double piL = pressureL + 0.5 * alphaEK * (rhoL * kappaPrimeL - kappaL) * phaseLeft->getVectorP().squaredNorm() +
               lambdaEK * rhoL / rhoSurEtaL * (1. - 1. / rhoSurEtaL);
  double tensorPxxL = piL + alphaEK * kappaL * vecPxL * vecPxL;
  //double tensorPyxL= kappaL * vecPxL * vecPyL;
  //double tensorPzxL= kappaL * vecPxL * vecPzL;

  //Compute maximal wave speed using Davis approximation and eigenvalues of the hyperbolic equations
  computeMaxWaveSpeed(
    rhoL, rhoL, uL, -uL, rhoSurEtaL, rhoSurEtaL, vecPxL, -vecPxL, kappaL, kappaL, kappaPrimeL, kappaPrimeL, kappaSecondL, kappaSecondL, c2L, c2L, sL, sR);
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));

  double mL            = rhoL * (uL - sL);
  double tensorPxxStar = tensorPxxL + mL * uL;
  double omegaStar     = omegaL - alphaEK / betaEK * kappaL * vecPxL / mL;

  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_mass = 0.;
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_momentum.setX(tensorPxxStar);
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqDensitySurEta = 0.;
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqOmega         = 0.;
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqVectorP.setX(-omegaStar);
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_energ = 0.;

  //Contact discontinuity velocity
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_sM = 0.;
}

//****************************************************************************

void ModExtendedNavierStokesKorteweg1D::solveRiemannNullFlux() const
{
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_mass = 0.;
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_momentum.setX(0.);
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqDensitySurEta = 0.;
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqOmega         = 0.;
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqVectorP.setX(0.);
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_energ = 0.;
}

//****************************************************************************

void ModExtendedNavierStokesKorteweg1D::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_momentum.reverseProjection(normal, tangent, binormal);
  static_cast<FluxExtendedNavierStokesKorteweg1D*>(fluxBuff)->m_eqVectorP.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************

void ModExtendedNavierStokesKorteweg1D::computeMaxWaveSpeed(const double& rhoL,
                                                            const double& rhoR,
                                                            const double& uL,
                                                            const double& uR,
                                                            const double& rhoSurEtaL,
                                                            const double& rhoSurEtaR,
                                                            const double& vecPxL,
                                                            const double& vecPxR,
                                                            const double& kappaL,
                                                            const double& kappaR,
                                                            const double& kappaPrimeL,
                                                            const double& kappaPrimeR,
                                                            const double& kappaSecondL,
                                                            const double& kappaSecondR,
                                                            const double& c2L,
                                                            const double& c2R,
                                                            double& sL,
                                                            double& sR) const
{
  double alphaEK = epsEK * epsEK;
  //Left
  double aL = c2L + alphaEK * vecPxL * vecPxL * (kappaL / rhoL + 2. * kappaPrimeL + 0.5 * rhoL * kappaSecondL);
  double bL = lambdaEK / rhoSurEtaL / rhoSurEtaL;
  double eL = alphaEK * kappaL / rhoL / betaEK;
  double dL = alphaEK * vecPxL * (rhoL * kappaPrimeL + kappaL) / rhoL;

  double first_sqrtL = std::sqrt((aL + bL - eL) * (aL + bL - eL) + 4. * dL * dL / betaEK);
  double first_termL = aL + eL + bL;

  double s1L = std::sqrt(0.5 * (first_termL + first_sqrtL));
  double s2L = std::sqrt(0.5 * (first_termL - first_sqrtL));

  //Right
  double aR = c2R + alphaEK * vecPxR * vecPxR * (kappaR / rhoR + 2. * kappaPrimeR + 0.5 * rhoR * kappaSecondR);
  double bR = lambdaEK / rhoSurEtaR / rhoSurEtaR;
  double eR = alphaEK * kappaR / rhoR / betaEK;
  double dR = alphaEK * vecPxR * (rhoR * kappaPrimeR + kappaR) / rhoR;

  double first_sqrtR = std::sqrt((aR + bR - eR) * (aR + bR - eR) + 4. * dR * dR / betaEK);
  double first_termR = aR + eR + bR;

  double s1R = std::sqrt(0.5 * (first_termR + first_sqrtR));
  double s2R = std::sqrt(0.5 * (first_termR - first_sqrtR));

  sL = std::min({uL - s1L, uR - s1R, uL - s2L, uR - s2R});
  sR = std::max({uL + s1L, uR + s1R, uL + s2L, uR + s2R});
}

//****************************************************************************
//******************************* Accessors **********************************
//****************************************************************************

double
ModExtendedNavierStokesKorteweg1D::selectScalar(Phase** phases, Mixture* /*mixture*/, Transport* transports, Variable nameVariable, int num) const
{
  switch (nameVariable) {
  case Variable::pressure:
    return phases[0]->getPressure();
    break;
  case Variable::density:
    return phases[0]->getDensity();
    break;
  case Variable::velocityU:
    return phases[0]->getU();
    break;
  case Variable::velocityV:
    return phases[0]->getV();
    break;
  case Variable::velocityW:
    return phases[0]->getW();
    break;
  case Variable::velocityMag:
    return phases[0]->getVelocity().norm();
    break;
  case Variable::transport:
    return transports[num].getValue();
    break;
  default:
    Errors::errorMessage("nameVariable unknown in selectScalar");
    return 0;
    break;
  }
}

//****************************************************************************
