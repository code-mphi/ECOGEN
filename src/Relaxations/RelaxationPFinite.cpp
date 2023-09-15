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

#include "RelaxationPFinite.h"

using namespace tinyxml2;

//Externalized for LSODA solver
double mu;                       //!< Relaxation coefficient. Herein, the relaxation coefficient is identical for all phase_k--phase_j combinations.

//***********************************************************************

RelaxationPFinite::RelaxationPFinite()
{
  m_LSODA = false;
  m_muFactor = 0.;
  mu = 0.;
}

//***********************************************************************

RelaxationPFinite::RelaxationPFinite(XMLElement* element, std::string fileName)
{
  //Read the relaxation-coefficient factor
  XMLError error;
  error = element->QueryDoubleAttribute("rate", &m_muFactor);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("rate", fileName, __FILE__, __LINE__);
  //Initialize the relaxation coefficient
  mu = m_muFactor;
  //Read the used solver
  std::string solver(element->Attribute("solver"));
  Tools::uppercase(solver);
  if      (solver == "EULER") m_LSODA = false;
  else if (solver == "LSODA") m_LSODA = true;
  else throw ErrorXMLAttribut("solver", fileName, __FILE__, __LINE__);
}

//***********************************************************************

RelaxationPFinite::~RelaxationPFinite(){}

//***********************************************************************

void RelaxationPFinite::relaxation(Cell* cell, const double& dt, Prim type)
{
  if (dt < 1.e-20) return; //Check if dt = 0 (for AMR)

  //Note that the relaxation coefficients are herein considered as identical for each phase_k--phase_j combination

  //--------------------------------------------------------------------------------
  //---------------------------------- Parameters ----------------------------------
  //--------------------------------------------------------------------------------
  double epsilon(1.e-15);                      //To avoid division by 0
  double thresholdAlpha(1.e-8);                //Threshold on volume fraction for option alpha = 0 and for the computation of rhoCSquareMix
  double criterionRelaxedPressures(1.e-5);     //Criterion for considering relative pressures as relaxed
  int iterMax(100);                            //Maximum number of iteration for the pressure relaxation

  //--------------------------------------------------------------------------------
  //---------------------------------- Integration ---------------------------------
  //--------------------------------------------------------------------------------

  //Is the pressure-relaxation procedure necessary?
  //-----------------------------------------------
  //If alpha = 0 is activated, a test is done to know if the relaxation procedure is necessary or not
  //Else, i.e. alpha = 0 is desactivated (alpha != 0), the relaxation procedure is always done (relax = true)
  bool relax(true);
  if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
    for (int k = 0; k < numberPhases; k++) {
      if (cell->getPhase(k, type)->getAlpha() > (1. - thresholdAlpha)) relax = false; //KS//TODO: Should be dynamic ; Shouldn't there be a smart way to choose this alpha threshold?
    }
  }

  if (relax) {
    bool pressuresRelaxed(true);

    //Save initial state
    //------------------
    Phase* phase(0);
    for (int k = 0; k < numberPhases; k++) {
      phase = cell->getPhase(k, type);
      TB->ak[k] = phase->getAlpha();
      TB->rhok[k] = phase->getDensity();
      TB->pk[k] = phase->getPressure();
      if (TB->ak[k] < thresholdAlpha) TB->alphaNull[k] = true;
      else TB->alphaNull[k] = false;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------- Euler solver ---------------------------------
    //--------------------------------------------------------------------------------
    //The relaxation is done with an Euler, variable-step scheme
    if (!m_LSODA) {

      bool alphaSmall(false);
      double dtlocal(dt), tlocal(0.), dtmax(0.), dtnew(0.), pI(0.), pMixture(0.), sourceAlpha(0.), sourcePressure(0.), alphak0(0.);
      int iter = 0;
      while ((dt - tlocal) > 1.e-15) {
        iter++;

        //Initial local time-step
        //-----------------------
        dtlocal = dt - tlocal;

        //Mixture pressure
        //----------------
        pMixture = 0.;
        for (int k = 0; k < numberPhases; k++) {
          pMixture += cell->getPhase(k, type)->getAlpha() * cell->getPhase(k, type)->getPressure();
        }

        //Determine if pressures are already relaxed
        //------------------------------------------
        pressuresRelaxed = true;
        for (int k = 0; k < numberPhases; k++) {
          if (std::fabs((cell->getPhase(k, type)->getPressure() - pMixture) / cell->getPhase(k, type)->getPressure()) > criterionRelaxedPressures) pressuresRelaxed = false;
        }

        if (!pressuresRelaxed) {

          //Pressure differences
          //--------------------
          computePressureDifferences(cell, type);

          //Interface pressure
          //------------------
          pI = computeInterfacePressure(cell, type);
          for (int k = 0; k < numberPhases; k++) {
            cell->getPhase(k, type)->getEos()->verifyAndModifyPressure(pI);
          }

          //Local relaxation-coefficient determination
          //------------------------------------------
          //Only a few minor tests
          // double deltaPmax(1.);
          // for (int k = 0; k < numberPhases; k++) {
          // if (std::fabs(TB->Deltapk[k]) > deltaPmax) deltaPmax = std::fabs(TB->Deltapk[k]);
          // }
          // mu = m_muFactor / sqrt(deltaPmax);
          // mu = m_muFactor / deltaPmax;
          // mu = cell->getPhase(0, type)->getAlpha() * cell->getPhase(1, type)->getAlpha() * m_muFactor / deltaPmax;
          // mu = (2.65e-2 * cell->getPhase(0, type)->getAlpha() + 5.55e-3 * cell->getPhase(1, type)->getAlpha()) / sqrt(deltaPmax);
          // mu = (3.74e-3 * cell->getPhase(0, type)->getAlpha() + 5.55e-3 * cell->getPhase(1, type)->getAlpha()) / sqrt(deltaPmax);
          // mu = (1.18 * cell->getPhase(0, type)->getAlpha() + 5.55 * cell->getPhase(1, type)->getAlpha()) / deltaPmax;
          // mu = cell->getPhase(0, type)->getAlpha() * cell->getPhase(1, type)->getAlpha() * m_muFactor / sqrt(deltaPmax);

          //Local time-step determination
          //-----------------------------
          if (mu > 1.e-20) {
            dtmax = 1.e10;
            //Find local time-step based on volume-fraction restrictions
            for (int k = 0; k < numberPhases; k++) {
              sourceAlpha = mu * TB->Deltapk[k];
              if      (sourceAlpha >  epsilon) dtmax = std::fmin(dtmax, (1. - cell->getPhase(k, type)->getAlpha()) / sourceAlpha);
              else if (sourceAlpha < -epsilon) dtmax = std::fmin(dtmax,     - cell->getPhase(k, type)->getAlpha()  / sourceAlpha);
            }
            //Find local time-step based on pressure restrictions
            for (int k = 0; k < numberPhases; k++) {
              phase = cell->getPhase(k, type);
              TB->rho_cIksquare[k] = phase->getEos()->computeDensityTimesInterfaceSoundSpeedSquare(phase->getDensity(), pI, phase->getPressure());
            }
            for (int k = 0; k < numberPhases; k++) {
              phase = cell->getPhase(k, type);
              if (phase->getAlpha() > thresholdAlpha) {
                sourcePressure = - TB->rho_cIksquare[k] * TB->Deltapk[k] / phase->getAlpha();
                for (int j = 0; j < numberPhases; j++) {
                  sourcePressure -= (cell->getPhase(j, type)->getPressure() - TB->rho_cIksquare[j]) * TB->Deltapk[j];
                }
                dtnew = (pMixture - phase->getPressure()) / mu / sourcePressure;
                if (dtnew > 0.) dtmax = std::fmin(dtmax, dtnew);
              }
            }
            dtmax *= 0.5;
            if (dtmax < dtlocal) dtlocal = dtmax;
          }

          //Local relaxation and complete necessary variables for the next local time step
          //------------------------------------------------------------------------------
          //Relaxation on phase volume-fraction and pressure equations
          for (int k = 0; k < numberPhases; k++) {
            if (!TB->alphaNull[k]) {
              phase = cell->getPhase(k, type);
              alphak0 = phase->getAlpha();
              phase->setEnergy(phase->getEnergy() - dtlocal * mu * pI * TB->Deltapk[k] / std::max(alphak0 * phase->getDensity(), epsilon));
              phase->setAlpha(alphak0 + dtlocal * mu * TB->Deltapk[k]);
              phase->setDensity(alphak0 * phase->getDensity() / std::max(phase->getAlpha(), epsilon));
              phase->setPressure(phase->getEos()->computePressure(phase->getDensity(), phase->getEnergy()));
              phase->verifyAndCorrectPhase();
              phase->extendedCalculusPhase(cell->getMixture(type)->getVelocity());
            }
          }
        }
        else {
          //Break the while loop if pressures considered as relaxed
          break;
        }

        //Compute local time
        //------------------
        tlocal = tlocal + dtlocal;

        //Stop the relaxation procedure if a volume fraction is almost equal to 1 or if the number of iterations reaches a threshold
        //--------------------------------------------------------------------------------------------------------------------------
        if (epsilonAlphaNull > 1.e-20) { // alpha = 0 is activated
          for (int k = 0; k < numberPhases; k++) {
            if (cell->getPhase(k, type)->getAlpha() > (1. - thresholdAlpha)) alphaSmall = true;
          }
          if (alphaSmall) break;
        }

        if (iter == iterMax) {
          pressuresRelaxed = true;
          break;
        }
      }

      //If pressures considered as relaxed, we do an infinite relaxation to guarantee a unique pressure and better estimate the solution
      //--------------------------------------------------------------------------------------------------------------------------------
      if (pressuresRelaxed && mu > 1.e-20) {
        //Initial state
        double pStar(0.);
        for (int k = 0; k < numberPhases; k++) {
          pStar += TB->ak[k] * TB->pk[k];
        }

        //Iterative process for relaxed pressure determination
        int iteration(0);
        NewtonRaphson(pStar, iteration);

        //Apply the relaxation procedure if it has converged to a solution.
        if (iteration < 100) {
          //Cell update
          for (int k = 0; k < numberPhases; k++) {
            phase = cell->getPhase(k, type);
            phase->setAlpha(TB->akS[k]);
            phase->setDensity(TB->rhokS[k]);
            phase->setPressure(pStar);
            phase->extendedCalculusPhase(cell->getMixture(type)->getVelocity());
          }
        }
        //Else cell is updated with interface pressure obtained with values from previous convergence with Euler method
        else {
          pI = computeInterfacePressure(cell, type);
          for (int k = 0; k < numberPhases; k++) { TB->eos[k]->verifyAndModifyPressure(pI); } //Physical pressure?

          for (int k = 0; k < numberPhases; k++) {
            cell->getPhase(k, type)->setPressure(pI);
            phase->extendedCalculusPhase(cell->getMixture(type)->getVelocity());
          }
        }
      }

      //Cell update
      //-----------
      cell->getMixture(type)->computeMixtureVariables(cell->getPhases(type));
    }

    //--------------------------------------------------------------------------------
    //--------------------------------- LSODA solver ---------------------------------
    //--------------------------------------------------------------------------------
    else {

      //Determine if pressures are already relaxed
      //------------------------------------------
      pressuresRelaxed = true;
      for (int k = 0; k < numberPhases; k++) {
        if (std::fabs((cell->getPhase(k, type)->getPressure() - cell->getMixture(type)->getPressure()) / cell->getPhase(k, type)->getPressure()) > criterionRelaxedPressures) pressuresRelaxed = false;
      }

      if (!pressuresRelaxed) {

        //Pressure differences
        //--------------------
        computePressureDifferences(cell, type);

        //LSODA relaxation
        //----------------
        double t0(0.);
        int neq(2 * numberPhases);
        std::vector<double> y(neq + 1);
        std::vector<double>::iterator it(y.begin());
        for (int k = 0; k < numberPhases; k++) {
          *(++it) = TB->ak[k];
          *(++it) = TB->pk[k];
        }
        int istate = 1;
        LSODA lsoda;
        lsoda.lsoda_update(system_relaxation, neq, y, &t0, dt, &istate, this);

        //Cell update
        //-----------
        it = y.begin();
        for (int k = 0; k < numberPhases; k++) {
          phase = cell->getPhase(k, type);
          phase->setAlpha(*(++it));
          phase->setDensity(TB->ak[k] * TB->rhok[k] / std::max(phase->getAlpha(), epsilon));
          phase->setPressure(*(++it));
          phase->verifyAndCorrectPhase();
          phase->extendedCalculusPhase(cell->getMixture(type)->getVelocity());
        }

        if (istate <= 0) //Number of iterations to converge is too consequent, which happens when pressures are considered as relaxed
        {
          //If pressures considered as relaxed, we do an infinite relaxation to guarantee a unique pressure and better estimate the solution
          //--------------------------------------------------------------------------------------------------------------------------------
          //Initial state
          double pStar(0.);
          for (int k = 0; k < numberPhases; k++) {
            pStar += TB->ak[k] * TB->pk[k];
          }

          //Iterative process for relaxed pressure determination
          int iteration(0);
          NewtonRaphson(pStar, iteration);

          //Apply the relaxation procedure only if it has converged to a solution.
          if (iteration < 100) {
            //Cell update
            for (int k = 0; k < numberPhases; k++) {
              phase = cell->getPhase(k, type);
              phase->setAlpha(TB->akS[k]);
              phase->setDensity(TB->rhokS[k]);
              phase->setPressure(pStar);
              phase->extendedCalculusPhase(cell->getMixture(type)->getVelocity());
            }
          }
          //Else cell is updated with interface pressure obtained with values from previous convergence with LSODA method
          else {
            double pI(0.);
            pI = computeInterfacePressure(cell, type);
            for (int k = 0; k < numberPhases; k++) { TB->eos[k]->verifyAndModifyPressure(pI); } //Physical pressure?

            for (int k = 0; k < numberPhases; k++) {
              cell->getPhase(k, type)->setPressure(pI);
              phase->extendedCalculusPhase(cell->getMixture(type)->getVelocity());
            }
          }
        }

        //Cell update
        //-----------
        cell->getMixture(type)->computeMixtureVariables(cell->getPhases(type));
      }
    }
  }
}

//***********************************************************************

void RelaxationPFinite::computePressureDifferences(Cell* cell, Prim type)
{
  //Deltapk = alpha_k sum_j alpha_j (p_k - p_j) ; where N is the number of phases, j is different from k and where p_k also considers solid terms
  for (int k = 0; k < numberPhases; k++) {
    TB->Deltapk[k] = 0.;
    if (!TB->alphaNull[k]) {
      for (int j = 0; j < numberPhases; j++) {
        if (!TB->alphaNull[j]) {
          if (j != k) TB->Deltapk[k] += cell->getPhase(j, type)->getAlpha() * (cell->getPhase(k, type)->getPressure() - TB->compactionPk[k]
                                                                             - cell->getPhase(j, type)->getPressure() + TB->compactionPk[j]);
        }
      }
      TB->Deltapk[k] *= cell->getPhase(k, type)->getAlpha();
    }
  }
}

//***********************************************************************

void RelaxationPFinite::system_relaxation(double /*t*/, double* y, double* ydot, void* /*data*/)
{
  double epsilon(1.e-15), pI(0.), sumZk(0.), sumZj(0.);

  //Initial state
  //-------------
  Phase* phase(0);
  int i(0);
  for (int k = 0; k < numberPhases; k++) {
    phase = bufferCellLeft->getPhase(k);
    phase->setAlpha(y[i++]);
    phase->setDensity(TB->ak[k] * TB->rhok[k] / std::max(phase->getAlpha(), epsilon));
    phase->setPressure(y[i++]);
    phase->verifyAndCorrectPhase();
  }

  //Pressure differences
  //--------------------
  //Deltapk = alpha_k sum_j alpha_j (p_k - p_j) ; where N is the number of phases and j is different from k
  for (int k = 0; k < numberPhases; k++) {
    TB->Deltapk[k] = 0.;
    if (!TB->alphaNull[k]) {
      for (int j = 0; j < numberPhases; j++) {
        if (!TB->alphaNull[j]) {
          if (j != k) TB->Deltapk[k] += bufferCellLeft->getPhase(j)->getAlpha() * (bufferCellLeft->getPhase(k)->getPressure() - bufferCellLeft->getPhase(j)->getPressure());
        }
      }
      TB->Deltapk[k] *= bufferCellLeft->getPhase(k)->getAlpha();
    }
  }

  //Interface pressure
  //------------------
  //pI = sum_k (p_k sum_j Z_j) / (N - 1) sum_k (Z_k) ; where j is different from k
  for (int k = 0; k < numberPhases; k++) {
    phase = bufferCellLeft->getPhase(k);
    TB->zk[k] = phase->getEos()->computeAcousticImpedance(phase->getDensity(), phase->getPressure());
  }
  pI = 0.;
  sumZk = 0.;
  for (int k = 0; k < numberPhases; k++) {
    sumZk += TB->zk[k];
    sumZj = 0.;
    for (int j = 0; j < numberPhases; j++) {
      if (j != k) sumZj += TB->zk[j];
    }
    pI += sumZj * bufferCellLeft->getPhase(k)->getPressure();
  }
  pI /= sumZk * (numberPhases - 1);

  //System of equations
  //-------------------
  i = 0;
  for (int k = 0; k < numberPhases; k++) {
    phase = bufferCellLeft->getPhase(k);
    ydot[i++] = mu * TB->Deltapk[k];
    ydot[i++] = - phase->getEos()->computeDensityTimesInterfaceSoundSpeedSquare(phase->getDensity(), pI, phase->getPressure()) / std::max(phase->getAlpha(), epsilon) * mu * TB->Deltapk[k];
  }
}