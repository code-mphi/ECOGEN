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

#ifndef MODEULERKORTEWEG_H
#define MODEULERKORTEWEG_H

#include "../Model.h"
#include "../../Order1/Cell.h"
#include "FluxEulerKorteweg.h"
#include "MixEulerKorteweg.h"

//! \class     ModEulerKorteweg
//! \brief     Model class for Augmented Euler--Korteweg mathematical system of equations (single phase)
class ModEulerKorteweg : public Model
{
  public:
    //! \brief     EulerKorteweg model constructor
    //! \param     numbTransports    number of additional transport equations
    ModEulerKorteweg(const int& numbTransports);
    //! \brief     EulerKorteweg model constructor
    //! \param     numbTransports       number of additional transport equations
    //! \param     alpha                parameter alpha of Euler-Korteweg equations
    //! \param     beta                 parameter beta of Euler-Korteweg equations
    //! \param     temperature          constant temperature of Euler-Korteweg equations
    //! \param     kappa                constant surface-tension coefficient of Euler-Korteweg equations
    ModEulerKorteweg(const int& numbTransports, const double& alpha, const double &beta, const double &temperature, const double &kappa);
    //! \brief     EulerKorteweg model constructor
    //! \param     name                 model name
    //! \param     numbTransports       number of additional transport equations
    //! \param     alpha                parameter alpha of Euler-Korteweg equations
    //! \param     beta                 parameter beta of Euler-Korteweg equations
    //! \param     temperature          constant temperature of Euler-Korteweg equations
    //! \param     kappa                constant surface-tension coefficient of Euler-Korteweg equations
    ModEulerKorteweg(const std::string& name, const int& numbTransports, const double& alpha, const double &beta, const double &temperature, const double &kappa);
    virtual ~ModEulerKorteweg();

    virtual void allocateCons(Flux** cons);
    virtual void allocatePhase(Phase** phase);
    virtual void allocateMixture(Mixture** mixture);

    //! \details    Does nothing for this model
    virtual void fulfillState(Phase** /*phases*/, Mixture* /*mixture*/) {};

    //! \details    Does nothing for this model
    virtual void fulfillStateRestart(Phase** /*phases*/, Mixture* /*mixture*/) {};
    
    //! \brief     Complete the augmented variables (such as the ones of Euler-Korteweg model)
    //! \param     cell           cell
    virtual void initializeAugmentedVariables(Cell* cell);

    //Hydrodynamic Riemann solvers
    //----------------------------
    virtual void solveRiemannIntern(Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double>& /*boundData*/ = DEFAULT_VEC_INTERFACE_DATA) const;

    virtual void reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const;

    //Accessors
    //---------
    //! \brief  Select a specific scalar variable
    //! \param  phases         phases array variables
    //! \param  mixture        mixture variables
    //! \param  vecTransports  vector of transport variables
    //! \param  nameVariables  Name of the variable to select
    //! \param  numPhases      Phases number's
    virtual double selectScalar(Phase** phases, Mixture* /*mixture*/, Transport* transports, Variable nameVariable, int num = 0) const;
    virtual const double& getSM() { Errors::errorMessage("getSM not available for required model"); return Errors::defaultDouble; };
    virtual const Coord& getVelocity(const Cell* cell) const { return cell->getPhase(0)->getVelocity(); };
    virtual Coord& getVelocity(Cell* cell) { return cell->getPhase(0)->getVelocity(); };

    virtual const std::string& whoAmI() const { return m_name; };

    //Methods specific to Euler-Korteweg
    //----------------------------------
    virtual double kappa(const double& /*density*/) const;
    virtual double kappaPrime(const double& /*density*/) const;
    virtual double kappaSecond(const double& /*density*/) const;
    virtual double epsilonPrime(Cell& cell, const double& density) const;
    virtual double epsilonSecond(Cell& cell, const double& density) const;
    //! \brief     Compute maximal wave speed using Davis approximation and eigenvalues of hyperbolic equations
    //! \param     cellLeft     cell left
    //! \param     cellRight    cell right
    //! \param     rhoL         density left
    //! \param     rhoR         density right
    //! \param     uL           velocity left
    //! \param     uR           velocity right
    //! \param     etaL         eta left
    //! \param     etaR         eta right
    //! \param     vecPxL       x-component of vectorP left
    //! \param     vecPxR       x-component of vectorP right
    //! \param     vecPyL       y-component of vectorP left
    //! \param     vecPyR       y-component of vectorP right
    virtual double computeMaxWaveSpeed(Cell& cellLeft, Cell& cellRight, const double& rhoL, const double& rhoR, const double& uL, const double& uR,
                                       const double& etaL, const double& etaR, const double& vecPxL, const double& vecPxR,
                                       const double& vecPyL, const double& vecPyR) const;

  private:
    static const std::string NAME;
};

#endif // MODEULERKORTEWEG_H
