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

#ifndef MODNONLINEARSCHRODINGER_H
#define MODNONLINEARSCHRODINGER_H

#include "../Model.h"
#include "../EulerKorteweg/ModEulerKorteweg.h"
#include "../../Order1/Cell.h"
#include "FluxNonLinearSchrodinger.h"
#include "MixNonLinearSchrodinger.h"

//! \class     ModNonLinearSchrodinger
//! \brief     Model class for Non-Linear Schrodinger mathematical system of equations (single phase)
class ModNonLinearSchrodinger : public ModEulerKorteweg
{
  public:
    //! \brief     NonLinearSchrodinger model constructor
    //! \param     numbTransports    number of additional transport equations
    //! \param     alpha               parameter alpha of Euler-Korteweg equations
    //! \param     beta                parameter beta of Euler-Korteweg equations
    ModNonLinearSchrodinger(const int& numbTransports, const double& alpha, const double &beta);
    virtual ~ModNonLinearSchrodinger();

    virtual void allocateCons(Flux** cons);
    virtual void allocatePhase(Phase** phase);
    virtual void allocateMixture(Mixture** mixture);

    //Methods specific to Euler-Korteweg
    //----------------------------------
    virtual double kappa(const double& density) const;
    virtual double kappaPrime(const double& density) const;
    virtual double kappaSecond(const double& density) const;
    virtual double epsilonPrime(Cell& /*cell*/, const double& /*density*/) const;
    virtual double epsilonSecond(Cell& /*cell*/, const double& /*density*/) const;

  private:
    static const std::string NAME;
};

#endif // MODNONLINEARSCHRODINGER_H
