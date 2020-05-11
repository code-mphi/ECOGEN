//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
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

#ifndef MIXEULER_H
#define MIXEULER_H

//! \file      MixEuler.h
//! \author    K. Schmidmayer, F. Petitpas
//! \version   1.1
//! \date      June 5 2019

#include <vector>
#include "../Mixture.h"

//! \class     MixEuler
//! \brief     Mixture variables for Euler equations (single phase)
class MixEuler : public Mixture
{
    public:
      MixEuler();
      virtual ~MixEuler();

      virtual void allocateAndCopyMixture(Mixture **mixture);
      virtual void copyMixture(Mixture &mixture) {};
      virtual double computeDensity(const double *alphak, const double *rhok, const int &numberPhases) { return 0.; };
      virtual double computePressure(const double *alphak, const double *pk, const int &numberPhases) { return 0.; };
      virtual double computeInternalEnergy(const double *Yk, const double *ek, const int &numberPhases) { return 0.; };
      virtual double computeFrozenSoundSpeed(const double *Yk, const double *ck, const int &numberPhases) { return 0.; };
      
      virtual void computeMixtureVariables(Phase **vecPhase, const int &numberPhases) {};
      virtual void internalEnergyToTotalEnergy(std::vector<QuantitiesAddPhys*> &vecGPA) {};
      virtual void totalEnergyToInternalEnergy(std::vector<QuantitiesAddPhys*> &vecGPA) {};

      virtual void localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal) {};
      virtual void reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal) {};

      //Specific methods for data printing
      //----------------------------------
      virtual int getNumberScalars() const { return 0; };
      virtual int getNumberVectors() const { return 0; };
      virtual double returnScalar(const int &numVar) const { return 0.; };
      virtual Coord returnVector(const int &numVar) const { return 0; };
      virtual std::string returnNameScalar(const int &numVar) const { return 0; };
      virtual std::string returnNameVector(const int &numVar) const { return 0; };

      //Specific methods for parallel computing
      //---------------------------------------
      virtual int numberOfTransmittedVariables() const { return 0; };
      virtual void fillBuffer(double *buffer, int &counter) const {};
      virtual void fillBuffer(std::vector<double> &dataToSend) const {};
      virtual void getBuffer(double *buffer, int &counter) {};
      virtual void getBuffer(std::vector<double> &dataToReceive, int &counter) {};

      //Specific methods for second order
      //---------------------------------
      virtual void computeSlopesMixture(const Mixture &sLeft, const Mixture &sRight, const double &distance) {};
      virtual void setToZero() {};
      virtual void extrapolate(const Mixture &slope, const double &distance) {};
      virtual void limitSlopes(const Mixture &slopeGauche, const Mixture &slopeDroite, Limiter &globalLimiter) {};

      //Specific methods for parallele computing at second order
      //--------------------------------------------------------
      virtual int numberOfTransmittedSlopes() const { return 0; };
      virtual void fillBufferSlopes(double *buffer, int &counter) const {};
      virtual void getBufferSlopes(double *buffer, int &counter) {};

      //Accessors
      //---------
      virtual const double& getDensity() const { return Errors::defaultDouble; };
      virtual const double& getPressure() const { return Errors::defaultDouble; };
      virtual const double& getU() const { return Errors::defaultDouble; };
      virtual const double& getV() const { return Errors::defaultDouble; };
      virtual const double& getW() const { return Errors::defaultDouble; };
      virtual const Coord& getVelocity() const { return Coord::defaultCoord; };
      virtual Coord& getVelocity() { return Coord::defaultCoordNonConst; };
      virtual const double& getEnergy() const { return Errors::defaultDouble; };
      virtual const double& getTotalEnergy() const { return Errors::defaultDouble; };
      virtual const double& getFrozenSoundSpeed() const { return Errors::defaultDouble; };
      virtual const double& getWoodSoundSpeed() const { return Errors::defaultDouble; };

      virtual void setPressure(const double &p) {};
      virtual void setVelocity(const double &u, const double &v, const double &w) {};
      virtual void setVelocity(const Coord &vit) {};
      virtual void setU(const double &u) {};
      virtual void setV(const double &v) {};
      virtual void setW(const double &w) {};
      virtual void setTotalEnergy(double &totalEnergy) {};

      //Operators
      //---------
      virtual void changeSign() {};
      virtual void multiplyAndAdd(const Mixture &slopesMixtureTemp, const double &coeff) {};
      virtual void divide(const double &coeff) {};

    protected:
    private:

};

#endif // MIXEULER_H
