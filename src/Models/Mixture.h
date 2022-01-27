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

#ifndef MIXTURE_H
#define MIXTURE_H

#include <vector>

class Mixture;

#include "../AdditionalPhysics/QuantitiesAddPhys.h"

//! \class     Mixture
//! \brief     Abstract class for mixture variables
class Mixture
{
    public:
      Mixture();
      virtual ~Mixture();
      //! \brief     Print mixture variables in file stream
      //! \param     fileStream      file stream to write in
      void printMixture(std::ofstream& fileStream) const;

      //! \brief     Compute saturation temperature for a liq/vapor couple of fluid at given pressure
      //! \param     eosLiq             pointer to equation of state of liquid phase
      //! \param     eosVap             pointer to equation of state of vapor phase
      //! \param     pressure           pressure
      //! \param     dTsat              temperature derivative as function of pressure
      //! \return    saturation temperature
      //virtual double computeTsat(const Eos* eosLiq, const Eos* eosVap, const double& pressure, double* dTsat=0) { Errors::errorMessage("computeTsat not available for required mixture"); return 0.; };
      double computeTsat(const Eos* eosLiq, const Eos* eosVap, const double& pressure, double* dTsat = 0);

      //! \brief     Copy mixture attributes in mixture
      //! \param     mixture      destination mixture variable 
      virtual void allocateAndCopyMixture(Mixture** /*mixture*/) { Errors::errorMessage("allocateAndCopyMixture not available for required mixture"); };
      //! \brief     Copy mixture in mixture attributes
      //! \param     mixture      source mixture to copy 
      virtual void copyMixture(Mixture& /*mixture*/) { Errors::errorMessage("copyMixture not available for required mixture"); };
      //! \brief     Compute mixture density
      //! \param     alphak             phase volume fraction array
      //! \param     rhok               phase density array
      //! \return    mixture density
      virtual double computeDensity(const double* /*alphak*/, const double* /*rhok*/) { Errors::errorMessage("computeDensity not available for required mixture"); return 0.; };
      //! \brief     Compute mixture pressure
      //! \param     alphak             phase volume fraction array
      //! \param     pk                 phase pressure array
      //! \return    mixture pressure
      virtual double computePressure(const double* /*alphak*/, const double* /*pk*/) { Errors::errorMessage("computePressure not available for required mixture"); return 0.; };
      virtual double computePressure(double* /*masses*/, const double& /*mixInternalEnerg*/, Phase** /*phases*/) { Errors::errorMessage("computePressure not available for required mixture"); return 0.; };
      virtual double computePressure(double /*mass*/, const double& /*internalEnergy*/, Phase** /*phases*/, Mixture* /*mixture*/, const int& /*liq*/, const int& /*vap*/) { Errors::errorMessage("computePressure not available for required mixture"); return 0.; };
      virtual double computeTemperature(double* /*masses*/, const double& /*pressure*/, Phase** /*phases*/) { Errors::errorMessage("computeTemperature not available for required mixture"); return 0.; };
      //! \brief     Compute mixture specific internal energy
      //! \param     Yk                 phase mass fraction array
      //! \param     ek                 phase specific internal energy array
      //! \return    mixture specific internal energy
      virtual double computeInternalEnergy(const double* /*Yk*/, const double* /*ek*/) { Errors::errorMessage("computeInternalEnergy not available for required mixture"); return 0.; };
      //! \brief     Compute mixture frozen speed of sound
      //! \param     Yk                 phase mass fraction array
      //! \param     ck                 phase speed of sound array
      //! \return    mixture frozen speed of sound
      virtual double computeFrozenSoundSpeed(const double* /*Yk*/, const double* /*ck*/) { Errors::errorMessage("computeFrozenSoundSpeed not available for required mixture"); return 0.; };

      //! \brief     Compute temperature for a mixture evolving at thermal equilibrium along mixture isentropic path
      //! \param     Yk                 array of mass fractions
      //! \param     p0                 initial pressure
      //! \param     T0                 initial temperature
      //! \param     p                  final pressure
      //! \param     dTdp               derivative according to pressure
      //! \return    temperature after isentropic path
      virtual double computeTemperatureIsentrope(const double* /*Yk*/, const double& /*p0*/, const double& /*T0*/, const double& /*p*/, double* /*dTdp*/ = 0) { Errors::errorMessage("computeTemperatureIsentrope not available for required mixture"); return 0.; };
      //! \brief     Compute mixture enthalpy for a mixture evolving at thermal equilibrium along mixture isentropic path
      //! \param     Yk                 array of mass fractions
      //! \param     p0                 initial pressure
      //! \param     T0                 initial temperature
      //! \param     p                  final pressure
      //! \param     dhdp               derivative according to pressure
      //! \return    enthalpy after isentropic path
      virtual double computeEnthalpyIsentrope(const double* /*Yk*/, const double& /*p0*/, const double& /*T0*/, const double& /*p*/, double* /*dhdp*/ = 0) { Errors::errorMessage("computeEnthalpyIsentrope not available for required mixture"); return 0.; };
      //! \brief     Compute mixture specific volume for a mixture evolving at thermal equilibrium along mixture isentropic path
      //! \param     Yk                 array of mass fractions
      //! \param     p0                 initial pressure
      //! \param     T0                 initial temperature
      //! \param     p                  final pressure
      //! \param     dvdp               derivative according to pressure
      //! \return    specific volume after isentropic path
      virtual double computeVolumeIsentrope(const double* /*Yk*/, const double& /*p0*/, const double& /*T0*/, const double& /*p*/, double* /*dvdp*/ = 0) { Errors::errorMessage("computeVolumeIsentrope not available for required mixture"); return 0.; };

      //! \brief     Fills some mixture attributes from a phase array
      //! \param     vecPhase           phase array
      virtual void computeMixtureVariables(Phase** /*vecPhase*/) { Errors::errorMessage("computeMixtureVariables not available for required mixture"); };
      //! \brief     Compute mixture total specific energy from internal one taking account for energies associated to extra physics
      //! \param     vecGPA             vector of additional physics variables
      virtual void internalEnergyToTotalEnergy(std::vector<QuantitiesAddPhys*>& /*vecGPA*/) { Errors::errorMessage("internalEnergyToTotalEnergy not available for required mixture"); };
      //! \brief     Compute mixture internal specific energy from total one taking account for energies associated to extra physics
      //! \param     vecGPA             vector of additional physics variables
      virtual void totalEnergyToInternalEnergy(std::vector<QuantitiesAddPhys*>& /*vecGPA*/) { Errors::errorMessage("totalEnergyToInternalEnergy not available for required mixture"); };
      
      //! \brief     velocity vector projection in a local Cartesian coordinate system
      //! \param     normal            normal vector associated to the cell interface
      //! \param     tangent           tangent vector associated to the cell interface
      //! \param     binormal          binormal vector associated to the cell interface
      virtual void localProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) { Errors::errorMessage("localProjection not available for required mixture"); };
      //! \brief     velocity vector reverse projection in the absolute Cartesian coordinate system
      //! \param     normal            normal vector associated to the cell interface
      //! \param     tangent           tangent vector associated to the cell interface
      //! \param     binormal          binormal vector associated to the cell interface
      virtual void reverseProjection(const Coord& /*normal*/, const Coord& /*tangent*/, const Coord& /*binormal*/) { Errors::errorMessage("reverseProjection not available for required mixture"); };

      //Specific methods for data printing
      //----------------------------------
      virtual int getNumberScalars() const { Errors::errorMessage("getNumberScalars non implemente pour mixture utilise"); return 0; };
      virtual int getNumberVectors() const { Errors::errorMessage("getNumberVectors non implemente pour mixture utilise"); return 0; };
      virtual double returnScalar(const int& /*numVar*/) const { Errors::errorMessage("returnScalar non implemente pour mixture utilise"); return 0.; };
      virtual Coord returnVector(const int& /*numVar*/) const { Errors::errorMessage("returnVector non implemente pour mixture utilise"); return 0; };
      virtual std::string returnNameScalar(const int& /*numVar*/) const { Errors::errorMessage("returnNameScalar non implemente pour mixture utilise"); return 0; };
      virtual std::string returnNameVector(const int& /*numVar*/) const { Errors::errorMessage("returnNameVector non implemente pour mixture utilise"); return 0; };

      //Specific method for reading from file
      //-------------------------------------
      virtual void setScalar(const int& /*numVar*/, const double& /*value*/) { Errors::errorMessage("setScalar not available for requested mixture type"); };
      virtual void setVector(const int& /*numVar*/, const Coord& /*value*/) { Errors::errorMessage("setVector not available for requested mixture type"); };

      //Specific methods for parallel computing
      //---------------------------------------
      virtual int numberOfTransmittedVariables() const { Errors::errorMessage("numberOfTransmittedVariables non implemente pour mixture utilise"); return 0; };
      virtual void fillBuffer(double* /*buffer*/, int& /*counter*/) const { Errors::errorMessage("fillBuffer non implemente pour mixture utilise"); };
      virtual void fillBuffer(std::vector<double>& /*dataToSend*/) const { Errors::errorMessage("fillBuffer non implemente pour mixture utilise"); };
      virtual void getBuffer(double* /*buffer*/, int& /*counter*/) { Errors::errorMessage("getBuffer non implemente pour mixture utilise"); };
      virtual void getBuffer(std::vector<double>& /*dataToReceive*/, int& /*counter*/) { Errors::errorMessage("getBuffer non implemente pour mixture utilise"); };

      //Specific methods for second order
      //---------------------------------
      virtual void computeSlopesMixture(const Mixture& /*sLeft*/, const Mixture& /*sRight*/, const double& /*distance*/) { Errors::errorMessage("computeSlopesMixture non implemente pour mixture utilise"); };
      virtual void setToZero() { Errors::errorMessage("setToZero non implemente pour mixture utilise"); };
      virtual void extrapolate(const Mixture& /*slope*/, const double& /*distance*/) { Errors::errorMessage("extrapolate non implemente pour mixture utilise"); };
      virtual void limitSlopes(const Mixture& /*slopeGauche*/, const Mixture& /*slopeDroite*/, Limiter& /*globalLimiter*/) { Errors::errorMessage("limitSlopes non implemente pour mixture utilise"); };

      //Specific methods for parallele computing at second order
      //--------------------------------------------------------
      virtual int numberOfTransmittedSlopes() const { Errors::errorMessage("numberOfTransmittedSlopes non implemente pour mixture utilise"); return 0; };
      virtual void fillBufferSlopes(double* /*buffer*/, int& /*counter*/) const { Errors::errorMessage("fillBufferSlopes non implemente pour mixture utilise"); };
      virtual void getBufferSlopes(double* /*buffer*/, int& /*counter*/) { Errors::errorMessage("getBufferSlopes non implemente pour mixture utilise"); };

      //Accessors
      //---------
      virtual const double& getDensity() const { Errors::errorMessage("getDensity not available for required mixture"); return Errors::defaultDouble; };
      virtual const double& getPressure() const { Errors::errorMessage("getPressure not available for required mixture"); return Errors::defaultDouble; };
      virtual const double& getTemperature() const { Errors::errorMessage("getTemperature not available for required mixture"); return Errors::defaultDouble; };
      virtual const double& getU() const { Errors::errorMessage("getU not available for required mixture"); return Errors::defaultDouble; };
      virtual const double& getV() const { Errors::errorMessage("getV not available for required mixture"); return Errors::defaultDouble; };
      virtual const double& getW() const { Errors::errorMessage("getW not available for required mixture"); return Errors::defaultDouble; };
      virtual const Coord& getVelocity() const { Errors::errorMessage("getVelocity not available for required mixture"); return Coord::defaultCoord; };
      virtual Coord& getVelocity() { Errors::errorMessage("getVelocity not available for required mixture"); return Coord::defaultCoordNonConst; };
      virtual const double& getEnergy() const { Errors::errorMessage("getEnergy not available for required mixture"); return Errors::defaultDouble; };
      virtual const double& getTotalEnergy() const { Errors::errorMessage("getTotalEnergy not available for required mixture"); return Errors::defaultDouble; };
      virtual const double& getFrozenSoundSpeed() const { Errors::errorMessage("getFrozenSoundSpeed not available for required mixture"); return Errors::defaultDouble; };
      virtual const double& getWoodSoundSpeed() const { Errors::errorMessage("getWoodSoundSpeed not available for required mixture"); return Errors::defaultDouble; };
      virtual const double& getMixSoundSpeed() const { Errors::errorMessage("getMixSoundSpeed not available for required mixture"); return Errors::defaultDouble; };

      virtual void setPressure(const double& /*p*/) { Errors::errorMessage("setPressure non implemente pour mixture utilise"); };
      virtual void setTemperature(const double& /*T*/) { Errors::errorMessage("setTemperature non implemente pour mixture utilise"); }
      virtual void setVelocity(const double& /*u*/, const double& /*v*/, const double& /*w*/) { Errors::errorMessage("setVelocity non implemente pour mixture utilise"); };
      virtual void setVelocity(const Coord& /*vit*/) { Errors::errorMessage("setVelocity non implemente pour mixture utilise"); };
      virtual void setU(const double& /*u*/) { Errors::errorMessage("setU non implemente pour mixture utilise"); };
      virtual void setV(const double& /*v*/) { Errors::errorMessage("setV non implemente pour mixture utilise"); };
      virtual void setW(const double& /*w*/) { Errors::errorMessage("setW non implemente pour mixture utilise"); };
      virtual void setTotalEnergy(double& /*totalEnergy*/) { Errors::errorMessage("setTotalEnergy not available for required mixture"); };

      //Operators
      //---------
      virtual void changeSign() { Errors::errorMessage("changeSign non implemente pour mixture utilise"); };
      virtual void multiplyAndAdd(const Mixture& /*slopesMixtureTemp*/, const double& /*coeff*/) { Errors::errorMessage("multiplyAndAdd non implemente pour mixture utilise"); };
      virtual void divide(const double& /*coeff*/) { Errors::errorMessage("divide non implemente pour mixture utilise"); };

    protected:
    private:

};

extern int numberScalarsMixture;

#endif // MIXTURE_H
