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

#ifndef CELLINTERFACEO2_H
#define CELLINTERFACEO2_H

//! \file      CellInterfaceO2.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "../Order1/CellInterface.h"

//class CellInterfaceO2; //Predeclaration de la classe CellInterfaceO2 pour pouvoir inclure CellO2.h
//
//#include "CellO2.h"

class CellInterfaceO2 : public CellInterface
{
  public:
    /** Default constructor */
    CellInterfaceO2();
    CellInterfaceO2(int lvl); //Pour AMR
    /** Default destructor */
    virtual ~CellInterfaceO2();

    virtual void allocateSlopes(const int &numberPhases, const int &numberTransports, int &allocateSlopeLocal);
    virtual void computeSlopes(const int &numberPhases, const int &numberTransports, Prim type = vecPhases);
    virtual void computeFlux(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type = vecPhases);
    void solveRiemann(const int &numberPhases, const int &numberTransports, double &ondeMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type = vecPhases); /*!< probleme de Riemann special ordre 2 */

    //Accesseurs
    virtual Phase* getSlopesPhase(const int &phaseNumber) const;
    virtual Mixture* getSlopesMixture() const;
    virtual Transport* getSlopesTransport(const int &numberTransport) const;
    //virtual Cell *getB(BO2 B) const;
    //virtual double getBeta(betaO2 beta) const;
    //virtual double getDistanceH(distanceHO2 dist) const;
    //virtual void setB(BO2 B, Cell *cell);
    //void setBeta(betaO2 beta, double &value);
    //virtual void setDistanceH(distanceHO2 dist, double &value);

    //Pour methode AMR
    virtual void creerCellInterfaceChild();                                                                        /*!< Creer un child cell interface (non initialize) */
    virtual void creerCellInterfaceChildInterne(const int &lvl, std::vector<CellInterface*> *childrenInternalCellInterfaces); /*!< Creer un intern child cell interface (non initialize) */

   protected:
     int m_numberPhases;
     Phase **m_vecPhasesSlopes;         /*!< vecteur des slopes des phases */
     Mixture *m_mixtureSlopes;          /*!< vecteur des slopes de mixture */
     Transport *m_vecTransportsSlopes;	/*!< vecteur des slopes des transports */

     //Stockage methode multislopes
     //Cell *m_BG1M; /*!< pointeurs vers cells Arrieres a gauche pour secondOrder  */
     //Cell *m_BG2M;
     //Cell *m_BG3M;
     //Cell *m_BG1P; /*!< pointeurs vers cells Avants a gauche pour secondOrder  */
     //Cell *m_BG2P;
     //Cell *m_BG3P;
     //Cell *m_BD1M; /*!< pointeurs vers cells Arrieres a droite pour secondOrder  */
     //Cell *m_BD2M;
     //Cell *m_BD3M;
     //Cell *m_BD1P; /*!< pointeurs vers cells Avants a droite pour secondOrder  */
     //Cell *m_BD2P;
     //Cell *m_BD3P;

     //double m_betaG1M;  /*!< ponderations pour secondOrder */
     //double m_betaG2M;
     //double m_betaG3M;
     //double m_betaG1P;
     //double m_betaG2P;
     //double m_betaG3P;
     //double m_betaD1M;
     //double m_betaD2M;
     //double m_betaD3M;
     //double m_betaD1P;
     //double m_betaD2P;
     //double m_betaD3P;

     //double m_distanceHGM;  /*!< distances au vertex geometrique pour le compute des slopes */
     //double m_distanceHGP;
     //double m_distanceHDM;
     //double m_distanceHDP;

   private:
};

extern Phase **slopesPhasesLocal1;
extern Phase **slopesPhasesLocal2;
extern Mixture *slopesMixtureLocal1;
extern Mixture *slopesMixtureLocal2;
extern double *slopesTransportLocal1;
extern double *slopesTransportLocal2;

#endif // CELLINTERFACEO2_H
