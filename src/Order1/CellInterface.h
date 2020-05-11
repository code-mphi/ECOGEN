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

#ifndef CELLINTERFACE_H
#define CELLINTERFACE_H

//! \file      CellInterface.h
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

class CellInterface; //Predeclaration de la classe CellInterface pour pouvoir inclure Cell.h

#include "Cell.h"
#include "../Models/Model.h"
#include "../Models/Flux.h"
#include "../Maths/Coord.h"
#include "../Meshes/Face.h"
#include "../Meshes/FaceCartesian.h"
#include "../AdditionalPhysics/AddPhys.h"

enum BO2 { BG1M, BG2M, BG3M, BG1P, BG2P, BG3P, BD1M, BD2M, BD3M, BD1P, BD2P, BD3P };
enum betaO2 { betaG1M, betaG2M, betaG3M, betaG1P, betaG2P, betaG3P, betaD1M, betaD2M, betaD3M, betaD1P, betaD2P, betaD3P };
enum distanceHO2 { distanceHGM, distanceHGP, distanceHDM, distanceHDP };

class CellInterface
{
  public:
    /** Default constructor */
    CellInterface();
    CellInterface(int lvl); //Pour AMR
    /** Default destructor */
    virtual ~CellInterface();

    void setFace(Face *face);

    virtual void computeFlux(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type = vecPhases);
    virtual void computeFluxAddPhys(const int &numberPhases, AddPhys &addPhys);
    virtual void solveRiemann(const int &numberPhases, const int &numberTransports, double &ondeMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type = vecPhases);
    virtual void initialize(Cell *cellLeft, Cell *cellRight);
    void initializeGauche(Cell *cellLeft);
    virtual void initializeDroite(Cell *cellRight);
    virtual void addFlux(const int &numberPhases, const int &numberTransports, const double &coefAMR);
    void subtractFlux(const int &numberPhases, const int &numberTransports, const double &coefAMR);
    double distance(Cell *c);

    void EffetsSurface1D(const int &numberPhases);

    void associeModel(Model *mod);

    virtual int whoAmI() const { return 0; };

    //Inutilise pour cell interfaces ordre 1
    virtual void allocateSlopes(const int &numberPhases, const int &numberTransports, int &allocateSlopeLocal) {};   /*!< Ne fait rien pour des cell interfaces ordre 1 */
    virtual void computeSlopes(const int &numberPhases, const int &numberTransports, Prim type = vecPhases) {};  /*!< Ne fait rien pour des cell interfaces ordre 1 */
    virtual Phase* getSlopesPhase(const int &phaseNumber) const { return 0; };                                   /*!< Ne fait rien pour des cell interfaces ordre 1 */
    virtual Mixture* getSlopesMixture() const { return 0; };                                                     /*!< Ne fait rien pour des cell interfaces ordre 1 */
    virtual Transport* getSlopesTransport(const int &numberTransport) const { return 0; };                       /*!< Ne fait rien pour des cell interfaces ordre 1 */
    //virtual Cell *getB(BO2 B) const { return 0; };                                                          /*!< Ne fait rien pour des cell interfaces ordre 1 */
    //virtual double getBeta(betaO2 beta) const { return 0.; };                                                  /*!< Ne fait rien pour des cell interfaces ordre 1 */
    //virtual double getDistanceH(distanceHO2 dist) const { return 0.; };                                        /*!< Ne fait rien pour des cell interfaces ordre 1 */
    //virtual void setB(BO2 B, Cell *cell) {};                                                             /*!< Ne fait rien pour des cell interfaces ordre 1 */
    //virtual void setBeta(betaO2 beta, double &value) {};                                                      /*!< Ne fait rien pour des cell interfaces ordre 1 */
    //virtual void setDistanceH(distanceHO2 dist, double &value) {};                                            /*!< Ne fait rien pour des cell interfaces ordre 1 */


    //Accesseurs
    Face *getFace();                                            /*!< Attention, getFace() non const */
    Model *getMod() const;
    Cell *getCellGauche() const;
    Cell *getCellDroite() const;
    virtual const int& getNumPhys() const { return Errors::defaultIntNeg; };
    //virtual double getDebit(int numPhase) const { Errors::errorMessage("getDebits non prevu pour CellInterface"); return 0.; }

    //Pour methode AMR
    virtual void computeXi(const double &criteriaVar, const bool &varRho, const bool &varP, const bool &varU, const bool &varAlpha);  /*!< Calcul de la variable Xi pour criteria de (de)raffinement a priori */
    void computeCritereAMR(const double &criteriaVar, Variable nameVariable, int num = 0);                                          /*!< Calcul de xi via le criteria de variation */
    virtual void computeFluxXi();                                      /*!< Calcul des flux de Xi (diffusion) pour smoothing */
    virtual void creerCellInterfaceChild();                                                                        /*!< Creer un child cell interface (non initialize) */
    virtual void creerCellInterfaceChildInterne(const int &lvl, std::vector<CellInterface*> *childrenInternalCellInterfaces); /*!< Creer un intern child cell interface (non initialize) */
    void creerFaceChild(CellInterface *cellInterfaceParent);           /*!< Creer une face enfant (non initialize) */
    virtual void raffineCellInterfaceExterne(const int &nbCellsY, const int &nbCellsZ, const double &dXParent, const double &dYParent, const double &dZParent, Cell *cellRef, const int &dim);      /*!< Raffinement du extern cell interface en creant si besoin des children cell interfaces + liaisons cells/cell interfaces */
    virtual void deraffineCellInterfaceExterne(Cell *cellRef);         /*!< Deraffinement du extern cell interface en supprimant si besoin ses children cell interfaces + liaisons cells/cell interfaces */
    void deraffineCellInterfacesChildren();                            /*!< Supprime les children cell interfaces */
    void constructionTableauCellInterfacesExternesLvl(std::vector<CellInterface *> *cellInterfacesLvl); /*!< Construction du nouveau tableau de cell interfaces du niveau (lvl + 1), cell interfaces externes ajoutes ici */
    bool getSplit() const;                                             /*!< Renvoie si oui ou non le cell interface est splitte */
    const int& getLvl() const { return m_lvl; };                       /*!< Renvoie le niveau du cell interface */
    int getNumberCellInterfacesChildren() const;                       /*!< Renvoie le number de children cell interfaces de ce cell interface*/
    CellInterface *getCellInterfaceChild(const int &numChild);         /*!< Renvoie le child cell interface correspondant au number */
    CellInterface *getCellInterfaceChildBack();                        /*!< Renvoie le child cell interface correspondant au number */
    void updatePointersInternalCellInterfaces();

   protected:
    Cell *m_cellLeft;
    Cell *m_cellRight;
    Model* m_mod;
    Face *m_face;
    
    //Attributs pour methode AMR
    int m_lvl;                                             /*!< Niveau dans l arbre AMR du cell interface */
    std::vector<CellInterface*> m_cellInterfacesChildren;  /*!< Tableau de children cell interfaces (taille : 1 en 1D, 2 en 2D et 4 en 3D) */

  private:
};

//Utile pour la resolution des problemes de Riemann
extern Cell *cellLeft;
extern Cell *cellRight;

#endif // CELLINTERFACE_H
