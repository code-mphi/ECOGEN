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
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      December 20 2017

class CellInterface; //Predeclaration de la classe CellInterface pour pouvoir inclure Cell.h

#include "Cell.h"
#include "Models/Model.h"
#include "Models/Flux.h"
#include "Maths/Coord.h"
#include "Meshes/Face.h"
#include "Meshes/FaceCartesian.h"
#include "AdditionalPhysics/AddPhys.h"

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

		virtual int whoAmI() const {	return 0; };

    //Inutilise pour Bord de Maille ordre 1
    virtual void allocateSlopes(const int &numberPhases, const int &numberTransports, int &allocateSlopeLocal) {};   /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    virtual void computeSlopes(const int &numberPhases, const int &numberTransports, Prim type = vecPhases) {};  /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    virtual Phase* getSlopesPhase(const int &phaseNumber) const { return 0; };                                   /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    virtual Mixture* getSlopesMixture() const { return 0; };                                                     /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    virtual Transport* getSlopesTransport(const int &numberTransport) const { return 0; };                       /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual Cell *getB(BO2 B) const { return 0; };                                                          /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual double getBeta(betaO2 beta) const { return 0.; };                                                  /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual double getDistanceH(distanceHO2 dist) const { return 0.; };                                        /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual void setB(BO2 B, Cell *cell) {};                                                             /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual void setBeta(betaO2 beta, double &value) {};                                                      /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual void setDistanceH(distanceHO2 dist, double &value) {};                                            /*!< Ne fait rien pour des Bord de Maille ordre 1 */


    //Accesseurs
    Face *getFace();                                            /*!< Attention, getFace() non const */
    Model *getMod() const;
    Cell *getCellGauche() const;
    Cell *getCellDroite() const;
    virtual int getNumPhys() const { return -1; };
    //virtual double getDebit(int numPhase) const { Errors::errorMessage("getDebits non prevu pour CellInterface"); return 0.; }

    //Pour methode AMR
    virtual void computeXi(const double &criteriaVar, const bool &varRho, const bool &varP, const bool &varU, const bool &varAlpha);  /*!< Calcul de la variable Xi pour criteria de (de)raffinement a priori */
    void computeCritereAMR(const double &criteriaVar, std::string nameVariable, int num = 0);                                          /*!< Calcul de xi via le criteria de variation */
    virtual void computeFluxXi();                                 /*!< Calcul des flux de Xi (diffusion) pour smoothing */
    virtual void creerBordChild();                                                                        /*!< Creer un bord enfant (non initialize) */
    virtual void creerBordChildInterne(const int &lvl, std::vector<CellInterface*> *childrenInternalBoundaries); /*!< Creer un bord enfant interne (non initialize) */
    void creerFaceChild(CellInterface *bordParent);             /*!< Creer une face enfant (non initialize) */
    virtual void raffineBordExterne(const int &nbCellsY, const int &nbCellsZ, const double &dXParent, const double &dYParent, const double &dZParent, Cell *cellRef, const int &dim);      /*!< Raffinement du bord externe en creant si besoin des boundaries enfants + liaisons cells/boundaries */
		void raffineBordExterneGhost(const int &nbCellsY, const int &nbCellsZ, const double &dXParent, const double &dYParent, const double &dZParent, Cell *cellRef, const int &dim);      /*!< Raffinement du bord externe pour les cells fantomes en creant si besoin des boundaries enfants + liaisons cells/boundaries */
    virtual void deraffineBordExterne(Cell *cellRef);        /*!< Deraffinement du bord externe en supprimant si besoin ses boundaries enfants + liaisons cells/boundaries */
    void finalizeFace();                                        /*!< Supprime la face correspondante au bord */
    void deraffineBordsChildren();                               /*!< Supprime les boundaries enfants */
    void constructionTableauBordsExternesLvl(std::vector<CellInterface *> *boundariesLvl); /*!< Construction du nouveau tableau de boundaries du niveau (lvl + 1), boundaries externes ajoutes ici */
    bool getSplit() const;                                      /*!< Renvoie si oui ou non le bord est splitte */
    int getLvl() const;                                         /*!< Renvoie le niveau du bord */
    int getNumberBordsChildren() const;                          /*!< Renvoie le number de boundaries enfants de ce bord*/
    CellInterface *getBordChild(const int &numChild);          /*!< Renvoie le bord enfant correspondant au number */

   protected:
    Cell *m_cellLeft;
    Cell *m_cellRight;
    Model* m_mod;
    Face *m_face;
    
    //Attributs pour methode AMR
    int m_lvl;                                             /*!< Niveau dans l arbre AMR du bord */
    std::vector<CellInterface*> m_boundariesChildren;             /*!< Tableau de boundaries enfants (taille : 1 en 1D, 2 en 2D et 4 en 3D) */

  private:
};

//Utile pour la resolution des problemes de Riemann
extern Cell *cellLeft;
extern Cell *cellRight;

#endif // CELLINTERFACE_H
