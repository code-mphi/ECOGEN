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

#ifndef CELL_H
#define CELL_H

//! \file      Cell.h
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      July 30 2018

#include <vector>
#include <fstream>
#include "Models/Phase.h"
#include "Maths/Coord.h"
#include "Transport/Transport.h"

class Cell; //Predeclaration of class to include following .h

#include "Models/Mixture.h"
#include "AdditionalPhysics/QuantitiesAddPhys.h"
#include "CellInterface.h"
#include "Models/Model.h"
#include "Models/Flux.h"
#include "Meshes/Element.h"
#include "Geometries/GeometricalDomain.h"
#include "Symmetries/Symmetry.h"

//! \class     Cell
//! \brief     Base class for a mesh cell
class Cell
{
    public:
        //! \brief     Basic Cell constructor for a non AMR cell
        Cell();
        //! \brief     Cell constructor for an AMR cell
        //! \param     lvl    level of current AMR cell
        Cell(int lvl); //Pour AMR
        virtual ~Cell();

        //!  \brief    Add a boundary to current cell
        //!  \param    bord   pointer to added cell boundary
        void addBoundary(CellInterface *bord);
        //!  \brief    Delete the boundary of current cell
        //!  \param    bord   pointer to deleted cell boundary
        void deleteBoundary(CellInterface *bord);

        //!  \brief    Memory allocation of cell attributes
        //!  \param    numberPhases       number of phases
        //!  \param    numberTransports   number of additional transport equations 
        //!  \param    addPhys            vector of additional physics
        //!  \param    model             hydrodynamical model used for memory allocation of the cell
        virtual void allocate(const int &numberPhases, const int &numberTransports, const std::vector<AddPhys*> &addPhys, Model *model);
        void allocateEos(const int &numberPhases, Model *model);
        //!  \brief    Filling cell properties using a physical domain
        //!  \param    domains           domain used for filling
        void fill(std::vector<GeometricalDomain*> &domains, const int &lvlMax);
        virtual void allocateAndCopyPhase(const int &phaseNumber, Phase *phase);
        virtual void copyPhase(const int &phaseNumber, Phase *phase);
        void copyMixture(Mixture *mixture);
        void setToZeroCons(const int &numberPhases, const int &numberTransports);
        void setToZeroConsGlobal(const int &numberPhases, const int &numberTransports);
        void setToZeroBufferFlux(const int &numberPhases);
        void timeEvolution(const double &dt, const int &numberPhases, const int &numberTransports, Symmetry *symmetry, Prim type = vecPhases);
        void timeEvolutionAddPhys(const double &dt, const int &numberPhases, const int &numberTransports);
        void buildPrim(const int &numberPhases);
        void buildCons(const int &numberPhases);
        void correctionEnergy(const int &numberPhases);
        void sourceTermIntegration(const double &dt, const int &numberPhases) {};
        void printPhasesMixture(const int &numberPhases, const int &numberTransports, std::ofstream &fileStream) const;
        virtual void completeFulfillState(Prim type = vecPhases);
        virtual void fulfillState(Prim type = vecPhases);
        virtual void localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal, const int &numberPhases, Prim type = vecPhases);
        virtual void reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal, const int &numberPhases, Prim type = vecPhases);
        virtual void copyInCell(Cell &cellSource, Prim type=vecPhases) const { Errors::errorMessage("methode copie non dispo pour cell"); };
        void copyVec(Phase **vecPhases, Mixture *mixture, Transport *vecTransports);
        //void printCut1Dde2D(std::ofstream &fileStream, std::string variableConstanteCut, const double &valueCut, const double &dL);    /*!< Ecriture de la cut 1D de la simulation 2D  */
        //void printCut1Dde3D(std::ofstream &fileStream, std::string variableConstanteCut1, std::string variableConstanteCut2, const double &valueCut1, const double &valueCut2, const double &dL1, const double &dL2);                                            /*!< Ecriture de la cut 1D de la simulation 3D  */

        //For additional physics
        //----------------------
        void prepareAddPhys();
        double selectScalar(std::string nameVariable, int num=0) const;
        void setScalar(std::string nameVariable, const double &value, int num = 0, int subscript = -1);
        Coord selectVector(std::string nameVector, int num=0, int subscript =-1) const;
        void setVector(std::string nameVector, const Coord &value, int num=0, int subscript =-1);
        
        Coord computeGradient(std::string nameVariable, int num=-1);
       
        QuantitiesAddPhys* getQPA(int &numQPA) const; //!< Allow to recover an additional physical quantity

        Coord getGradTk(int &numPhase, int &numAddPhys) const;
        void setGradTk(int &numPhase, int &numAddPhys, double *buffer, int &counter);
        void addNonConsAddPhys(const int &numberPhases, AddPhys &addPhys, Symmetry *symmetry);

				void reinitializeColorFunction(const int &numTransport, const int &numPhase); //!< Re-initialize the color function (transport) with alpha
        
        //Accessors
        //---------
        int getBordsSize() const;
        CellInterface* getBord(int &b);
        virtual Phase* getPhase(const int &phaseNumber, Prim type=vecPhases) const;
        virtual Phase** getPhases(Prim type=vecPhases) const;
        virtual Mixture* getMixture(Prim type = vecPhases) const;
        Flux* getCons() const;
        void setCons(Flux *cons);
        Coord getPosition() const;
        Coord getSize() const;
        double getSizeX() const;
        double getSizeY() const;
        double getSizeZ() const;
        void setElement(Element *element, const int &numCell);
        Element* getElement();
        virtual void setTransport(double value, int &numTransport, Prim type = vecPhases);
        virtual Transport& getTransport(const int &numTransport, Prim type = vecPhases) const;
        virtual Transport* getTransports(Prim type = vecPhases) const;
        Transport* getConsTransport(const int &numTransport) const;
        void setConsTransport(double value, const int &numTransport);
        std::vector<QuantitiesAddPhys*>& getVecQuantitiesAddPhys();
        int getNumberPhases() const;
        int getNumberTransports() const;
        double getXi() const { return m_xi; };
        double getGradient();
        Model *getModel();
        Coord getVelocity();

        //Not used for first order cells
        //------------------------------
        virtual void computeLocalSlopes(const int &numberPhases, const int &numberTransports, CellInterface &bord, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, double &alphaCellAfterOppositeSide, double &alphaCell, double &alphaCellOtherInterfaceSide, double &epsInterface) {};  /*!< Do nothing for first order cells */
        virtual void computeLocalSlopesLimite(const int &numberPhases, const int &numberTransports, CellInterface &bord, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter) {};  /*!< Do nothing for first order cells */
        virtual void computeMultiSlope(const int &numberPhases, CellInterface *bord, Limiter *globalLimiter) {};                            /*!< Do nothing for first order cells */
        virtual Phase* getSlopes(const int &phaseNumber) const { return 0; };                                                               /*!< Do nothing for first order cells */
        virtual Transport* getSlopesTransport(const int &numberTransport) const { return 0; };                                              /*!< Do nothing for first order cells */
        virtual void saveCons(const int &numberPhases, const int &numberTransports) {};                                                     /*!< Do nothing for first order cells */
        virtual void recuperationCons(const int &numberPhases, const int &numberTransports) {};                                             /*!< Do nothing for first order cells */
        virtual void predictionOrdre2(const double &dt, const int &numberPhases, const int &numberTransports, Symmetry *symmetry) {};       /*!< Do nothing for first order cells */

        void printInfo() const;

        //methods for distance to an other object (Cell or CellBoundary)
        //--------------------------------------------------------------
        double distance(Cell *c);            /*!< Distance totale  */
        double distanceX(Cell *c);           /*!< Distance selon x */
        double distanceY(Cell *c);           /*!< Distance selon y */
        double distanceZ(Cell *c);           /*!< Distance selon z */
        double distance(CellInterface *b);   /*!< Distance totale  */
        double distanceX(CellInterface *b);  /*!< Distance selon x */
        double distanceY(CellInterface *b);  /*!< Distance selon y */
        double distanceZ(CellInterface *b);  /*!< Distance selon z */

        bool traverseObjet(const GeometricObject &objet) const;

        //Printing
        //--------
        bool printGnuplotAMR(std::ofstream &fileStream, const int &dim, GeometricObject *objet = 0);
        void computeIntegration(double &integration);
        void lookForPmax(double *pMax, double *pMaxWall);

        //Specific for AMR method
        //-----------------------
        void setToZeroXi();                                              /*!< set m_xi to zero */
        void setToZeroConsXi();                                          /*!< set m_consXi to zero */
        void timeEvolutionXi();                                          /*!< time evolution of Xi for smoothing */
        void chooseRefine(const double &xiSplit, const int &nbCellsY, const int &nbCellsZ,
          const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR); /*!< Choice for refinement of parent cell */
        void chooseUnrefine(const double &xiJoin, int &nbCellsTotalAMR); /*!< Choice for unrefinement of parent cell */
        void refineCellAndBoundaries(const int &nbCellsY, const int &nbCellsZ, const std::vector<AddPhys*> &addPhys, Model *model);           /*!< Refine parent cell by creation of children cells */
        virtual void createChildCell(const int &num, const int &lvl);    /*!< Create a child cell (not initialized) */
        void unrefineCellAndBoundaries();                                /*!< Unrefine parent cell by destruction of children cells */
        void averageChildrenInParent();                                  /*!< Average variables of children cells in the parent cell, needed for computation of xi. */
        bool lvlNeighborTooHigh();                                       /*!< Look for AMR level of neighboring cells if too high to unrefine*/
        bool lvlNeighborTooLow();                                        /*!< Look for AMR level of neighboring cells if too low to refine*/
        void buildLvlCellsAndLvlInternalBoundariesArrays(std::vector<Cell *> *cellsLvl, std::vector<CellInterface *> *boundariesLvl);      /*!< Build new arrays of cells and boundaries for level (lvl+1), only internal boundaries are added here */
        int getLvl();                                                    /*!< Get the cell AMR level in the AMR tree */
        bool getSplit();                                                 /*!< Return true if the cells is plit, false otherwise */
        double getXi();                                                  /*!< Return Xi cell value */
        void setXi(double value);                                        /*!< Set the Xi cell value */
        void addFluxXi(double value);                                    /*!< Add xi cell flux */
        void subtractFluxXi(double value);                               /*!< Substract xi cell flux */
        int getNumberCellsChildren();                                    /*!< Return the number of children cells */
        Cell* getCellChild(const int &num);                              /*!< Return pointer to the corresponding child cell number */
        std::vector<Cell *>* getChildVector();                           /*!< Return pointer to chil vector */

        //For parallel computing (no AMR)
        //-------------------------------
        void fillBufferPrimitives(double *buffer, int &counter, Prim type = vecPhases) const;
        void getBufferPrimitives(double *buffer, int &counter, Eos **eos, Prim type = vecPhases);
        virtual void fillBufferSlopes(double *buffer, int &counter, std::string whichCpuAmIForNeighbour) const {};  /*!< Do nothing for first order cells */
        virtual void getBufferSlopes(double *buffer, int &counter) {};                                              /*!< Do nothing for first order cells */
        void fillBufferVector(double *buffer, int &counter, const int &dim, std::string nameVector, int num = 0, int index = -1) const;
        void getBufferVector(double *buffer, int &counter, const int &dim, std::string nameVector, int num = 0, int index = -1);
        void fillBufferTransports(double *buffer, int &counter) const;
        void getBufferTransports(double *buffer, int &counter);
        virtual bool isCellO2Ghost() const { return false; };

        //For parallel AMR computing
        //--------------------------
        void fillBufferPrimitivesAMR(double *buffer, int &counter, const int &lvl, std::string whichCpuAmIForNeighbour, Prim type = vecPhases) const;
        void getBufferPrimitivesAMR(double *buffer, int &counter, const int &lvl, Eos **eos, Prim type = vecPhases);
        virtual void fillBufferSlopesAMR(double *buffer, int &counter, const int &lvl, std::string whichCpuAmIForNeighbour) const {};   /*!< Do nothing for first order cells */
        virtual void getBufferSlopesAMR(double *buffer, int &counter, const int &lvl) {};                                               /*!< Do nothing for first order cells */
        void fillBufferVectorAMR(double *buffer, int &counter, const int &lvl, std::string whichCpuAmIForNeighbour, const int &dim, std::string nameVector, int num = 0, int index = -1) const;
        void getBufferVectorAMR(double *buffer, int &counter, const int &lvl, const int &dim, std::string nameVector, int num = 0, int index = -1);
        void fillBufferTransportsAMR(double *buffer, int &counter, const int &lvl, std::string whichCpuAmIForNeighbour) const;
        void getBufferTransportsAMR(double *buffer, int &counter, const int &lvl);
        void chooseRefineDeraffineGhost(const int &nbCellsY, const int &nbCellsZ, const std::vector<AddPhys*> &addPhys, Model *model, std::vector<Cell *> *cellsLvlGhost);           /*!< Choice for refinement, unrefinement of the ghost parent cell + Update of ghost cell vector for lvl+1 */
        void refineCellAndBoundariesGhost(const int &nbCellsY, const int &nbCellsZ, const std::vector<AddPhys*> &addPhys, Model *model);                                               /*!< Refinement of parent ghost cell by creation of children ghost cells */
        void unrefineCellAndBoundariesGhost();                                                               /*!< Unrefinement of parent ghost cell by destruction of children ghost cells */
        void fillBufferXi(double *buffer, int &counter, const int &lvl, std::string whichCpuAmIForNeighbour) const;
        void getBufferXi(double *buffer, int &counter, const int &lvl);
        void fillBufferSplit(bool *buffer, int &counter, const int &lvl, std::string whichCpuAmIForNeighbour) const;
        void getBufferSplit(bool *buffer, int &counter, const int &lvl);
        void fillNumberElementsToSendToNeighbour(int &numberNumberElementsToSendToNeighbor, const int &lvl, std::string whichCpuAmIForNeighbour);

    protected:
      int m_numberPhases;
      int m_numberTransports;
      Phase **m_vecPhases;
      Mixture *m_mixture;
      Transport *m_vecTransports;                                 /*!< Array of passive variables advected in the flow */
      Flux *m_cons;                                               /*!< Conservative variables */
      Transport *m_consTransports;                                /*!< Array of fluxes for advected variables */
      Element *m_element;                                         /*!< Pointer to corresponding geometrical mesh element */
      std::vector<CellInterface*> m_boundaries;                   /*!< Vector of boundaries cell pointers */
      std::vector<QuantitiesAddPhys*> m_vecQuantitiesAddPhys;     /*!< Vector of Pointors to the Quantities of Additional Physics of the cell */
      Model *m_model;                                             /*!< Pointer to hydrodynamic model */
     
      //Attributs pour methode AMR
      int m_lvl;                                                  /*!< Cell AMR level in the AMR tree */
      double m_xi;                                                /*!< Criteria for refine/unrefine cell */
      double m_consXi;                                            /*!< Buffer variable for Xi fluxes */
	  bool m_split;                                               /*!< Indicator for splitted cell (Do I possess children ?) */
      std::vector<Cell*> m_childrenCells;                         /*!< Vector of children cells pointers */
      std::vector<CellInterface*> m_childrenInternalBoundaries;   /*!< Vector of Internal children boundaries pointers of the cell */

    private:
};

#endif // CELL_H
