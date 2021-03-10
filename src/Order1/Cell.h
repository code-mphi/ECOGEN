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

#ifndef CELL_H
#define CELL_H

#include <vector>
#include <fstream>
#include "../Models/Phase.h"
#include "../Maths/Coord.h"
#include "../Maths/Tensor.h"
#include "../Transport/Transport.h"

enum Variable { transport, pressure, density, alpha, velocityMag, velocityU, velocityV, velocityW, temperature, QPA };

class Cell; //Predeclaration of class to include following .h

#include "../Models/Mixture.h"
#include "../AdditionalPhysics/QuantitiesAddPhys.h"
#include "CellInterface.h"
#include "../Models/Model.h"
#include "../Models/Flux.h"
#include "../Meshes/Element.h"
#include "../Geometries/GeometricalDomain.h"
#include "../Symmetries/Symmetry.h"

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

        //!  \brief    Add a cell interface to current cell
        //!  \param    cellInterface   pointer to added cell interface
        void addCellInterface(CellInterface* cellInterface);
        //!  \brief    Delete the cell interface of current cell
        //!  \param    cellInterface   pointer to deleted cell interface
        void deleteCellInterface(CellInterface* cellInterface);

        //!  \brief    Memory allocation of cell attributes
        //!  \param    numberPhases       number of phases
        //!  \param    numberTransports   number of additional transport equations 
        //!  \param    addPhys            vector of additional physics
        //!  \param    model             hydrodynamical model used for memory allocation of the cell
        virtual void allocate(const int& numberPhases, const int& numberTransports, const std::vector<AddPhys*>& addPhys, Model* model);
        void allocateEos(const int& numberPhases, Model* model);
        //!  \brief    Filling cell properties using a physical domain
        //!  \param    domains           domain used for filling
        void fill(std::vector<GeometricalDomain*>& domains, const int& /*lvlMax*/);
        virtual void copyPhase(const int& phaseNumber, Phase* phase);
        void copyMixture(Mixture* mixture);
        void setToZeroCons(const int& numberPhases, const int& numberTransports);
        void setToZeroConsGlobal(const int& numberPhases, const int& numberTransports);
        void setToZeroBufferFlux(const int& numberPhases);
        void timeEvolution(const double& dt, const int& numberPhases, const int& numberTransports, Symmetry* symmetry, Prim type = vecPhases);
        void timeEvolutionAddPhys(const double& dt, const int& numberPhases);
        void buildPrim(const int& numberPhases);
        void buildCons(const int& numberPhases);
        void correctionEnergy(const int& numberPhases);
        void sourceTermIntegration(const double& /*dt*/, const int& /*numberPhases*/) {};
        void printPhasesMixture(const int& numberPhases, const int& numberTransports, std::ofstream& fileStream) const;
        virtual void completeFulfillState(Prim type = vecPhases);
        virtual void fulfillState(Prim type = vecPhases);
        void localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal, const int& numberPhases, Prim /*type*/ = vecPhases);
        virtual void reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal, const int& numberPhases);
        virtual void copyInCell(Cell& /*cellSource*/) const { Errors::errorMessage("methode copie non dispo pour cell"); };
        void copyVec(Phase** vecPhases, Mixture* mixture, Transport* vecTransports);
        //void printCut1Dde2D(std::ofstream& fileStream, std::string variableConstanteCut, const double& valueCut, const double& dL);    /*!< Ecriture de la cut 1D de la simulation 2D  */
        //void printCut1Dde3D(std::ofstream& fileStream, std::string variableConstanteCut1, std::string variableConstanteCut2, const double& valueCut1, const double& valueCut2, const double& dL1, const double& dL2);                                            /*!< Ecriture de la cut 1D de la simulation 3D  */

        //For additional physics
        //----------------------
        void prepareAddPhys();
        double selectScalar(Variable nameVariable, int num = 0) const;
        void setScalar(Variable nameVariable, const double& value, int num = 0);
        Coord selectVector(Variable nameVector, int num = 0, int subscript = -1) const;
        void setVector(Variable nameVector, const Coord& value, int num = 0, int subscript = -1);
        
        Coord computeGradient(Variable nameVariable, int num = -1);
        void computeGradient(std::vector<Coord>& grads, std::vector<Variable>& nameVariables, std::vector<int>& numPhases);
       
        const QuantitiesAddPhys* getQPA(int& numQPA) const { return m_vecQuantitiesAddPhys[numQPA]; }; //!< Allow to recover an additional physical quantity

        const Coord& getGradTk(int& numPhase, int& numAddPhys) const;
        void setGradTk(int& numPhase, int& numAddPhys, double* buffer, int& counter);
        void addNonConsAddPhys(const int& numberPhases, AddPhys& addPhys, Symmetry* symmetry);

        void reinitializeColorFunction(const int& numTransport, const int& numPhase); //!< Re-initialize the color function (transport) with alpha
        
        //Accessors
        //---------
        int getCellInterfacesSize() const;
        CellInterface* getCellInterface(const int& b);
        virtual Phase* getPhase(const int& phaseNumber, Prim /*type*/ = vecPhases) const;
        virtual Phase** getPhases(Prim /*type*/ = vecPhases) const;
        virtual Mixture* getMixture(Prim /*type*/ = vecPhases) const;
        Flux* getCons() const;
        void setCons(Flux* cons);
        const Coord& getPosition() const { return m_element->getPosition(); };
        const Coord& getSize() const { return m_element->getSize(); };
        const double& getSizeX() const { return m_element->getSizeX(); };
        const double& getSizeY() const { return m_element->getSizeY(); };
        const double& getSizeZ() const { return m_element->getSizeZ(); };
        void setElement(Element* element, const int& numCell);
        Element* getElement() const;
        virtual void setTransport(double value, int& numTransport, Prim /*type*/ = vecPhases);
        virtual Transport& getTransport(const int& numTransport, Prim /*type*/ = vecPhases) const;
        virtual Transport* getTransports(Prim /*type*/ = vecPhases) const;
        Transport* getConsTransport(const int& numTransport) const;
        void setConsTransport(double value, const int& numTransport);
        std::vector<QuantitiesAddPhys*>& getVecQuantitiesAddPhys();
        const int& getNumberPhases() const { return m_numberPhases; };
        const int& getNumberTransports() const { return m_numberTransports; };
        double getDensityGradient();
        Model* getModel();
        Coord& getVelocity();
        const Coord& getVelocity() const;

        //Not used for first order cells
        //------------------------------
        virtual void computeLocalSlopes(const int& /*numberPhases*/, const int& /*numberTransports*/, CellInterface& /*cellInterface*/,
            Limiter& /*globalLimiter*/, Limiter& /*interfaceLimiter*/, Limiter& /*globalVolumeFractionLimiter*/, Limiter& /*interfaceVolumeFractionLimiter*/,
            double& /*alphaCellAfterOppositeSide*/, double& /*alphaCell*/, double& /*alphaCellOtherInterfaceSide*/, double& /*epsInterface*/) {};           /*!< Does nothing for first order cells */
        virtual void computeLocalSlopesLimite(const int& /*numberPhases*/, const int& /*numberTransports*/, CellInterface& /*cellInterface*/,
            Limiter& /*globalLimiter*/, Limiter& /*interfaceLimiter*/, Limiter& /*globalVolumeFractionLimiter*/, Limiter& /*interfaceVolumeFractionLimiter*/,
            double& /*epsInterface*/) {};                                                                                                       /*!< Does nothing for first order cells */
        virtual Phase* getSlopes(const int& /*phaseNumber*/) const { return 0; };                                                               /*!< Does nothing for first order cells */
        virtual Transport* getSlopesTransport(const int& /*numberTransport*/) const { return 0; };                                              /*!< Does nothing for first order cells */
        virtual void saveCons(const int& /*numberPhases*/, const int& /*numberTransports*/) {};                                                     /*!< Does nothing for first order cells */
        virtual void recuperationCons(const int& /*numberPhases*/, const int& /*numberTransports*/) {};                                             /*!< Does nothing for first order cells */
        virtual void predictionOrdre2(const double& /*dt*/, const int& /*numberPhases*/, const int& /*numberTransports*/, Symmetry* /*symmetry*/) {};       /*!< Does nothing for first order cells */

        //CellInterface* getCellInterface(); //FP//TODO// A FAIRE...

        void printInfo() const;

        //methods for distance to an other object (Cell or CellInterface)
        //---------------------------------------------------------------
        double distance(Cell* c);            /*!< Distance totale  */
        double distanceX(Cell* c);           /*!< Distance selon x */
        double distanceY(Cell* c);           /*!< Distance selon y */
        double distanceZ(Cell* c);           /*!< Distance selon z */
        double distance(CellInterface* b);   /*!< Distance totale  */
        double distanceX(CellInterface* b);  /*!< Distance selon x */
        double distanceY(CellInterface* b);  /*!< Distance selon y */
        double distanceZ(CellInterface* b);  /*!< Distance selon z */

        bool traverseObjet(const GeometricObject& objet) const;

        //Printing
        //--------
        bool printGnuplotAMR(std::ofstream& fileStream, const int& dim, GeometricObject* objet = 0);
        void computeVolumePhaseK(double& integration, const int& numPhase);
        void computeMass(double& mass, double& alphaRef);
        void computeTotalMass(double& mass);
        void computeTotalEnergy(double& totalEnergy);
        void lookForPmax(double* pMax, double* pMaxWall);

        //Specific for AMR method
        //-----------------------
        void setToZeroXi();                                              /*!< set m_xi to zero */
        void setToZeroConsXi();                                          /*!< set m_consXi to zero */
        void timeEvolutionXi();                                          /*!< time evolution of Xi for smoothing */
        void chooseRefine(const double& xiSplit, const int& nbCellsY, const int& nbCellsZ,
          const std::vector<AddPhys*>& addPhys, Model* model, int& nbCellsTotalAMR); /*!< Choice for refinement of parent cell */
        void chooseUnrefine(const double& xiJoin, int& nbCellsTotalAMR); /*!< Choice for unrefinement of parent cell */
        void refineCellAndCellInterfaces(const int& nbCellsY, const int& nbCellsZ, const std::vector<AddPhys*>& addPhys, Model* model, const bool& refineExternalCellInterfaces);           /*!< Refine parent cell by creation of children cells */
        virtual void createChildCell(const int& lvl);                    /*!< Create a child cell (not initialized) */
        void unrefineCellAndCellInterfaces();                            /*!< Unrefine parent cell by destruction of children cells */
        void averageChildrenInParent();                                  /*!< Average variables of children cells in the parent cell, needed for computation of xi. */
        bool lvlNeighborTooHigh();                                       /*!< Look for AMR level of neighboring cells if too high to unrefine*/
        bool lvlNeighborTooLow();                                        /*!< Look for AMR level of neighboring cells if too low to refine*/
        void buildLvlCellsAndLvlInternalCellInterfacesArrays(std::vector<Cell*>* cellsLvl, std::vector<CellInterface*>* cellInterfacesLvl);      /*!< Build new arrays of cells and cell interfaces for level (lvl+1), only internal cell interfaces are added here */
        const int& getLvl() const { return m_lvl; };                     /*!< Get the cell AMR level in the AMR tree */
        const bool& getSplit() const { return m_split; };                /*!< Return true if the cells is plit, false otherwise */
        const double& getXi() const { return m_xi; };                    /*!< Return Xi cell value */
        void setXi(double value);                                        /*!< Set the Xi cell value */
        void addFluxXi(double value);                                    /*!< Add xi cell flux */
        int getNumberCellsChildren();                                    /*!< Return the number of children cells */
        Cell* getCellChild(const int& num);                              /*!< Return pointer to the corresponding child cell number */
        std::vector<Cell*>* getChildVector();                           /*!< Return pointer to child vector */

        //For parallel computing (no AMR)
        //-------------------------------
        virtual void pushBackSlope() {};                                 /*!< Does nothing for non-ghost O2 cells */
        virtual int getRankOfNeighborCPU() const { return -1; };
        virtual void setRankOfNeighborCPU(int /*rank*/) {};                  /*!< Does nothing for non-ghost cells */
        void fillBufferPrimitives(double* buffer, int& counter, const int& lvl, const int& neighbour, Prim type = vecPhases) const;
        void getBufferPrimitives(double* buffer, int& counter, const int& lvl, Eos** eos, Prim type = vecPhases);
        void fillBufferVector(double* buffer, int& counter, const int& lvl, const int& neighbour, const int& dim, Variable nameVector, int num = 0, int index = -1) const;
        void getBufferVector(double* buffer, int& counter, const int& lvl, const int& dim, Variable nameVector, int num = 0, int index = -1);
        void fillBufferTransports(double* buffer, int& counter, const int& lvl, const int& neighbour) const;
        void getBufferTransports(double* buffer, int& counter, const int& lvl);
        virtual void fillBufferSlopes(double* /*buffer*/, int& /*counter*/, const int& /*lvl*/, const int& /*neighbour*/) const {}; /*!< Does nothing for first order cells */
        virtual void getBufferSlopes(double* /*buffer*/, int& /*counter*/, const int& /*lvl*/) {};                              /*!< Does nothing for first order cells */
        virtual bool isCellGhost() const { return false; };
        bool hasNeighboringGhostCellOfCPUneighbour(const int& neighbour) const;                      /*!< Return a bool that is true if the cell has a neighboring ghost cell corresponding to CPU "neighbour" */
        int numberOfNeighboringGhostCellsOfCPUneighbour(const int& neighbour) const;                 /*!< Return the number of neighboring ghost cells corresponding to CPU "neighbour" this cell has */

        //For parallel AMR computing
        //--------------------------
        void chooseRefineDeraffineGhost(const int& nbCellsY, const int& nbCellsZ, const std::vector<AddPhys*>& addPhys, Model* model, std::vector<Cell*>* cellsLvlGhost); /*!< Choice for refinement, unrefinement of the ghost parent cell + Update of ghost cell vector for lvl+1 */
        void refineCellAndCellInterfacesGhost(const int& nbCellsY, const int& nbCellsZ, const std::vector<AddPhys*>& addPhys, Model* model);                               /*!< Refinement of parent ghost cell by creation of children ghost cells */
        void unrefineCellAndCellInterfacesGhost();                                                                                      /*!< Unrefinement of parent ghost cell by destruction of children ghost cells */
        void fillBufferXi(double* buffer, int& counter, const int& lvl, const int& neighbour) const;
        void getBufferXi(double* buffer, int& counter, const int& lvl);
        void fillBufferSplit(bool* buffer, int& counter, const int& lvl, const int& neighbour) const;
        void getBufferSplit(bool* buffer, int& counter, const int& lvl);
        void fillNumberElementsToSendToNeighbour(int& numberElementsToSendToNeighbor, int& numberSlopesToSendToNeighbor, const int& lvl, const int& neighbour, int numberNeighboursOfCPUneighbour);
        void fillDataToSend(std::vector<double>& dataToSend, std::vector<int>& dataSplitToSend, const int& lvl) const;
        void getDataToReceiveAndRefine(std::vector<double>& dataToReceive, std::vector<int>& dataSplitToReceive, const int& lvl, Eos** eos, int& counter, int& counterSplit,
            const int& nbCellsY, const int& nbCellsZ, const std::vector<AddPhys*>& addPhys, Model* model);
        void computeLoad(double& load, int lvl) const;
        void computeLvlMax(int& lvlMax) const;
        void clearExternalCellInterfaces(const int& nbCellsY, const int& nbCellsZ);
        void updatePointersInternalCellInterfaces();
        void updateNbCellsTotalAMR(int& nbCellsTotalAMR);

    protected:
      int m_numberPhases;
      int m_numberTransports;
      Phase** m_vecPhases;
      Mixture* m_mixture;
      Transport* m_vecTransports;                                 /*!< Array of passive variables advected in the flow */
      Flux* m_cons;                                               /*!< Conservative variables */
      Transport* m_consTransports;                                /*!< Array of fluxes for advected variables */
      Element *m_element;                                         /*!< Pointer to corresponding geometrical mesh element */
      std::vector<CellInterface*> m_cellInterfaces;               /*!< Vector of cell-interface pointers */
      std::vector<QuantitiesAddPhys*> m_vecQuantitiesAddPhys;     /*!< Vector of pointers to the Quantities of Additional Physics of the cell */
      Model* m_model;                                             /*!< Pointer to hydrodynamic model */
     
      //Attributs pour methode AMR
      int m_lvl;                                                  /*!< Cell AMR level in the AMR tree */
      double m_xi;                                                /*!< Criteria for refine/unrefine cell */
      double m_consXi;                                            /*!< Buffer variable for Xi fluxes */
	  bool m_split;                                               /*!< Indicator for splitted cell (Do I possess children ?) */
      std::vector<Cell*> m_childrenCells;                         /*!< Vector of children cells pointers */
      std::vector<CellInterface*> m_childrenInternalCellInterfaces; /*!< Vector of Internal children cell-interface pointers of the cell */

    private:
};

#endif // CELL_H
