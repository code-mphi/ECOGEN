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

#include <fstream>
#include "../Models/Phase.h"
#include "../Maths/Coord.h"
#include "../Maths/Tensor.h"
#include "../Transport/Transport.h"
#include "../Transport/GradTransport.h"

class Cell; //Predeclaration of class to include following .h

#include "../Models/Mixture.h"
#include "../AdditionalPhysics/QuantitiesAddPhys.h"
#include "CellInterface.h"
#include "../Models/Model.h"
#include "../Gradients/Gradient.h"
#include "../Models/Flux.h"
#include "../Meshes/Element.h"
#include "../Geometries/GeometricalDomain.h"
#include "../Symmetries/Symmetry.h"

class GradPhase;
class GradMixture;

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

        //!  \brief    Associate external variables
        //!  \param    mod       pointer to model
        //!  \param    grad      pointer to gradient method
        void associateExtVar(Model* mod, Gradient* grad);
        //!  \brief    Memory allocation of cell attributes
        //!  \param    addPhys            vector of additional physics
        virtual void allocate(const std::vector<AddPhys*>& addPhys);
        void allocateEos();
        //!  \brief    Filling cell properties using a physical domain
        //!  \param    domains           domain used for filling
        void fill(std::vector<GeometricalDomain*>& domains, const int& /*lvlMax*/);
        virtual void copyPhase(const int& phaseNumber, Phase* phase);
        void copyMixture(Mixture* mixture);
        void setToZeroCons();
        void setToZeroConsGlobal();
        void timeEvolution(const double& dt, Symmetry* symmetry);
        void timeEvolutionAddPhys(const double& dt);
        void buildPrim();
        void buildCons();
        void correctionEnergy();
        void sourceTermIntegration(const double& /*dt*/) {};
        void printPhasesMixture(std::ofstream& fileStream) const;
        virtual void completeFulfillState();
        virtual void fulfillState(Prim /*type*/ = vecPhases);
        void fulfillStateRestart();
        void localProjection(const Coord& normal, const Coord& tangent, const Coord& binormal, Prim /*type*/ = vecPhases);
        virtual void reverseProjection(const Coord& normal, const Coord& tangent, const Coord& binormal);
        virtual void copyInCell(Cell& /*cellSource*/) const { Errors::errorMessage("methode copie non dispo pour cell"); };
        void copyVec(Phase** vecPhases, Mixture* mixture, Transport* vecTransports);
        //void printCut1Dde2D(std::ofstream& fileStream, std::string variableConstanteCut, const double& valueCut, const double& dL);    /*!< Write de la cut 1D de la simulation 2D  */
        //void printCut1Dde3D(std::ofstream& fileStream, std::string variableConstanteCut1, std::string variableConstanteCut2, const double& valueCut1, const double& valueCut2, const double& dL1, const double& dL2);                                            /*!< Write de la cut 1D de la simulation 3D  */
        void deleteInterface(const int &b);

        //For additional physics
        //----------------------
        void prepareAddPhys();
        Coord selectVector(Variable nameVector, int num = 0, int subscript = -1) const;
        void setVector(Variable nameVector, const Coord& value, int num = 0, int subscript = -1);

        const QuantitiesAddPhys* getQPA(int& numQPA) const { return m_vecQuantitiesAddPhys[numQPA]; }; //!< Allow to recover an additional physical quantity

        const Coord& getGradTk(int& numPhase, int& numAddPhys) const;
        void setGradTk(int& numPhase, int& numAddPhys, double* buffer, int& counter);
        void addNonConsAddPhys(AddPhys& addPhys, Symmetry* symmetry);

        void reinitializeColorFunction(const int& numTransport, const int& numPhase); //!< Re-initialize the color function (transport) with alpha

        //Gradients
        //---------
        //! \brief  Compute gradients (temperature, velocity, density) of a cell
        //! \param  grads          Array of desired gradients, e.g. for temperature each component represents phase temperature and for velocity each component represents the gradient of a velocity component (grad(u), grad(v), grad(w))
        //! \param  nameVariables  Name of the variable for which the gradient is calculated
        //! \param  numPhases      Phases number's
        void computeGradients(std::vector<Coord>& grads, std::vector<Variable>& nameVariables, std::vector<int>& numPhases);
        
        //Accessors
        //---------
        int getCellInterfacesSize() const;
        CellInterface* getCellInterface(const int& b);
        void setCellInterface(const int& b, CellInterface* cellInterface);
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
        const int& getNumberPhases() const;
        const int& getNumberTransports() const;
        double getDensityGradient();
        Model* getModel();
        Coord& getVelocity();
        const Coord& getVelocity() const;
        void setWall(bool wall);
        bool getWall() const { return m_wall; };
        //! \brief  Select a specific scalar variable
        //! \param  nameVariables  Name of the variable to select
        //! \param  numPhases      Phases number's
        double selectScalar(Variable nameVariable, int num = 0) const;

        //Second order (not used for first-order cells)
        //---------------------------------------------
        //! \brief  Compute global variable buffers (min, max, etc.) and initialize speficic gradient vectors for 2nd-order scheme on unstructured mesh
        virtual void allocateSecondOrderBuffersAndGradientVectors(Phase** /*phases*/, Mixture* /*mixture*/) {};

        //! \brief  Compute gradients for 2nd-order scheme on unstructured mesh
        virtual void computeGradientsO2() {};
        virtual void limitGradientsO2(Limiter& /*globalLimiter*/) {};

        //! \brief  Compute slopes for 2nd-order scheme on unstructured mesh
        virtual void computeLocalSlopes(CellInterface& /*cellInterfaceRef*/) {};
        virtual void computeLocalSlopes(CellInterface& /*cellInterface*/, Limiter& /*globalLimiter*/, Limiter& /*interfaceLimiter*/,
            Limiter& /*globalVolumeFractionLimiter*/, Limiter& /*interfaceVolumeFractionLimiter*/,
            double& /*alphaCellAfterOppositeSide*/, double& /*alphaCell*/, double& /*alphaCellOtherInterfaceSide*/,
            double& /*epsInterface*/) {};                                               /*!< Does nothing for first order cells */
        virtual void computeLocalSlopesLimite(CellInterface& /*cellInterface*/, Limiter& /*globalLimiter*/,
            Limiter& /*interfaceLimiter*/, Limiter& /*globalVolumeFractionLimiter*/, Limiter& /*interfaceVolumeFractionLimiter*/,
            double& /*epsInterface*/) {};                                               /*!< Does nothing for first order cells */
        virtual Phase* getSlopes() const { return 0; };                                 /*!< Does nothing for first order cells */
        virtual Transport* getSlopesTransport() const { return 0; };                    /*!< Does nothing for first order cells */
        virtual void saveCons() {};                                                     /*!< Does nothing for first order cells */
        virtual void getBackCons() {};                                             /*!< Does nothing for first order cells */
        virtual void predictionOrdre2(const double& /*dt*/, Symmetry* /*symmetry*/) {}; /*!< Does nothing for first order cells */

        //Methods for distance to an other object (Cell or CellInterface)
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
        void printInfo() const;
        //! \brief  Compute saturation pressure for a liquid/vapor fluid mixture (first phase is considered predominant -> not generalized be careful)
        double getPsat();
        bool printGnuplotAMR(std::ofstream& fileStream, const int& dim, GeometricObject* objet = 0, bool recordPsat = false);
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
          const std::vector<AddPhys*>& addPhys, int& nbCellsTotalAMR);   /*!< Choice for refinement of parent cell */
        void chooseUnrefine(const double& xiJoin, int& nbCellsTotalAMR); /*!< Choice for unrefinement of parent cell */
        void refineCellAndCellInterfaces(const int& nbCellsY, const int& nbCellsZ, const std::vector<AddPhys*>& addPhys, const bool& refineExternalCellInterfaces);  /*!< Refine parent cell by creation of children cells */
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
        std::vector<Cell*>* getChildVector();                            /*!< Return pointer to child vector */

        //For parallel computing (no AMR)
        //-------------------------------
        virtual void pushBackSlope() {};                                 /*!< Does nothing for non-ghost O2 cells */
        virtual void popBackSlope() {};                                  /*!< Does nothing for non-ghost O2 cells */
        virtual int getRankOfNeighborCPU() const { return -1; };
        virtual void setRankOfNeighborCPU(int /*rank*/) {};              /*!< Does nothing for non-ghost cells */
        void fillBufferPrimitives(double* buffer, int& counter, const int& lvl, const int& neighbour, Prim type = vecPhases) const;
        void getBufferPrimitives(double* buffer, int& counter, const int& lvl, Eos** eos, Prim type = vecPhases);
        void fillBufferVector(double* buffer, int& counter, const int& lvl, const int& neighbour, const int& dim, Variable nameVector, int num = 0, int index = -1) const;
        void getBufferVector(double* buffer, int& counter, const int& lvl, const int& dim, Variable nameVector, int num = 0, int index = -1);
        void fillBufferTransports(double* buffer, int& counter, const int& lvl, const int& neighbour) const;
        void getBufferTransports(double* buffer, int& counter, const int& lvl);
        virtual void fillBufferSlopes(double* /*buffer*/, int& /*counter*/, const int& /*lvl*/, const int& /*neighbour*/) const { Errors::errorMessage("fillBufferSlopes not available for Cell"); }; /*!< Does nothing for first order cells */
        virtual void getBufferSlopes(double* /*buffer*/, int& /*counter*/, const int& /*lvl*/) { Errors::errorMessage("getBufferSlopes not available for Cell"); };                                  /*!< Does nothing for first order cells */
        virtual bool isCellGhost() const { return false; };
        bool hasNeighboringGhostCellOfCPUneighbour(const int& neighbour) const;                      /*!< Return a bool that is true if the cell has a neighboring ghost cell corresponding to CPU "neighbour" */
        int numberOfNeighboringGhostCellsOfCPUneighbour(const int& neighbour) const;                 /*!< Return the number of neighboring ghost cells corresponding to CPU "neighbour" this cell has */
        virtual GradPhase* getGradPhase(const int& /*phaseNumber*/) const { Errors::errorMessage("getGradPhase not available for Cell"); return nullptr; };
        virtual GradMixture* getGradMixture() const { Errors::errorMessage("getGradMixture not available for Cell"); return nullptr; };
        virtual GradTransport* getGradTransport(const int& /*transportNumber*/) const { Errors::errorMessage("getGradTransport not available for Cell"); return nullptr; };

        //For parallel AMR computing
        //--------------------------
        void chooseRefineDeraffineGhost(const int& nbCellsY, const int& nbCellsZ, const std::vector<AddPhys*>& addPhys, std::vector<Cell*>* cellsLvlGhost); /*!< Choice for refinement, unrefinement of the ghost parent cell + Update of ghost cell vector for lvl+1 */
        void refineCellAndCellInterfacesGhost(const int& nbCellsY, const int& nbCellsZ, const std::vector<AddPhys*>& addPhys);                              /*!< Refinement of parent ghost cell by creation of children ghost cells */
        void unrefineCellAndCellInterfacesGhost();                                                                                                          /*!< Unrefinement of parent ghost cell by destruction of children ghost cells */
        void fillBufferXi(double* buffer, int& counter, const int& lvl, const int& neighbour) const;
        void getBufferXi(double* buffer, int& counter, const int& lvl);
        void fillBufferSplit(bool* buffer, int& counter, const int& lvl, const int& neighbour) const;
        void getBufferSplit(bool* buffer, int& counter, const int& lvl);
        void fillNumberElementsToSendToNeighbour(int& numberElementsToSendToNeighbor, int& numberSlopesToSendToNeighbor, const int& lvl, const int& neighbour, int numberNeighboursOfCPUneighbour);
        void fillDataToSend(std::vector<double>& dataToSend, std::vector<int>& dataSplitToSend, const int& lvl) const;
        void getDataToReceiveAndRefine(std::vector<double>& dataToReceive, std::vector<int>& dataSplitToReceive, const int& lvl, Eos** eos, int& counter, int& counterSplit,
            const int& nbCellsY, const int& nbCellsZ, const std::vector<AddPhys*>& addPhys);
        void computeLoad(double& load, int lvl) const;
        void computeLvlMax(int& lvlMax) const;
        void clearExternalCellInterfaces(const int& nbCellsY, const int& nbCellsZ);
        void updatePointersInternalCellInterfaces();
        void updateNbCellsTotalAMR(int& nbCellsTotalAMR);

    protected:
      bool m_wall;                                                /*!< Bool indicating if the cell is a solid boundary (for immersed boundaries) */
      Phase** m_vecPhases;                                        /*!< Array of phases */
      Mixture* m_mixture;                                         /*!< Mixture */
      Transport* m_vecTransports;                                 /*!< Array of passive variables advected in the flow */
      Flux* m_cons;                                               /*!< Conservative variables */
      Transport* m_consTransports;                                /*!< Array of fluxes for advected variables */
      Element *m_element;                                         /*!< Pointer to corresponding geometrical mesh element */
      std::vector<CellInterface*> m_cellInterfaces;               /*!< Vector of cell-interface pointers */
      std::vector<QuantitiesAddPhys*> m_vecQuantitiesAddPhys;     /*!< Vector of pointers to the Quantities of Additional Physics of the cell */
     
      //Attributs pour methode AMR
      int m_lvl;                                                  /*!< Cell AMR level in the AMR tree */
      double m_xi;                                                /*!< Criteria for refine/unrefine cell */
      double m_consXi;                                            /*!< Buffer variable for Xi fluxes */
      bool m_split;                                               /*!< Indicator for splitted cell (Do I possess children ?) */
      std::vector<Cell*> m_childrenCells;                         /*!< Vector of children cells pointers */
      std::vector<CellInterface*> m_childrenInternalCellInterfaces; /*!< Vector of Internal children cell-interface pointers of the cell */

    private:
};

extern Model* model;                                      /*!< Pointer to model */
extern Gradient* gradient;                                /*!< Pointer to gradient method */
extern std::vector<Coord> gradRho;                        /*!< Gradient of density */
extern std::vector<Variable> variableDensity;             /*!< Variable name for density gradients */
extern std::vector<int> numeratorDefault;                 /*!< Default numerator (used for density gradients) */

#endif // CELL_H
