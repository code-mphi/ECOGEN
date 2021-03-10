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

#ifndef MESH_H
#define MESH_H

#include <ctime>
#include <fstream>
#include <vector>
#include <list>
#include "../libTierces/tinyxml2.h"
#include "../Order1/CellGhost.h"
#include "../Order1/CellInterface.h"
#include "../Order2/CellInterfaceO2.h"
#include "../Order2/CellO2.h"
#include "../Order2/CellO2Ghost.h"
#include "../BoundConds/HeaderBoundCond.h"
#include "../Maths/Coord.h"
#include "../Parallel/Parallel.h"
#include "../AdditionalPhysics/HeaderQuantitiesAddPhys.h"
#include "../Maths/GeometricObject.h"
#include "../Sources/HeaderSources.h"

//! \class     Mesh
//! \brief     Abstract class for a mesh
//! \details   Can not be instanciated, depend on the mesh properties
class Mesh
{
public:
  Mesh();
  virtual ~Mesh();

  virtual void attributLimites(std::vector<BoundCond*>& boundCond) = 0;
  virtual int initializeGeometrie(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<Cell*>& cellsGhost, TypeMeshContainer<CellInterface*>& cellInterfaces,
    const int& restartSimulation, bool pretraitementParallele = true, std::string ordreCalcul = "FIRSTORDER") = 0; //!< renvoi le number de dimensions (1,2 ou 3)
  virtual std::string whoAmI() const { Errors::errorMessage("whoAmI pas prevu pour le mesh demande"); return 0; };

  //Accessors
  //---------
  const int& getGeometrie() const { return m_geometrie; };
  const int& getNumberCells() const { return m_numberCellsCalcul; };
  const int& getNumberCellsTotal() const { return m_numberCellsTotal; };
  const int& getNumberFaces() const { return m_numberFacesTotal; };
  const int& getNumFichier() const { return m_numFichier; };
  virtual double getdX() const { return 0; };
  virtual double getdY() const { return 0; };
  virtual double getdZ() const { return 0; };
  const TypeM& getType() const { return m_type; };
  virtual int getNumberCellsY() { return 0; };
  virtual int getNumberCellsZ() { return 0; };
  virtual int getLvlMax() const { return 0; };

  //Printing
  //--------
  void ecritSolutionGnuplot(std::vector<Cell*>* cellsLvl, std::ofstream& fileStream, GeometricObject* objet = 0) const;
  virtual void ecritHeaderPiece(std::ofstream& /*fileStream*/, std::vector<Cell*>* /*cellsLvl*/) const { Errors::errorMessage("ecritHeaderPiece non prevu pour mesh considere"); };
  virtual std::string recupereChaineExtent(bool /*global*/ = false) const { Errors::errorMessage("recupereChaineExtent non prevu pour mesh considere"); return 0; };
  virtual void recupereCoord(std::vector<double>& /*jeuDonnees*/, Axis /*axis*/) const { Errors::errorMessage("recupereCoord non prevu pour mesh considere"); };
  virtual void recupereNoeuds(std::vector<double>& /*jeuDonnees*/, std::vector<Cell*>* /*cellsLvl*/) const { Errors::errorMessage("recupereNoeuds non prevu pour mesh considere"); };
  virtual void recupereConnectivite(std::vector<double>& /*jeuDonnees*/, std::vector<Cell*>* /*cellsLvl*/) const { Errors::errorMessage("recupereConnectivite non prevu pour mesh considere"); };
  virtual void recupereOffsets(std::vector<double>& /*jeuDonnees*/, std::vector<Cell*>* /*cellsLvl*/) const { Errors::errorMessage("recupereOffsets non prevu pour mesh considere"); };
  virtual void recupereTypeCell(std::vector<double>& /*jeuDonnees*/, std::vector<Cell*>* /*cellsLvl*/) const { Errors::errorMessage("recupereTypeCell non prevu pour mesh considere"); };
  //! \brief     Extracting data for printing results
  //! \details   This method enable to extract a set of data for mixture or phase, scalar or vetor
  //! \param     cellsLvl         data structure containing pointer to cells
  //! \param     var              number of requested varaible to extract (>0 for scalar, <0 for vector)
  //! \param     phase            number of requested phase (-1 for mixture, -2 for transport, -3 for xi, -4 for gradient density mixture)
  //! \param     jeuDonnees       double vector containing the extracted data
  virtual void recupereDonnees(std::vector<Cell*>* /*cellsLvl*/, std::vector<double>& /*jeuDonnees*/, const int /*var*/, int /*phase*/) const { Errors::errorMessage("recupereDonnees non prevu pour mesh considere"); };
  //! \brief     Extracting data for printing results
  //! \details   This method enable to extract a set of data for mixture or phase, scalar or vetor
  //! \param     cellsLvl         data structure containing pointer to cells
  //! \param     var              number of requested varaible to extract (>0 for scalar, <0 for vector)
  //! \param     phase            number of requested phase (-1 for mixture, -2 for transport, -3 for xi, -4 for gradient density mixture)
  //! \param     jeuDonnees       double vector containing the extracted data
  virtual void setDataSet(std::vector<double>& /*jeuDonnees*/, std::vector<Cell*>* /*cellsLvl*/, const int /*var*/, int /*phase*/) const { Errors::errorMessage("setDataSet not available for requested mesh"); };
  virtual void refineCellAndCellInterfaces(Cell* /*cell*/, const std::vector<AddPhys*>& /*addPhys*/, Model* /*model*/, int& /*nbCellsTotalAMR*/) { Errors::errorMessage("refineCellAndCellInterfaces not available for requested mesh"); };
  //! \brief     Extracting absolute velocity for specific Moving Reference Frame computations
  //! \param     cellsLvl         data structure containing pointer to cells
  //! \param     sourceMRF        pointer to the corresponding MRF source
  //! \param     jeuDonnees       double vector containing the extracted data
  virtual void extractAbsVelocityMRF(std::vector<Cell*>* /*cellsLvl*/, std::vector<double>& /*jeuDonnees*/, Source* /*sourceMRF*/) const { Errors::errorMessage("extractAbsVeloxityMRF non prevu pour mesh considere"); };
  virtual void printDomainDecomposition(std::ofstream& /*fileStream*/) {};
  virtual void readDomainDecomposition(std::ifstream& /*fileStream*/) {};
  
  //Specific to AMR method
  //----------------------
  virtual void procedureRaffinementInitialization(std::vector<Cell*>* /*cellsLvl*/, TypeMeshContainer<Cell*>* /*cellsLvlGhost*/, std::vector<CellInterface*>* /*cellInterfacesLvl*/,
    const std::vector<AddPhys*>& /*addPhys*/, Model* /*model*/, int& nbCellsTotalAMR, std::vector<GeometricalDomain*>& /*domains*/, Eos** /*eos*/, const int& /*restartSimulation*/, std::string /*ordreCalcul*/,
    const int& /*numberPhases*/, const int& /*numberTransports*/) { nbCellsTotalAMR = m_numberCellsCalcul; };
  virtual void procedureRaffinement(std::vector<Cell*>* /*cellsLvl*/, TypeMeshContainer<Cell*>* /*cellsLvlGhost*/, std::vector<CellInterface*>* /*cellInterfacesLvl*/, const int& /*lvl*/,
    const std::vector<AddPhys*>& /*addPhys*/, Model* /*model*/, int& /*nbCellsTotalAMR*/, Eos** /*eos*/) {};

	//Specific for parallel
  //---------------------
	virtual void initializePersistentCommunications(const int numberPhases, const int numberTransports, const TypeMeshContainer<Cell*>& cells, std::string ordreCalcul);
	virtual void finalizeParallele(const int& lvlMax);
  virtual void parallelLoadBalancingAMR(std::vector<Cell*>* /*cellsLvl*/, TypeMeshContainer<Cell*>* /*cellsLvlGhost*/, std::vector<CellInterface*>* /*cellInterfacesLvl*/, std::string /*ordreCalcul*/,
    const int& /*numberPhases*/, const int& /*numberTransports*/, const std::vector<AddPhys*>& /*addPhys*/, Model* /*model*/, Eos** /*eos*/, int& /*nbCellsTotalAMR*/, bool /*init*/ = false) {};

  
protected:
  mutable int m_numFichier;

  int m_geometrie;                         /*indicator 1D/2D/3D*/
  int m_numberElements;                    /*Number d'elements au total (cells de computes internes de dimension n + elements limites de dimension n-1 + ghost cells de dimensions n)*/
  int m_numberFacesTotal;                  /*Number de faces entre deux cells ou entre une cell et une limite*/
  int m_numberCellsCalcul;                 /*Number de cells de compute internes au domain*/
	int m_numberCellsTotal;                  /*Cells de compute internes + cells fantomes dediees aux communications parallele*/

  TypeM m_type;

	int m_numberPhases;
	int m_numberTransports;
};
#endif // MESH_H
