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

#ifndef MESHCARTESIANAMR_H
#define MESHCARTESIANAMR_H

#include "MeshCartesian.h"

class MeshCartesianAMR : public MeshCartesian
{
public:
  MeshCartesianAMR(double lX, int numberCellsX, double lY, int numberCellsY, double lZ, int numberCellsZ,
    std::vector<stretchZone> stretchX, std::vector<stretchZone> stretchY, std::vector<stretchZone> stretchZ,
		int lvlMax = 0, double criteriaVar = 1.e10, bool varRho = false, bool varP = false, bool varU = false, 
    bool varAlpha = false, double xiSplit = 1., double xiJoin = 1.);
  virtual ~MeshCartesianAMR();

  virtual int initializeGeometrie(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<Cell*>& cellsGhost, TypeMeshContainer<CellInterface*>& cellInterfaces,
    const int& restartSimulation, bool /*pretraitementParallele*/, std::string ordreCalcul);
  void initializeGeometrieAMR(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<Cell*>& cellsGhost, TypeMeshContainer<CellInterface*>& cellInterfaces, const int& restartSimulation, std::string ordreCalcul);
  void assignElementProperties(TypeMeshContainer<Cell*>& cells, std::vector<decomposition::Key<3>>& keys);
  void createCellInterfacesFacesAndGhostCells(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<Cell*>& cellsGhost, TypeMeshContainer<CellInterface*>& cellInterfaces, std::string ordreCalcul);
  virtual void procedureRaffinementInitialization(TypeMeshContainer<Cell*>* cellsLvl, TypeMeshContainer<Cell*>* cellsLvlGhost, TypeMeshContainer<CellInterface*>* cellInterfacesLvl,
		const std::vector<AddPhys*>& addPhys, int& nbCellsTotalAMR, std::vector<GeometricalDomain*>& domains,
		Eos** eos, const int& restartSimulation, std::string ordreCalcul);
  virtual void procedureRaffinement(TypeMeshContainer<Cell*>* cellsLvl, TypeMeshContainer<Cell*>* cellsLvlGhost, TypeMeshContainer<CellInterface*>* cellInterfacesLvl, const int& lvl,
    const std::vector<AddPhys*>& addPhys, int& nbCellsTotalAMR, Eos** eos);
  virtual std::string whoAmI() const;

  //Printing / Reading
  virtual void ecritHeaderPiece(std::ofstream& fileStream, TypeMeshContainer<Cell*>* cellsLvl) const;
  virtual void recupereNoeuds(std::vector<double>& jeuDonnees, std::vector<Cell*>* cellsLvl) const;
  virtual void recupereConnectivite(std::vector<double>& jeuDonnees, std::vector<Cell*>* cellsLvl) const;
  virtual void recupereOffsets(std::vector<double>& jeuDonnees, std::vector<Cell*>* cellsLvl) const;
  virtual void recupereTypeCell(std::vector<double>& jeuDonnees, std::vector<Cell*>* cellsLvl) const;
  virtual void recupereDonnees(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& jeuDonnees, const int var, int phase) const;
  virtual void setDataSet(std::vector<double>& jeuDonnees, TypeMeshContainer<Cell*>* cellsLvl, const int var, int phase) const;
  virtual void refineCellAndCellInterfaces(Cell* cell, const std::vector<AddPhys*>& addPhys, int& nbCellsTotalAMR);
  virtual void printDomainDecomposition(std::ofstream& fileStream);
  virtual void readDomainDecomposition(std::ifstream& fileStream);

  //Accesseurs
  virtual int getLvlMax() const { return m_lvlMax; };

	//Pour parallele
  virtual void initializePersistentCommunications(const TypeMeshContainer<Cell*>& cells, std::string ordreCalcul);
  virtual void finalizeParallele(const int& lvlMax);
  virtual void parallelLoadBalancingAMR(TypeMeshContainer<Cell*>* cellsLvl, TypeMeshContainer<Cell*>* cellsLvlGhost,
    TypeMeshContainer<CellInterface*>* cellInterfacesLvl, std::string ordreCalcul,
    const std::vector<AddPhys*>& addPhys, Eos** eos, int& nbCellsTotalAMR, bool init = false);
  virtual void computePotentialBalancing(TypeMeshContainer<Cell*>* cellsLvl, bool init, int lvl, bool& balance,
    std::vector<typename decomposition::Key<3>::value_type>& indicesSendStartGlobal, std::vector<typename decomposition::Key<3>::value_type>& indicesSendEndGlobal,
    std::vector<typename decomposition::Key<3>::value_type>& indicesReceiveStartGlobal, std::vector<typename decomposition::Key<3>::value_type>& indicesReceiveEndGlobal);
  virtual void balance(TypeMeshContainer<Cell*>* cellsLvl, TypeMeshContainer<Cell*>* cellsLvlGhost,
    TypeMeshContainer<CellInterface*>* cellInterfacesLvl, std::string ordreCalcul,
    const std::vector<AddPhys*>& addPhys, Eos** eos, int& nbCellsTotalAMR,
    std::vector<typename decomposition::Key<3>::value_type>& indicesSendStartGlobal, std::vector<typename decomposition::Key<3>::value_type>& indicesSendEndGlobal,
    std::vector<typename decomposition::Key<3>::value_type>& indicesReceiveStartGlobal, std::vector<typename decomposition::Key<3>::value_type>& indicesReceiveEndGlobal);

private:
  int m_lvlMax;                               //!<Niveau maximal sur l arbre AMR (si m_lvlMax = 0, pas d AMR)
	double m_criteriaVar;                       //!<Valeur du criteria a depasser sur la variation d'une variable pour le (de)raffinement (met xi=1.)
	bool m_varRho, m_varP, m_varU, m_varAlpha;  //!<Choix sur quelle variation on (de)raffine
	double m_xiSplit, m_xiJoin;                 //!<Valeur de xi pour split ou join les mailles
  decomposition::Decomposition m_decomp;      //!<Parallel domain decomposition based on keys

};

#endif // MESHCARTESIANAMR_H
