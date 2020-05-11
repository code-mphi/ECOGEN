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

#ifndef MESHCARTESIAN_H
#define MESHCARTESIAN_H

//! \file      MeshCartesian.h
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include "Mesh.h"
#include "ElementCartesian.h"
#include "FaceCartesian.h"
#include "stretchZone.h"
#include "../Parallel/decomposition.hpp"

class MeshCartesian : public Mesh
{
public:
  MeshCartesian(double lX, int numberCellsX, double lY, int numberCellsY, double lZ, int numberCellsZ,
    std::vector<stretchZone> stretchX, std::vector<stretchZone> stretchY, std::vector<stretchZone> stretchZ);
  virtual ~MeshCartesian();

  virtual void attributLimites(std::vector<BoundCond*> &boundCond);
  void recupereIJK(const int &index, int &i, int &j, int &k) const;
  void construitIGlobal(const int &i, const int &j, const int &k, int &index) const;
  virtual int initializeGeometrie(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost, TypeMeshContainer<CellInterface *> &cellInterfaces,
    const int &restartSimulation, bool pretraitementParallele, std::string ordreCalcul);
  void meshStretching();
  void initializeGeometrieMonoCpu(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<CellInterface *> &cellInterfaces, std::string ordreCalcul);
  void initializeGeometrieParallele(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost, TypeMeshContainer<CellInterface *> &cellInterfaces, std::string ordreCalcul);
  void decoupageParallele(std::string ordreCalcul, TypeMeshContainer<Cell *> &cells);
  virtual std::string whoAmI() const;

  //Accessors
  //---------
  virtual int getNumberCellsY() { return m_numberCellsY; };
  virtual int getNumberCellsZ() { return m_numberCellsZ; };

  //Printing / Reading
  //------------------
  virtual std::string recupereChaineExtent(int localRank, bool global = false) const;
  virtual void recupereCoord(TypeMeshContainer<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, Axis axis) const;
  virtual void recupereDonnees(TypeMeshContainer<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, const int var, int phase) const;
  virtual void setDataSet(std::vector<double> &jeuDonnees, TypeMeshContainer<Cell *> *cellsLvl, const int var, int phase) const;

protected:
  TypeMeshContainer<Element *> m_elements; //!<Vector of element objects: Contains geometrical attributes
  TypeMeshContainer<Face *> m_faces;       //!<Vector of face objects (between two elements or at boundaries): Contains geometrical attributes

  double m_lX;
  double m_lY;
  double m_lZ;
  int m_numberCellsX;
  int m_numberCellsY;
  int m_numberCellsZ;
  int m_numberCellsXGlobal;
  int m_numberCellsYGlobal;
  int m_numberCellsZGlobal;
  std::vector<double> m_dXi;                 /*!< Array of the lengths of the cells in the x-direction */
  std::vector<double> m_dYj;                 /*!< Array of the lengths of the cells in the y-direction */
  std::vector<double> m_dZk;                 /*!< Array of the lengths of the cells in the z-direction */
  std::vector<double> m_posXi;               /*!< Array of the positions of the cells in the x-direction */
  std::vector<double> m_posYj;               /*!< Array of the positions of the cells in the y-direction */
  std::vector<double> m_posZk;               /*!< Array of the positions of the cells in the z-direction */
  std::vector<stretchZone> m_stretchX;   /*!< Array of stretch zones in the x-direction */
  std::vector<stretchZone> m_stretchY;   /*!< Array of stretch zones in the y-direction */
  std::vector<stretchZone> m_stretchZ;   /*!< Array of stretch zones in the z-direction */
  int m_numberCpuX;              /*!< Optimal number of processors in the x-direction */
  int m_numberCpuY;              /*!< Optimal number of processors in the y-direction */
  int m_numberCpuZ;              /*!< Optimal number of processors in the z-direction */
  int m_CpuCoordX;               /*!< X-coordinate of the current CPU */
  int m_CpuCoordY;               /*!< Y-coordinate of the current CPU */
  int m_CpuCoordZ;               /*!< Z-coordinate of the current CPU */
  int m_offsetX;                 /*!< Offset in the x-direction of the current CPU for the array of cell lenghts and cell positions */
  int m_offsetY;                 /*!< Offset in the y-direction of the current CPU for the array of cell lenghts and cell positions */
  int m_offsetZ;                 /*!< Offset in the z-direction of the current CPU for the array of cell lenghts and cell positions */

  int m_numberBoundCondInit;
  BoundCond *m_limXm;
  BoundCond *m_limXp;
  BoundCond *m_limYm;
  BoundCond *m_limYp;
  BoundCond *m_limZm;
  BoundCond *m_limZp;
};

#endif // MESHCARTESIAN_H
