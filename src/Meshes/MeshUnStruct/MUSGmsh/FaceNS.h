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

#ifndef FACENS_H
#define FACENS_H

#include "../../Face.h"

class ElementNS; //Predeclaration de la classe Element pour pouvoir inclure Element.h
#include "ElementNS.h"

class FaceNS : public Face
{
public:
  FaceNS();
  FaceNS(const int& numberNodes);
  virtual ~FaceNS();

  void construitFace(const Coord* nodes, const int& numNodeOther, ElementNS* elementNeighbor);
  bool faceExists(FaceNS** faces, const int& indexMaxFaces, int& indexFaceExiste);
  void addElementNeighbor(ElementNS* elementNeighbor);
  void addElementNeighborLimite(ElementNS* elementNeighbor);

  //Accesseurs
  ElementNS *getElementGauche() const;
  ElementNS *getElementDroite() const;
  void setEstLimite(const bool& estLimite);
  void setEstComm(const bool& estComm);
  const bool& getEstComm() const { return m_comm; };
  const int& getSumNumNodes() const { return m_sumNumNodes; };

  const int& getNumberNodes() const { return m_numberNodes; };
  const int& getNumNode(const int& numNode) const { return m_numNodes[numNode]; };
  void getInfoNodes(int* numNodes, int& sumNumNodes) const;
  const bool& getEstLimite() const { return m_limite; };
  void printNodes() const;
  virtual void printInfo() const;
  static int searchFace(int* face, int& sumNodes, int** arrayFaces, int* arraySumNodes, int numberNodes, int& indexMaxFaces); // Recherche si face appartient au array arrayFaces : renvoi le number ou -1 si absence
  static int searchFace(int* face, int& sumNodes, std::vector<int*> arrayFaces, std::vector<int> arraySumNodes, int numberNodes, int& indexMaxFaces); // Recherche si face appartient au array arrayFaces : renvoi le number ou -1 si absence. Utilise seulement sur ancienne version

protected:
  virtual void computeSurface(const Coord* /*nodes*/){};
  virtual void computeRepere(const Coord* /*nodes*/, const int& /*numNodeOther*/, ElementNS* /*elementNeighbor*/){};

  int m_numberNodes;   /*!< Number de nodes de la face */
  int* m_numNodes;     /*!< Number des node relatif au array de nodes global */
  int m_sumNumNodes;
  bool m_limite;
  bool m_comm;
  ElementNS* m_elementGauche;
  ElementNS* m_elementDroite;

};

#endif // FACENS_H

//Surcharge operateur externe a la classe car prends deux arguments
bool operator==(const FaceNS& a, const FaceNS& b);
bool operator!=(const FaceNS& a, const FaceNS& b);