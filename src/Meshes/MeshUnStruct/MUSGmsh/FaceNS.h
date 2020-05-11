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

#ifndef FACENS_H
#define FACENS_H

//! \file      FaceNS.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "../../Face.h"

class ElementNS; //Predeclaration de la classe Element pour pouvoir inclure Element.h
#include "ElementNS.h"

class FaceNS : public Face
{
public:
  FaceNS();
  FaceNS(const int &numberNoeuds);
  virtual ~FaceNS();

  void construitFace(const Coord *noeuds, const int &numNoeudAutre, ElementNS *elementVoisin);
  bool faceExiste(FaceNS **faces, const int &indexMaxFaces, int &indexFaceExiste);
  void ajouteElementVoisin(ElementNS *elementVoisin);
  void ajouteElementVoisinLimite(ElementNS *elementVoisin);

  //Accesseurs
  ElementNS *getElementGauche() const;
  ElementNS *getElementDroite() const;
  void setEstLimite(const bool &estLimite);
  void setEstComm(const bool &estComm);
  const bool& getEstComm() const { return m_comm; };
  const int& getSommeNumNoeuds() const { return m_sommeNumNoeuds; };

  const int& getNumberNoeuds() const { return m_numberNoeuds; };
  const int& getNumNoeud(const int &numNoeud) const { return m_numNoeuds[numNoeud]; };
  void getInfoNoeuds(int *numNoeuds, int &sommeNumNoeuds) const;
  const bool& getEstLimite() const { return m_limite; };
  void afficheNoeuds() const;
  virtual void printInfo() const;
  static int rechercheFace(int *face, int &sommeNoeuds, int **tableauFaces, int *tableauSommeNoeuds, int numberNoeuds, int &indexMaxFaces); // Recherche si face appartient au tableau tableauFaces : renvoi le number ou -1 si absence
  static int rechercheFace(int *face, int &sommeNoeuds, std::vector<int*> tableauFaces, std::vector<int> tableauSommeNoeuds, int numberNoeuds, int &indexMaxFaces); // Recherche si face appartient au tableau tableauFaces : renvoi le number ou -1 si absence. Utilise seulement sur ancienne version

protected:
  virtual void computeSurface(const Coord *noeuds){};
  virtual void computeRepere(const Coord *noeuds, const int &numNoeudAutre, ElementNS *elementVoisin){};

  int m_numberNoeuds;   /*!< Number de noeuds de la face */
  int *m_numNoeuds;     /*!< Number des noeud relatif au tableau de noeuds global */
  int m_sommeNumNoeuds;
  bool m_limite;
  bool m_comm;
  ElementNS *m_elementGauche;
  ElementNS *m_elementDroite;

};

#endif // FACENS_H

//Surcharge operateur externe a la classe car prends deux arguments
bool operator==(const FaceNS &a, const FaceNS &b);
bool operator!=(const FaceNS &a, const FaceNS &b);