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

#ifndef ELEMENTNS_H
#define ELEMENTNS_H

//! \file      ElementNS.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "../../Element.h"

class FaceNS; //Predeclaration de la classe Face pour pouvoir inclure Face.h
#include "FaceNS.h"

class ElementNS : public Element
{
public:
  ElementNS();
  ElementNS(const int &typeGmsh, const int &numberNoeuds, const int &numberFaces, const int &typeVTK);
  virtual ~ElementNS();

  void construitElement(const int *numNoeuds, const Coord *noeuds, const int numberEntitePhysique, const int numberEntiteGeometrique, int &indexElement); /*Calcul des proprietes de l element*/
  void construitElementParallele(const Coord *noeuds); /*Calcul des proprietes de l element*/
  virtual void attributFaceLimite(const Coord *noeuds, FaceNS **faces, const int &indexMaxFaces){ Errors::errorMessage("attributFaceLimite non prevu pour le type d element demande"); };
  virtual void attributFaceCommunicante(const Coord *noeuds, FaceNS **faces, const int &indexMaxFaces, const int &numberNoeudsInternes){ Errors::errorMessage("attributFaceCommunicante non prevu pour le type d element demande"); };
  virtual void construitFaces(const Coord *noeuds, FaceNS **faces, int &indexMaxFaces, int** facesTemp, int* sommeNoeudsTemp) { Errors::errorMessage("construitFaces non prevu pour le type d element demande"); }; //Pour tests
  virtual void construitFacesSimplifie(int &iMax, int** facesTemp, int* sommeNoeudsTemp){ Errors::errorMessage("construitFacesSimplifie non prevu pour le type d element demande"); };
  virtual int compteFaceCommunicante(std::vector<int*> &faces, std::vector<int> &sommeNoeudsTemp){ Errors::errorMessage("compteFaceCommunicante non prevu pour le type d element demande"); return 0; };
  virtual int compteFaceCommunicante(int &iMax, int **faces, int *sommeNoeudsTemp){ Errors::errorMessage("compteFaceCommunicante non prevu pour le type d element demande"); return 0; };
  /*  void enleveCPUAutres(const int &numCPU);*/
  void enleveCPUAutres(std::vector<int> &numCPU);

  //Accesseurs
  void setIndex(int &index);
  void setAppartenancePhysique(int &appartenancePhysique);
  void setNumNoeud(int *numNoeuds);
  void setNumNoeud(int &noeud, int &numNoeud);
  void setIsFantome(bool isFantome);
  void setIsCommunicant(bool isCommunicant);
  void setAppartenanceCPU(const int *numCPU, const int &numberCPU);

  virtual const int& getIndex() const { return m_index; };
  const int& getNumberNoeuds() const { return m_numberNoeuds; };
  const int& getNumberFaces() const { return m_numberFaces; };
  const int& getTypeGmsh() const { return m_typeGmsh; };
  const int& getTypeVTK() const { return m_typeVTK; };
  const int& getNumNoeud(int &noeud) const { return m_numNoeuds[noeud]; };
  virtual const int& getAppartenancePhysique() const { return m_appartenancePhysique; };
  const int& getAppartenanceGeometrique() const { return m_appartenanceGeometrique; };
  const int& getCPU() const { return m_CPU; };
  const int& getNumberAutresCPU() const { return m_numberautresCPU; };
  const int& getAutreCPU(const int &autreCPU) const;
  void printInfo() const;

  const bool& isFantome() const { return m_isFantome; };
  const bool& isCommunicant() const { return m_isCommunicant; };

protected:
  virtual void computeVolume(const Coord *noeuds){};
  virtual void computeLCFL(const Coord *noeuds){};

  int m_index;

  int m_typeGmsh;
  int m_typeVTK;
  int m_numberNoeuds;
  int m_numberFaces;
  int m_appartenancePhysique;
  int m_appartenanceGeometrique;
  bool m_isFantome;
  bool m_isCommunicant;
  int m_CPU;              /*Number du CPU sur lequel l'element est present physiquement*/
  int m_numberautresCPU;
  int *m_autresCPU;       /*Number des CPU sur lesquels l'element est present en tant que fantome*/
  int *m_numNoeuds;       /*Correspondance avec le tableau de noeuds du mesh*/

};

#endif // ELEMENTNS_H