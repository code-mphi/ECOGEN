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

#ifndef ELEMENT_H
#define ELEMENT_H

//! \file      Element.h
//! \author    F. Petitpas, K.Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      December 20 2017

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../Maths/Coord.h"
#include "../Maths/GeometricObject.h"
#include "../Errors.h"
#include "../Tools.h"

class Element;

#include "Face.h"

class Element
{
public:
  Element();
  virtual ~Element();

  //Accesseurs
  void setCellAssociee(const int &numCell);
  Coord getPosition() const;
  double getLCFL() const;
  double getVolume() const;
  int getNumCellAssociee() const;

  virtual int getIndex() const { Errors::errorMessage("getIndex not available for requested element"); return 0; };
  virtual int getAppartenancePhysique() const { return 0; }; //!< default
  virtual void setVolume(const double &volume){ Errors::errorMessage("setVolume not available for requested element"); };
  virtual void setLCFL(const double &lCFL){ Errors::errorMessage("setlCFL not available for requested element"); };
  virtual void setPos(const double &X, const double &Y, const double &Z){ Errors::errorMessage("setPos not available for requested element"); };
  virtual void setPos(const Coord &pos) { Errors::errorMessage("setPos not available for requested element"); };
  virtual void setPosX(const double &X) { Errors::errorMessage("setPosX not available for requested element"); };
  virtual void setPosY(const double &Y) { Errors::errorMessage("setPosY not available for requested element"); };
  virtual void setPosZ(const double &Z) { Errors::errorMessage("setPosZ not available for requested element"); };
  virtual void setSize(const double &sizeX, const double &sizeY, const double &sizeZ) { Errors::errorMessage("setSize not available for requested element"); };
  virtual void setSize(const Coord &size) { Errors::errorMessage("setSize not available for requested element"); };
  
  void ecritPos(std::ofstream &fileStream, Axe axe);

  virtual void printInfo() const{ Errors::errorMessage("AfficheInfos not available for requested element"); };
  
  Coord vecteur(const Element *e); /*!< Cree un vecteur a partir des centers d elements */
  Coord vecteur(const Face *f);    /*!< Cree un vecteur entre center element et center d une face */

  double distance(const Element *e);  /*!< Calcul de la distance entre center et center d un autre element */
  double distanceX(const Element *e); /*!< Calcul de la distance selon x entre center et center d un autre element */
  double distanceY(const Element *e); /*!< Calcul de la distance selon y entre center et center d un autre element */
  double distanceZ(const Element *e); /*!< Calcul de la distance selon z entre center et center d un autre element */
  double distance(const Face *f);     /*!< Calcul de la distance entre center et center d une face */
  double distanceX(const Face *f);    /*!< Calcul de la distance selon x entre center et center d une face */
  double distanceY(const Face *f);    /*!< Calcul de la distance selon y entre center et center d une face */
  double distanceZ(const Face *f);    /*!< Calcul de la distance selon z entre center et center d une face */

  virtual double getSizeX() { Errors::errorMessage("getSizeX not available for requested element"); return 0; };
  virtual double getSizeY() { Errors::errorMessage("getSizeY not available for requested element"); return 0; };
  virtual double getSizeZ() { Errors::errorMessage("getSizeZ not available for requested element"); return 0; };
  virtual Coord getSize() { Errors::errorMessage("getSize not available for requested element"); return 0; };

  bool traverseObjet(const GeometricObject &objet) const;

  //Pour methode AMR
  virtual void creerElementChild() { Errors::errorMessage("creerElementsChildren not available for requested element"); };
  virtual Element* getElementChild(const int &numberChild) { Errors::errorMessage("getElementChild not available for requested element"); return 0; };
  virtual void finalizeElementsChildren() { Errors::errorMessage("finalizeElementsChildren not available for requested element"); };

protected:

  Coord m_position;       /*!< Position du center de l'element*/
  double m_volume;        /*!< Volume pour les elements 3D, Aire pour les elements 2D, Longueur pour les elements 1D*/
  double m_lCFL;          /*!< Longueur utile pour le compute du pas de temps*/
  int m_numCellAssociee;
};

#endif // ELEMENT_H