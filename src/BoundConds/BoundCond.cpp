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

//! \file      BoundCond.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "BoundCond.h"
#include <iostream>

using namespace std;

//***********************************************************************

BoundCond::BoundCond(){}

//***********************************************************************

BoundCond::BoundCond(int numPhysique) : m_numPhysique(numPhysique)
{}

//***********************************************************************

BoundCond::BoundCond(const BoundCond &Source) : m_numPhysique(Source.m_numPhysique)
{}

//***********************************************************************

BoundCond::~BoundCond() {}

//***********************************************************************

void BoundCond::initialize(Cell *cellLeft, Cell *cellRight)
{
  m_cellLeft = cellLeft;
  m_cellRight = NULL;
}

//***********************************************************************

void BoundCond::computeFlux(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type)
{
  this->solveRiemann(numberPhases, numberTransports, dtMax, globalLimiter, interfaceLimiter, globalVolumeFractionLimiter, interfaceVolumeFractionLimiter, type);
  this->subtractFlux(numberPhases, numberTransports, 1.); //Retrait du flux sur maille gauche
}

//***********************************************************************

void BoundCond::computeFluxAddPhys(const int &numberPhases, AddPhys &addPhys)
{
  addPhys.computeFluxAddPhysBoundary(this, numberPhases);
}

//***********************************************************************

void BoundCond::solveRiemann(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type)
{
  cellLeft->copyVec(m_cellLeft->getPhases(type), m_cellLeft->getMixture(type), m_cellLeft->getTransports(type));
  //Projection des velocities sur repere attache a la face
  cellLeft->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal(), numberPhases);
  //Calcul des variables etendus (Phases, Mixture, AddPhys)
  cellLeft->fulfillState();

  //Probleme de Riemann
  double dxLeft(m_cellLeft->getElement()->getLCFL());
  dxLeft = dxLeft*pow(2., (double)m_lvl);
  this->solveRiemannLimite(*cellLeft, numberPhases, dxLeft, dtMax);
  //Traitement des fonctions de transport (m_Sm connu : doit etre place apres l appel au Solveur de Riemann)
  if (numberTransports > 0) { this->solveRiemannTransportLimite(*cellLeft, numberTransports); }

  //Projection du flux sur le repere absolu
  m_mod->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
}

//****************************************************************************

int BoundCond::getNumPhys() const
{
  return m_numPhysique;
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void BoundCond::computeFluxXi()
{
  if ((m_cellLeft->getXi() > 0.05)) {
    m_cellLeft->addFluxXi(0.1);
  }
}

//****************************************************************************

void BoundCond::raffineBordExterne(const int &nbCellsY, const int &nbCellsZ, const double &dXParent, const double &dYParent,
  const double &dZParent, Cell *cellRef, const int &dim)
{
  //Le bord est une CL -> Creation des boundaries enfants
  double surfaceChild(pow(0.5, dim - 1.)*m_face->getSurface());
  double epsilon(1.e-6);
  int allocateSlopeLocal = 1;

  if (nbCellsZ == 1) {
    if (nbCellsY == 1) {

      //--------------------------------------------------
      //--------------------- Cas 1D ---------------------
      //--------------------------------------------------

      this->creerBordChild();
      m_boundariesChildren[0]->creerFaceChild(this);
      m_boundariesChildren[0]->getFace()->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
      m_boundariesChildren[0]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY(), m_face->getPos().getZ());
      m_boundariesChildren[0]->getFace()->setSize(m_face->getSize());
      if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
        //Bord number 1 (gauche)
        m_boundariesChildren[0]->initializeGauche(cellRef->getCellChild(0));
        cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[0]);
      }
      else {
        //Bord number 2 (droite)
        m_boundariesChildren[0]->initializeGauche(cellRef->getCellChild(1));
        cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[0]);
      }
      m_boundariesChildren[0]->associeModel(m_mod);
      m_boundariesChildren[0]->allocateSlopes(cellRef->getNumberPhases(), cellRef->getNumberTransports(), allocateSlopeLocal);
    }
    else {

      //--------------------------------------------------
      //--------------------- Cas 2D ---------------------
      //--------------------------------------------------

      //Creation des boundaries et faces enfants avec premiere initialization
      //----------------------------------------------------------------
      for (int i = 0; i < 2; i++) {
        this->creerBordChild();
        m_boundariesChildren[i]->creerFaceChild(this);
        m_boundariesChildren[i]->getFace()->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
      }

      //Face selon X
      //------------
      if (abs(m_face->getNormal().getX()) > epsilon) {
        //Cote gauche
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          for (int i = 0; i < 2; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(0));
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else {
              m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(2));
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
            }
            m_boundariesChildren[i]->getFace()->setSize(0., 0.5*m_face->getSizeY(), m_face->getSizeZ());
          }
        }
        //Cote droite
        else {
          for (int i = 0; i < 2; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else {
              m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
            }
            m_boundariesChildren[i]->getFace()->setSize(0., 0.5*m_face->getSizeY(), m_face->getSizeZ());
          }
        }
      }

      //Face selon Y
      //------------
      else {
        //Cote bas
        if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
          for (int i = 0; i < 2; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(0));
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else {
              m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
            }
            m_boundariesChildren[i]->getFace()->setSize(0.5*m_face->getSizeX(), 0., m_face->getSizeZ());
          }
        }
        //Cote haut
        else {
          for (int i = 0; i < 2; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(2));
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else {
              m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
            }
            m_boundariesChildren[i]->getFace()->setSize(0.5*m_face->getSizeX(), 0., m_face->getSizeZ());
          }
        }
      }

      //Association du model et des slopes
      //-----------------------------------
      for (int i = 0; i < 2; i++) {
        m_boundariesChildren[i]->associeModel(m_mod);
        m_boundariesChildren[i]->allocateSlopes(cellRef->getNumberPhases(), cellRef->getNumberTransports(), allocateSlopeLocal);
      }

    }
  }
  else {

    //--------------------------------------------------
    //--------------------- Cas 3D ---------------------
    //--------------------------------------------------

    //Creation des boundaries et faces enfants avec premiere initialization
    //----------------------------------------------------------------
    for (int i = 0; i < 4; i++) {
      this->creerBordChild();
      m_boundariesChildren[i]->creerFaceChild(this);
      m_boundariesChildren[i]->getFace()->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
      m_boundariesChildren[i]->getFace()->setSize(0.5*m_face->getSize());
    }

    //Face selon X
    //------------
    if (abs(m_face->getNormal().getX()) > epsilon) {
      //Cote gauche
      if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(0));
            cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(4));
            cellRef->getCellChild(4)->addBoundary(m_boundariesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(2));
            cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
          }
          //Fourth face
          else {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(6));
            cellRef->getCellChild(6)->addBoundary(m_boundariesChildren[i]);
          }
        }
      }
      //Cote droite
      else {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
            cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(5));
            cellRef->getCellChild(5)->addBoundary(m_boundariesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
            cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
          }
          //Fourth face
          else {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(7));
            cellRef->getCellChild(7)->addBoundary(m_boundariesChildren[i]);
          }
        }
      }
    }

    //Face selon Y
    //------------
    else if (abs(m_face->getNormal().getY()) > epsilon) {
      //Cote bas
      if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(0));
            cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
            cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(4));
            cellRef->getCellChild(4)->addBoundary(m_boundariesChildren[i]);
          }
          //Fourth face
          else {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(5));
            cellRef->getCellChild(5)->addBoundary(m_boundariesChildren[i]);
          }
        }
      }
      //Cote haut
      else {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(2));
            cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
            cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(6));
            cellRef->getCellChild(6)->addBoundary(m_boundariesChildren[i]);
          }
          //Fourth face
          else {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(7));
            cellRef->getCellChild(7)->addBoundary(m_boundariesChildren[i]);
          }
        }
      }
    }

    //Face selon Z
    //------------
    else {
      //Cote devant
      if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(0));
            cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
            cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(2));
            cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
          }
          //Fourth face
          else {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
            cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
          }
        }
      }
      //Cote derriere
      else {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(4));
            cellRef->getCellChild(4)->addBoundary(m_boundariesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(5));
            cellRef->getCellChild(5)->addBoundary(m_boundariesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(6));
            cellRef->getCellChild(6)->addBoundary(m_boundariesChildren[i]);
          }
          //Fourth face
          else {
            m_boundariesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(7));
            cellRef->getCellChild(7)->addBoundary(m_boundariesChildren[i]);
          }
        }
      }
    }

    //Association du model et des slopes
    //-----------------------------------
    for (int i = 0; i < 4; i++) {
      m_boundariesChildren[i]->associeModel(m_mod);
      m_boundariesChildren[i]->allocateSlopes(cellRef->getNumberPhases(), cellRef->getNumberTransports(), allocateSlopeLocal);
    }

  }
}

//***********************************************************************

void BoundCond::deraffineBordExterne(Cell *cellRef)
{
  //Dans le cas des CL, on parcourt toujours les boundaries parents mais on peut directement les deraffiner et mettre a jour la cell gauche.

  //Parcours les boundaries parents
  if (cellRef->getLvl() == m_lvl) {
    //cellRef est forcement la cell gauche, on deraffine le bord parent et on met a jour la cell gauche (de reference)
    for (unsigned int bordChild = 0; bordChild < m_boundariesChildren.size(); bordChild++) {
      cellRef->deleteBoundary(m_boundariesChildren[bordChild]);
    }
    this->deraffineBordsChildren();
  }
}

//****************************************************************************