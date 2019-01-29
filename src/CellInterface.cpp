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

//! \file      CellInterface.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      December 20 2017

#include "CellInterface.h"
#include <iostream>

using namespace std;

//Utile pour la resolution des problemes de Riemann
Cell *cellLeft;
Cell *cellRight;

//***********************************************************************

CellInterface::CellInterface() : m_mod(0), m_cellLeft(0), m_cellRight(0), m_face(0), m_boundariesChildren(0)
{
  m_lvl = 0;
}

//***********************************************************************

CellInterface::CellInterface(int lvl) : m_mod(0), m_cellLeft(0), m_cellRight(0), m_face(0), m_boundariesChildren(0)
{
  m_lvl = lvl;
}

//***********************************************************************

CellInterface::~CellInterface()
{
  for (unsigned int i = 0; i < m_boundariesChildren.size(); i++) {
    m_boundariesChildren[i]->finalizeFace();
    delete m_boundariesChildren[i];
  }
  m_boundariesChildren.clear();
}

//***********************************************************************

void CellInterface::initialize(Cell *cellLeft, Cell *cellRight)
{
  m_cellLeft = cellLeft;
  m_cellRight = cellRight;
}

//***********************************************************************

void CellInterface::initializeGauche(Cell *cellLeft)
{
  m_cellLeft = cellLeft;
}

//***********************************************************************

void CellInterface::initializeDroite(Cell *cellRight)
{
  m_cellRight = cellRight;
}

//***********************************************************************

void CellInterface::setFace(Face *face)
{
  m_face = face;
}

//***********************************************************************

void CellInterface::computeFlux(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type)
{
  this->solveRiemann(numberPhases, numberTransports, dtMax, globalLimiter, interfaceLimiter, globalVolumeFractionLimiter, interfaceVolumeFractionLimiter, type);

  if (m_cellLeft->getLvl() == m_cellRight->getLvl()) {     //CoefAMR = 1 pour les deux
    this->addFlux(numberPhases, numberTransports, 1.);      //Ajout du flux sur maille droite
    this->subtractFlux(numberPhases, numberTransports, 1.);     //Retrait du flux sur maille gauche
  }
  else if (m_cellLeft->getLvl() > m_cellRight->getLvl()) { //CoefAMR = 1 pour la gauche et 0.5 pour la droite
    this->addFlux(numberPhases, numberTransports, 0.5);     //Ajout du flux sur maille droite
    this->subtractFlux(numberPhases, numberTransports, 1.);     //Retrait du flux sur maille gauche
  }
  else {                                                      //CoefAMR = 0.5 pour la gauche et 1 pour la droite
    this->addFlux(numberPhases, numberTransports, 1.);      //Ajout du flux sur maille droite
    this->subtractFlux(numberPhases, numberTransports, 0.5);    //Retrait du flux sur maille gauche
  }
}

//***********************************************************************

void CellInterface::computeFluxAddPhys(const int &numberPhases, AddPhys &addPhys)
{
  addPhys.computeFluxAddPhys(this, numberPhases);
}

//***********************************************************************

void CellInterface::solveRiemann(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type)
{
  //Projection des velocities sur repere attache a la face
  m_cellLeft->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal(), numberPhases);
  m_cellRight->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal(), numberPhases);

  //Probleme de Riemann
  double dxLeft(m_cellLeft->getElement()->getLCFL());
  double dxRight(m_cellRight->getElement()->getLCFL());
  dxLeft = dxLeft*pow(2., (double)m_lvl);
  dxRight = dxRight*pow(2., (double)m_lvl);
  m_mod->solveRiemannIntern(*m_cellLeft, *m_cellRight, numberPhases, dxLeft, dxRight, dtMax);
  //Traitement des fonctions de transport (m_Sm connu : doit etre place apres l appel au Solveur de Riemann)
  if (numberTransports > 0) { m_mod->solveRiemannTransportIntern(*m_cellLeft, *m_cellRight, numberTransports); }

  //Projection du flux sur le repere absolu
  m_mod->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
  m_cellLeft->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal(), numberPhases);
  m_cellRight->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal(), numberPhases);
}

//***********************************************************************

void CellInterface::addFlux(const int &numberPhases, const int &numberTransports, const double &coefAMR)
{
  double volume(m_cellRight->getElement()->getVolume());
  double surface(m_face->getSurface());
  double coefA(surface / volume); //pas de "pas de temps"
  coefA = coefA*coefAMR;
  m_cellRight->getCons()->addFlux(coefA, numberPhases);
  m_cellRight->getCons()->addNonCons(coefA, m_cellRight, numberPhases);
  double sM(m_mod->getSM());
  for (int k = 0; k < numberTransports; k++) {
    m_cellRight->getConsTransport(k)->addFlux(coefA, k);
    m_cellRight->getConsTransport(k)->addNonCons(coefA, m_cellRight->getTransport(k).getValue(), sM);
  }
}

//***********************************************************************

void CellInterface::subtractFlux(const int &numberPhases, const int &numberTransports, const double &coefAMR)
{
  double volume(m_cellLeft->getElement()->getVolume());
  double surface(m_face->getSurface());
  double coefA(surface / volume); //pas de "pas de temps"
  coefA = coefA*coefAMR;
  m_cellLeft->getCons()->subtractFlux(coefA, numberPhases);
  m_cellLeft->getCons()->subtractNonCons(coefA, m_cellLeft, numberPhases);
  double sM(m_mod->getSM());
  for (int k = 0; k < numberTransports; k++) {
    m_cellLeft->getConsTransport(k)->subtractFlux(coefA, k);
    m_cellLeft->getConsTransport(k)->subtractNonCons(coefA, m_cellLeft->getTransport(k).getValue(), sM);
  }
}

//***********************************************************************

double CellInterface::distance(Cell *c)
{
  return m_face->distance(c->getElement());
}

//***********************************************************************
void CellInterface::EffetsSurface1D(const int &numberPhases)
{
  if (m_cellRight != NULL){ m_cellRight->getCons()->addTuyere1D(m_face->getNormal(), m_face->getSurface(), m_cellRight, numberPhases); }
  m_cellLeft->getCons()->subtractTuyere1D(m_face->getNormal(), m_face->getSurface(), m_cellLeft, numberPhases);
}

//***********************************************************************

void CellInterface::associeModel(Model *model)
{
  m_mod = model;
}

//***********************************************************************

Face *CellInterface::getFace()
{
  return m_face;
}

//***********************************************************************

Model *CellInterface::getMod() const
{
  return m_mod;
}

//***********************************************************************

Cell *CellInterface::getCellGauche() const
{
  return m_cellLeft;
}

//***********************************************************************

Cell *CellInterface::getCellDroite() const
{
  return m_cellRight;
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void CellInterface::computeXi(const double &criteriaVar, const bool &varRho, const bool &varP, const bool &varU, const bool &varAlpha)
{
  if (varRho) { this->computeCritereAMR(criteriaVar, "RHO"); }
  if (varP) {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { this->computeCritereAMR(criteriaVar, "P"); }
  }
  if (varU) {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { this->computeCritereAMR(criteriaVar, "u"); }
  }
  if (varAlpha) {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { this->computeCritereAMR(criteriaVar, "ALPHA", 1); }
  }
}

//***********************************************************************

void CellInterface::computeCritereAMR(const double &criteriaVar, string nameVariable, int num)
{
  double valueMin, variation, cd, cg;
  // Recuperation des values de la variable en question a gauche et a droite
  cg = m_cellLeft->selectScalar(nameVariable, num);
  cd = m_cellRight->selectScalar(nameVariable, num);

  // Valeur de la variation
  valueMin = min(abs(cd), abs(cg));
  if (valueMin < 1.e-2) { valueMin = 1.e-2; } //Utile pour alpha (quasi-seulement) ou velocity
  variation = abs(cd - cg) / valueMin;

  //Mise a jour de xi si la variation est superieure au criteria
  if (variation >= criteriaVar) {
    m_cellLeft->setXi(1.);
    m_cellRight->setXi(1.);
  }
}

//***********************************************************************

void CellInterface::computeFluxXi()
{
  if (m_cellLeft->getLvl() > m_cellRight->getLvl()) {
    if ((m_cellLeft->getXi() > 0.05)) {
      m_cellLeft->addFluxXi(0.1);
      m_cellRight->addFluxXi(0.1);
    }
  }
  else if (m_cellLeft->getLvl() < m_cellRight->getLvl()) {
    if ((m_cellRight->getXi() > 0.05)) {
      m_cellLeft->addFluxXi(0.1);
      m_cellRight->addFluxXi(0.1);
    }
  }
  else {
    if ((m_cellLeft->getXi() > 0.05) || (m_cellRight->getXi() > 0.05)) {
      m_cellLeft->addFluxXi(0.1);
      m_cellRight->addFluxXi(0.1);
    }
  }
}

//***********************************************************************

void CellInterface::creerBordChild()
{
  m_boundariesChildren.push_back(new CellInterface(m_lvl + 1));
}

//***********************************************************************

void CellInterface::creerBordChildInterne(const int &lvl, vector<CellInterface*> *childrenInternalBoundaries)
{
  (*childrenInternalBoundaries).push_back(new CellInterface(lvl + 1));
}

//***********************************************************************

void CellInterface::creerFaceChild(CellInterface *bordParent)
{
  m_face = bordParent->m_face->creerNouvelleFace();
}

//***********************************************************************

void CellInterface::raffineBordExterne(const int &nbCellsY, const int &nbCellsZ, const double &dXParent, const double &dYParent, const double &dZParent, Cell *cellRef, const int &dim)
{
  //La creation des boundaries enfants n'est pas systematique, on regarde d'abord si ces boundaries enfants ne sont pas deja crees.
  //Dans tous les cas on re-attribut les liaisons cells/boundaries.

  double epsilon(1.e-6);
  int allocateSlopeLocal(1);
  double surfaceChild(pow(0.5,dim-1.)*m_face->getSurface());

  if (nbCellsZ == 1) {
    if (nbCellsY == 1) {

      //--------------------------------------------------
      //--------------------- Cas 1D ---------------------
      //--------------------------------------------------

      //Bord pas encore split -> Creation des boundaries enfants
      //---------------------------------------------------
      if (m_lvl == cellRef->getLvl()) {

        this->creerBordChild();
        m_boundariesChildren[0]->m_face = m_face->creerNouvelleFace();
        m_boundariesChildren[0]->m_face->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
        m_boundariesChildren[0]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY(), m_face->getPos().getZ());
        m_boundariesChildren[0]->m_face->setSize(m_face->getSize());
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          //Bord number 1 (gauche)
          m_boundariesChildren[0]->initializeGauche(m_cellLeft);
          m_boundariesChildren[0]->initializeDroite(cellRef->getCellChild(0));
          m_cellLeft->addBoundary(m_boundariesChildren[0]);
          cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[0]);
        }
        else {
          //Bord number 2 (droite)
          m_boundariesChildren[0]->initializeGauche(cellRef->getCellChild(1));
          m_boundariesChildren[0]->initializeDroite(m_cellRight);
          cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[0]);
          m_cellRight->addBoundary(m_boundariesChildren[0]);
        }
        m_boundariesChildren[0]->associeModel(m_mod);
        m_boundariesChildren[0]->allocateSlopes(cellRef->getNumberPhases(), cellRef->getNumberTransports(), allocateSlopeLocal);
      }

      //Bord deja split -> on met seulement a jour les liaisons cells/boundaries
      //----------------------------------------------------------------------
      else {
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          //Bord number 1 (gauche)
          m_cellRight = cellRef->getCellChild(0);
          cellRef->getCellChild(0)->addBoundary(this);
        }
        else {
          //Bord number 2 (droite)
          m_cellLeft = cellRef->getCellChild(1);
          cellRef->getCellChild(1)->addBoundary(this);
        }
      }
    }
    else {

      //--------------------------------------------------
      //--------------------- Cas 2D ---------------------
      //--------------------------------------------------

      //Bord pas encore split -> Creation des boundaries enfants
      //---------------------------------------------------
      if (m_lvl == cellRef->getLvl()) {

        //Creation des boundaries et faces enfants avec premiere initialization
        //----------------------------------------------------------------
        for (int i = 0; i < 2; i++) {
          this->creerBordChild();
          m_boundariesChildren[i]->m_face = m_face->creerNouvelleFace();
          m_boundariesChildren[i]->m_face->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
        }

        //Face selon X
        //------------
        if (abs(m_face->getNormal().getX()) > epsilon) {
          //Cote gauche
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(m_cellLeft);
                m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(0));
                m_cellLeft->addBoundary(m_boundariesChildren[i]);
                cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
              }
              //Second face
              else {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(m_cellLeft);
                m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(2));
                m_cellLeft->addBoundary(m_boundariesChildren[i]);
                cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
              }
              m_boundariesChildren[i]->m_face->setSize(0., 0.5*m_face->getSizeY(), m_face->getSizeZ());
            }
          }
          //Cote droite
          else {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
                m_boundariesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
                m_cellRight->addBoundary(m_boundariesChildren[i]);
              }
              //Second face
              else {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
                m_boundariesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
                m_cellRight->addBoundary(m_boundariesChildren[i]);
              }
              m_boundariesChildren[i]->m_face->setSize(0., 0.5*m_face->getSizeY(), m_face->getSizeZ());
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
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(m_cellLeft);
                m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(0));
                m_cellLeft->addBoundary(m_boundariesChildren[i]);
                cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
              }
              //Second face
              else {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(m_cellLeft);
                m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(1));
                m_cellLeft->addBoundary(m_boundariesChildren[i]);
                cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
              }
              m_boundariesChildren[i]->m_face->setSize(0.5*m_face->getSizeX(), 0., m_face->getSizeZ());
            }
          }
          //Cote haut
          else {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(2));
                m_boundariesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
                m_cellRight->addBoundary(m_boundariesChildren[i]);
              }
              //Second face
              else {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
                m_boundariesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
                m_cellRight->addBoundary(m_boundariesChildren[i]);
              }
              m_boundariesChildren[i]->m_face->setSize(0.5*m_face->getSizeX(), 0., m_face->getSizeZ());
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

      //Bord deja split -> on met seulement a jour les liaisons cells/boundaries
      //----------------------------------------------------------------------
      else {

        //Face selon X
        //------------
        if (abs(m_face->getNormal().getX()) > epsilon) {
          //Cote gauche
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            //First face
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            //Second face
            else {
              m_cellRight = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
          }
          //Cote droite
          else {
            //First face
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              m_cellLeft = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
            //Second face
            else {
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
          }
        }

        //Face selon Y
        //------------
        else {
          //Cote bas
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            //First face
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            //Second face
            else {
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
          }
          //Cote haut
          else {
            //First face
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              m_cellLeft = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
            //Second face
            else {
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
          }
        }
      }
    }
  }
  else {

    //--------------------------------------------------
    //--------------------- Cas 3D ---------------------
    //--------------------------------------------------

    //Bord pas encore split -> Creation des boundaries enfants
    //---------------------------------------------------
    if (m_lvl == cellRef->getLvl()) {

      //Creation des boundaries et faces enfants avec premiere initialization
      //----------------------------------------------------------------
      for (int i = 0; i < 4; i++) {
        this->creerBordChild();
        m_boundariesChildren[i]->m_face = m_face->creerNouvelleFace();
        m_boundariesChildren[i]->m_face->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
        m_boundariesChildren[i]->m_face->setSize(0.5*m_face->getSize());
      }

      //Face selon X
      //------------
      if (abs(m_face->getNormal().getX()) > epsilon) {
        //Cote gauche
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(0));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(4));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(4)->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(2));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(6));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(6)->addBoundary(m_boundariesChildren[i]);
            }
          }
        }
        //Cote droite
        else {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(5));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(5)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(7));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(7)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
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
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(0));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(1));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(4));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(4)->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(5));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(5)->addBoundary(m_boundariesChildren[i]);
            }
          }
        }
        //Cote haut
        else {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(2));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(6));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(6)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(7));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(7)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
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
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(0));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(1));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(2));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(3));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
            }
          }
        }
        //Cote derriere
        else {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(4));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(4)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(5));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(5)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(6));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(6)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(7));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(7)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
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

    //Bord deja split -> on met seulement a jour les liaisons cells/boundaries
    //----------------------------------------------------------------------
    else {

      //Face selon X
      //------------
      if (abs(abs(m_face->getNormal().getX()) - 1.) < epsilon) {
        //Cote gauche
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //First face
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            else {
              //Second face
              m_cellRight = cellRef->getCellChild(4);
              cellRef->getCellChild(4)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Third face
              m_cellRight = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellRight = cellRef->getCellChild(6);
              cellRef->getCellChild(6)->addBoundary(this);
            }
          }
        }
        //Cote droite
        else {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //First face
              m_cellLeft = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
            else {
              //Second face
              m_cellLeft = cellRef->getCellChild(5);
              cellRef->getCellChild(5)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Third face
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellLeft = cellRef->getCellChild(7);
              cellRef->getCellChild(7)->addBoundary(this);
            }
          }
        }
      }

      //Face selon Y
      //------------
      else if (abs(abs(m_face->getNormal().getY()) - 1.) < epsilon) {
        //Cote bas
        if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
          if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //First face
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            else {
              //Second face
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Third face
              m_cellRight = cellRef->getCellChild(4);
              cellRef->getCellChild(4)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellRight = cellRef->getCellChild(5);
              cellRef->getCellChild(5)->addBoundary(this);
            }
          }
        }
        //Cote haut
        else {
          if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //First face
              m_cellLeft = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
            else {
              //Second face
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Third face
              m_cellLeft = cellRef->getCellChild(6);
              cellRef->getCellChild(6)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellLeft = cellRef->getCellChild(7);
              cellRef->getCellChild(7)->addBoundary(this);
            }
          }
        }
      }

      //Face selon Z
      //------------
      else {
        //Cote devant
        if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //First face
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            else {
              //Second face
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Third face
              m_cellRight = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellRight = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
          }
        }
        //Cote derriere
        else {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //First face
              m_cellLeft = cellRef->getCellChild(4);
              cellRef->getCellChild(4)->addBoundary(this);
            }
            else {
              //Second face
              m_cellLeft = cellRef->getCellChild(5);
              cellRef->getCellChild(5)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Third face
              m_cellLeft = cellRef->getCellChild(6);
              cellRef->getCellChild(6)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellLeft = cellRef->getCellChild(7);
              cellRef->getCellChild(7)->addBoundary(this);
            }
          }
        }
      }
    }
  }
}


//***********************************************************************

void CellInterface::raffineBordExterneGhost(const int &nbCellsY, const int &nbCellsZ, const double &dXParent, const double &dYParent,
	const double &dZParent, Cell *cellRef, const int &dim)
{
	//La creation des boundaries enfants n'est pas systematique, on regarde d'abord si ces boundaries enfants ne sont pas deja crees.
	//Dans tous les cas on re-attribut les liaisons cells/boundaries.

	int allocateSlopeLocal(1);
  double surfaceChild(pow(0.5, dim - 1.)*m_face->getSurface());

	if (nbCellsZ == 1) {
		if (nbCellsY == 1) {

			//--------------------------------------------------
			//--------------------- Cas 1D ---------------------
			//--------------------------------------------------

			//Bord pas encore split -> Creation des boundaries enfants
			//---------------------------------------------------
			if (m_lvl == cellRef->getLvl()) {

				this->creerBordChild();
				m_boundariesChildren[0]->m_face = m_face->creerNouvelleFace();
				m_boundariesChildren[0]->m_face->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
				m_boundariesChildren[0]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY(), m_face->getPos().getZ());
        m_boundariesChildren[0]->m_face->setSize(m_face->getSize());
				if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
					//Bord number 1 (gauche)
					m_boundariesChildren[0]->initializeGauche(m_cellLeft);
					m_boundariesChildren[0]->initializeDroite(cellRef->getCellChild(0));
					m_cellLeft->addBoundary(m_boundariesChildren[0]);
					cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[0]);
				}
				else {
					//Bord number 2 (droite)
					m_boundariesChildren[0]->initializeGauche(cellRef->getCellChild(0));
					m_boundariesChildren[0]->initializeDroite(m_cellRight);
					cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[0]);
					m_cellRight->addBoundary(m_boundariesChildren[0]);
				}
				m_boundariesChildren[0]->associeModel(m_mod);
				m_boundariesChildren[0]->allocateSlopes(cellRef->getNumberPhases(), cellRef->getNumberTransports(), allocateSlopeLocal);
			}

			//Bord deja split -> on met seulement a jour les liaisons cells/boundaries
			//----------------------------------------------------------------------
			else {
				if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
					//Bord number 1 (gauche)
					m_cellRight = cellRef->getCellChild(0);
					cellRef->getCellChild(0)->addBoundary(this);
				}
				else {
					//Bord number 2 (droite)
					m_cellLeft = cellRef->getCellChild(0);
					cellRef->getCellChild(0)->addBoundary(this);
				}
			}
		}
		else {

			//--------------------------------------------------
			//--------------------- Cas 2D ---------------------
			//--------------------------------------------------

			//Bord pas encore split -> Creation des boundaries enfants
			//---------------------------------------------------
			if (m_lvl == cellRef->getLvl()) {

				//Creation des boundaries et faces enfants avec premiere initialization
				//----------------------------------------------------------------
				for (int i = 0; i < 2; i++) {
					this->creerBordChild();
					m_boundariesChildren[i]->m_face = m_face->creerNouvelleFace();
					m_boundariesChildren[i]->m_face->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
				}

				//Face in the x-direction
				//-----------------------
        if (abs(m_face->getNormal().getX()) > 0.99) {
          //Left side
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(m_cellLeft);
                m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(0));
                m_cellLeft->addBoundary(m_boundariesChildren[i]);
                cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
              }
              //Second face
              else {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(m_cellLeft);
                m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(1));
                m_cellLeft->addBoundary(m_boundariesChildren[i]);
                cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
              }
              m_boundariesChildren[i]->m_face->setSize(0., 0.5*m_face->getSizeY(), m_face->getSizeZ());
            }
          }
          //Right side
          else {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(0));
                m_boundariesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
                m_cellRight->addBoundary(m_boundariesChildren[i]);
              }
              //Second face
              else {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
                m_boundariesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
                m_cellRight->addBoundary(m_boundariesChildren[i]);
              }
              m_boundariesChildren[i]->m_face->setSize(0., 0.5*m_face->getSizeY(), m_face->getSizeZ());
            }
          }
        }
        //Face in the y-direction
        //-----------------------
        else if (abs(m_face->getNormal().getY()) > 0.99) {
          //Bottom side
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(m_cellLeft);
                m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(0));
                m_cellLeft->addBoundary(m_boundariesChildren[i]);
                cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
              }
              //Second face
              else {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(m_cellLeft);
                m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(1));
                m_cellLeft->addBoundary(m_boundariesChildren[i]);
                cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
              }
              m_boundariesChildren[i]->m_face->setSize(0.5*m_face->getSizeX(), 0., m_face->getSizeZ());
            }
          }
          //Top side
          else {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(0));
                m_boundariesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
                m_cellRight->addBoundary(m_boundariesChildren[i]);
              }
              //Second face
              else {
                m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
                m_boundariesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
                m_cellRight->addBoundary(m_boundariesChildren[i]);
              }
              m_boundariesChildren[i]->m_face->setSize(0.5*m_face->getSizeX(), 0., m_face->getSizeZ());
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

			//Bord deja split -> on met seulement a jour les liaisons cells/boundaries
			//----------------------------------------------------------------------
			else {

				//Face in the x-direction
				//-----------------------
        if (abs(m_face->getNormal().getX()) > 0.99) {
          //Left side
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            //First face
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            //Second face
            else {
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
          }
          //Right side
          else {
            //First face
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              m_cellLeft = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            //Second face
            else {
              m_cellLeft = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
          }
        }
        //Face in the y-direction
        //-----------------------
        else if (abs(m_face->getNormal().getY()) > 0.99) {
          //Bottom side
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            //First face
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            //Second face
            else {
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
          }
          //Top side
          else {
            //First face
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              m_cellLeft = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            //Second face
            else {
              m_cellLeft = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
          }
        }
			}
		}
	}
	else {

		//--------------------------------------------------
		//--------------------- Cas 3D ---------------------
		//--------------------------------------------------

		//Bord pas encore split -> Creation des boundaries enfants
		//---------------------------------------------------
		if (m_lvl == cellRef->getLvl()) {

			//Creation des boundaries et faces enfants avec premiere initialization
			//----------------------------------------------------------------
			for (int i = 0; i < 4; i++) {
				this->creerBordChild();
				m_boundariesChildren[i]->m_face = m_face->creerNouvelleFace();
				m_boundariesChildren[i]->m_face->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
        m_boundariesChildren[i]->m_face->setSize(0.5*m_face->getSize());
			}

      //Face in the x-direction
      //-----------------------
      if (abs(m_face->getNormal().getX()) > 0.99) {
        //Left side
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(0));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(2));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(1));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(3));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
            }
          }
        }
        //Right side
        else {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(0));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(2));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
          }
        }
      }
      //Face in the y-direction
      //-----------------------
      else if (abs(m_face->getNormal().getY()) > 0.99) {
        //Bottom side
        if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(0));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(2));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(1));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(3));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
            }

          }
        }
        //Top side
        else {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(0));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(2));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
          }
        }
      }
      //Face in the z-direction
      //-----------------------
      else if (abs(m_face->getNormal().getZ()) > 0.99) {
        //Back side
        if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(0));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(2));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(1));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(m_cellLeft);
              m_boundariesChildren[i]->initializeDroite(cellRef->getCellChild(3));
              m_cellLeft->addBoundary(m_boundariesChildren[i]);
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
            }
          }
        }
        //Front side
        else {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(0));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(0)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(2));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(2)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(1));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(1)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
            }
            //Fourth face
            else {
              m_boundariesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_boundariesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              m_boundariesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(3)->addBoundary(m_boundariesChildren[i]);
              m_cellRight->addBoundary(m_boundariesChildren[i]);
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

		//Bord deja split -> on met seulement a jour les liaisons cells/boundaries
		//----------------------------------------------------------------------
		else {

      //Face in the x-direction
      //-----------------------
      if (abs(m_face->getNormal().getX()) > 0.99) {
        //Left side
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //First face
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            else {
              //Second face
              m_cellRight = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Third face
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellRight = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
          }
        }
        //Right side
        else {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //First face
              m_cellLeft = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            else {
              //Second face
              m_cellLeft = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Third face
              m_cellLeft = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
          }
        }
      }
      //Face in the y-direction
      //-----------------------
      else if (abs(m_face->getNormal().getY()) > 0.99) {
        //Bottom side
        if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //First face
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            else {
              //Second face
              m_cellRight = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Third face
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellRight = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
          }
        }
        //Top side
        else {
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //First face
              m_cellLeft = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            else {
              //Second face
              m_cellLeft = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Third face
              m_cellLeft = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
          }
        }
      }
      //Face in the z-direction
      //-----------------------
      else if (abs(m_face->getNormal().getZ()) > 0.99) {
        //Back side
        if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              //First face
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            else {
              //Second face
              m_cellRight = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              //Third face
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellRight = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
          }
        }
        //Front side
        else {
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              //First face
              m_cellLeft = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addBoundary(this);
            }
            else {
              //Second face
              m_cellLeft = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addBoundary(this);
            }
          }
          else {
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              //Third face
              m_cellLeft = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addBoundary(this);
            }
            else {
              //Fourth face
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addBoundary(this);
            }
          }
        }
      }
		}
	}
}

//***********************************************************************

void CellInterface::deraffineBordExterne(Cell *cellRef)
{
  //On parcourt seulement les boundaries parents pour regarder si la cell voisine de celle de reference a des enfants,
  //si oui (enfants), on ne peut pas deraffiner le bord, on reaffecte donc les liaisons cells/boundaries des boundaries enfants,
  //si non (pas enfants), on peut deraffiner le bord et mettre a jour les neighbouring cells.
  //Plus, si je suis un bord parent, qui a donc des enfants, mais que la cell de reference ne les connait pas encore, on les ajoute a ses boundaries.

  //Parcours les boundaries parents
  if (cellRef->getLvl() == m_lvl) {
    //cellRef est la cell gauche
    if (m_cellLeft == cellRef) {
      //Si la cell voisine (droite) a des enfants, on reaffecte les liaisons cells/boundaries des boundaries enfants
      if (m_cellRight->getSplit()) {
        for (unsigned int bordChild = 0; bordChild < m_boundariesChildren.size(); bordChild++) {
          m_boundariesChildren[bordChild]->initializeGauche(cellRef);
        }
      }
      //La cell voisine (droite) n'a pas d'enfants, on deraffine le bord et met a jour les neighbouring cells
      else {
        //Il faut aussi enlever ses boundaries enfants des boundaries de la cell droite et de mes boundaries
        for (unsigned int bordChild = 0; bordChild < m_boundariesChildren.size(); bordChild++) {
          m_cellRight->deleteBoundary(m_boundariesChildren[bordChild]);
          cellRef->deleteBoundary(m_boundariesChildren[bordChild]);
        }
        this->deraffineBordsChildren();
      }
    }
    //cellRef est la cell droite
    else {
      //Si la cell voisine (gauche) a des enfants, on reaffecte les liaisons cells/boundaries des boundaries enfants
      if (m_cellLeft->getSplit()) {
        for (unsigned int bordChild = 0; bordChild < m_boundariesChildren.size(); bordChild++) {
          m_boundariesChildren[bordChild]->initializeDroite(cellRef);
        }
      }
      //La cell voisine (gauche) n'a pas d'enfants, on deraffine le bord et met a jour les neighbouring cells
      else {
        //Il faut aussi enlever ses boundaries enfants des boundaries de la cell gauche et de mes boundaries
        for (unsigned int bordChild = 0; bordChild < m_boundariesChildren.size(); bordChild++) {
          m_cellLeft->deleteBoundary(m_boundariesChildren[bordChild]);
          cellRef->deleteBoundary(m_boundariesChildren[bordChild]);
        }
        this->deraffineBordsChildren();
      }
    }

    //Si je suis un bord parent, qui a donc des enfants, mais que la cell de reference ne les connait pas encore, on les ajoute a ses boundaries.
    for (unsigned int bordChild = 0; bordChild < m_boundariesChildren.size(); bordChild++) {
      bool ajoutChildAuxBordsCellRef(true);
      for (int i = 0; i < cellRef->getBordsSize(); i++) {
        if (cellRef->getBord(i) == m_boundariesChildren[bordChild]) { ajoutChildAuxBordsCellRef = false; break; }
      }
      if (ajoutChildAuxBordsCellRef) { cellRef->addBoundary(m_boundariesChildren[bordChild]); }
    }
  }
}

//***********************************************************************

void CellInterface::finalizeFace()
{
  delete m_face;
}

//***********************************************************************

void CellInterface::deraffineBordsChildren()
{
  for (unsigned int i = 0; i < m_boundariesChildren.size(); i++) {
    m_boundariesChildren[i]->finalizeFace();
    delete m_boundariesChildren[i];
  }
  m_boundariesChildren.clear();
}

//***********************************************************************

void CellInterface::constructionTableauBordsExternesLvl(vector<CellInterface *> *boundariesLvl)
{
  for (unsigned int i = 0; i < m_boundariesChildren.size(); i++) {
    boundariesLvl[m_lvl + 1].push_back(m_boundariesChildren[i]);
  }
}

//***********************************************************************

bool CellInterface::getSplit() const
{
  bool split = false;
  if (m_boundariesChildren.size() > 0) { split = true; }
  return split;
}

//***********************************************************************

int CellInterface::getLvl() const
{
  return m_lvl;
}

//***********************************************************************

int CellInterface::getNumberBordsChildren() const
{
  return m_boundariesChildren.size();
}

//***********************************************************************

CellInterface* CellInterface::getBordChild(const int &numChild)
{
  return m_boundariesChildren[numChild];
}

//***********************************************************************