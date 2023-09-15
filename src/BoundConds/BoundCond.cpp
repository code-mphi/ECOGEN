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

#include "BoundCond.h"

//***********************************************************************

BoundCond::BoundCond(){}

//***********************************************************************

BoundCond::BoundCond(int numPhysique) : m_numPhysique(numPhysique), m_boundData(std::vector<double>(5,0))
{}

//***********************************************************************

BoundCond::BoundCond(const BoundCond &Source, const int& lvl) : CellInterface(lvl), m_numPhysique(Source.m_numPhysique), m_boundData(std::vector<double>(5,0))
{}

//***********************************************************************

BoundCond::~BoundCond() {}

//***********************************************************************

void BoundCond::initialize(Cell* cellLeft, Cell* /*cellRight*/)
{
  m_cellLeft = cellLeft;
  m_cellRight = NULL;
}

//***********************************************************************

void BoundCond::computeFlux(double& dtMax, Limiter& globalLimiter, Limiter& interfaceLimiter,
  Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter, Prim type)
{
  this->solveRiemann(dtMax, globalLimiter, interfaceLimiter, globalVolumeFractionLimiter, interfaceVolumeFractionLimiter, type);
  this->subtractFlux(1.); //Substract flux on left cell
}

//***********************************************************************

void BoundCond::computeFluxAddPhys(AddPhys &addPhys)
{
  addPhys.computeFluxAddPhysBoundary(this);
}

//***********************************************************************

void BoundCond::solveRiemann(double& dtMax, Limiter& /*globalLimiter*/, Limiter& /*interfaceLimiter*/, Limiter& /*globalVolumeFractionLimiter*/, Limiter& /*interfaceVolumeFractionLimiter*/, Prim type)
{
  bufferCellLeft->copyVec(m_cellLeft->getPhases(type), m_cellLeft->getMixture(type), m_cellLeft->getTransports(type));
  //Computation of extended variables (Phases, Mixture, AddPhys)
  bufferCellLeft->fulfillState();

  //Vector and tensor projections on geometric reference frame of the face
  bufferCellLeft->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());

  //Riemann problem
  double dxLeft(m_cellLeft->getElement()->getLCFL());
  dxLeft = dxLeft*std::pow(2., (double)m_lvl);
  this->solveRiemannBoundary(*bufferCellLeft, dxLeft, dtMax);

  //Handling of transport functions (m_Sm known: need to be called after Riemann solver)
  if (numberTransports > 0) { this->solveRiemannTransportBoundary(*bufferCellLeft); }

  //Flux projection on absolute reference frame
  model->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
}

//***********************************************************************

double BoundCond::getBoundData(VarBoundary var) const
{
  return m_boundData[var];
}

//***************************************************************************
//******************************AMR Method***********************************
//***************************************************************************

void BoundCond::computeFluxXi()
{
  if ((m_cellLeft->getXi() > 0.05)) {
    m_cellLeft->addFluxXi(0.1);
  }
}

//****************************************************************************

void BoundCond::raffineCellInterfaceExterne(const int& nbCellsY, const int& nbCellsZ, const double& dXParent, const double& dYParent,
  const double& dZParent, Cell* cellRef, const int& dim)
{
  //Le cell interface est une CL -> Creation des children cell interfaces
  double surfaceChild(std::pow(0.5, dim - 1.)*m_face->getSurface());
  double epsilon(1.e-6);
  int allocateSlopeLocal = 1;

  if (nbCellsZ == 1) {
    if (nbCellsY == 1) {

      //---------------------------------------------------
      //--------------------- 1D case ---------------------
      //---------------------------------------------------

      this->creerCellInterfaceChild();
      m_cellInterfacesChildren[0]->creerFaceChild(this);
      m_cellInterfacesChildren[0]->getFace()->initializeOthers(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
      m_cellInterfacesChildren[0]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY(), m_face->getPos().getZ());
      m_cellInterfacesChildren[0]->getFace()->setSize(m_face->getSize());
      if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
        //Cell interface number 1 (left)
        m_cellInterfacesChildren[0]->initializeGauche(cellRef->getCellChild(0));
        cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[0]);
      }
      else {
        //Cell interface number 2 (right)
        m_cellInterfacesChildren[0]->initializeGauche(cellRef->getCellChild(1));
        cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[0]);
      }
      m_cellInterfacesChildren[0]->allocateSlopes(allocateSlopeLocal);
    }
    else {

      //---------------------------------------------------
      //--------------------- 2D case ---------------------
      //---------------------------------------------------

      //Creation des cell interfaces et faces enfants avec premiere initialization
      //--------------------------------------------------------------------------
      for (int i = 0; i < 2; i++) {
        this->creerCellInterfaceChild();
        m_cellInterfacesChildren[i]->creerFaceChild(this);
        m_cellInterfacesChildren[i]->getFace()->initializeOthers(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
      }

      //Face selon X
      //------------
      if (std::fabs(m_face->getNormal().getX()) > epsilon) {
        //Cote gauche
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          for (int i = 0; i < 2; i++) {
            //First face
            if (i == 0) {
              m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(0));
              cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Second face
            else {
              m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(2));
              cellRef->getCellChild(2)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            m_cellInterfacesChildren[i]->getFace()->setSize(0., 0.5*m_face->getSizeY(), m_face->getSizeZ());
          }
        }
        //Cote droite
        else {
          for (int i = 0; i < 2; i++) {
            //First face
            if (i == 0) {
              m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(1));
              cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Second face
            else {
              m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              cellRef->getCellChild(3)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            m_cellInterfacesChildren[i]->getFace()->setSize(0., 0.5*m_face->getSizeY(), m_face->getSizeZ());
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
              m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(0));
              cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Second face
            else {
              m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(1));
              cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            m_cellInterfacesChildren[i]->getFace()->setSize(0.5*m_face->getSizeX(), 0., m_face->getSizeZ());
          }
        }
        //Cote haut
        else {
          for (int i = 0; i < 2; i++) {
            //First face
            if (i == 0) {
              m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(2));
              cellRef->getCellChild(2)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Second face
            else {
              m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              cellRef->getCellChild(3)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            m_cellInterfacesChildren[i]->getFace()->setSize(0.5*m_face->getSizeX(), 0., m_face->getSizeZ());
          }
        }
      }

      //Association of model and slopes
      //-------------------------------
      for (int i = 0; i < 2; i++) {
        m_cellInterfacesChildren[i]->allocateSlopes(allocateSlopeLocal);
      }

    }
  }
  else {

    //---------------------------------------------------
    //--------------------- 3D case ---------------------
    //---------------------------------------------------

    //Creation des cell interfaces et faces enfants avec premiere initialization
    //--------------------------------------------------------------------------
    for (int i = 0; i < 4; i++) {
      this->creerCellInterfaceChild();
      m_cellInterfacesChildren[i]->creerFaceChild(this);
      m_cellInterfacesChildren[i]->getFace()->initializeOthers(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
      m_cellInterfacesChildren[i]->getFace()->setSize(0.5*m_face->getSize());
    }

    //Face selon X
    //------------
    if (std::fabs(m_face->getNormal().getX()) > epsilon) {
      //Cote gauche
      if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(0));
            cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(4));
            cellRef->getCellChild(4)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(2));
            cellRef->getCellChild(2)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Fourth face
          else {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(6));
            cellRef->getCellChild(6)->addCellInterface(m_cellInterfacesChildren[i]);
          }
        }
      }
      //Cote droite
      else {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(1));
            cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(5));
            cellRef->getCellChild(5)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(3));
            cellRef->getCellChild(3)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Fourth face
          else {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(7));
            cellRef->getCellChild(7)->addCellInterface(m_cellInterfacesChildren[i]);
          }
        }
      }
    }

    //Face selon Y
    //------------
    else if (std::fabs(m_face->getNormal().getY()) > epsilon) {
      //Cote bas
      if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(0));
            cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(1));
            cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(4));
            cellRef->getCellChild(4)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Fourth face
          else {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(5));
            cellRef->getCellChild(5)->addCellInterface(m_cellInterfacesChildren[i]);
          }
        }
      }
      //Cote haut
      else {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(2));
            cellRef->getCellChild(2)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(3));
            cellRef->getCellChild(3)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(6));
            cellRef->getCellChild(6)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Fourth face
          else {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(7));
            cellRef->getCellChild(7)->addCellInterface(m_cellInterfacesChildren[i]);
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
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(0));
            cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(1));
            cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(2));
            cellRef->getCellChild(2)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Fourth face
          else {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(3));
            cellRef->getCellChild(3)->addCellInterface(m_cellInterfacesChildren[i]);
          }
        }
      }
      //Cote derriere
      else {
        for (int i = 0; i < 4; i++) {
          //First face
          if (i == 0) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(4));
            cellRef->getCellChild(4)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Second face
          else if (i == 1) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(5));
            cellRef->getCellChild(5)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Third face
          else if (i == 2) {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(6));
            cellRef->getCellChild(6)->addCellInterface(m_cellInterfacesChildren[i]);
          }
          //Fourth face
          else {
            m_cellInterfacesChildren[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(7));
            cellRef->getCellChild(7)->addCellInterface(m_cellInterfacesChildren[i]);
          }
        }
      }
    }

    //Association du model et des slopes
    //-----------------------------------
    for (int i = 0; i < 4; i++) {
      m_cellInterfacesChildren[i]->allocateSlopes(allocateSlopeLocal);
    }

  }
}

//***********************************************************************

void BoundCond::deraffineCellInterfaceExterne(Cell* cellRef)
{
  //Dans le cas des CL, on parcourt toujours les parent cell interfaces mais on peut directement les deraffiner et mettre a jour la cell gauche.

  //Parcours les parent cell interfaces
  if (cellRef->getLvl() == m_lvl) {
    //cellRef est forcement la cell gauche, on deraffine le parent cell interface et on met a jour la cell gauche (de reference)
    for (unsigned int cellInterfaceChild = 0; cellInterfaceChild < m_cellInterfacesChildren.size(); cellInterfaceChild++) {
      cellRef->deleteCellInterface(m_cellInterfacesChildren[cellInterfaceChild]);
    }
    this->deraffineCellInterfacesChildren();
  }
}

//****************************************************************************