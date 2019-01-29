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

//! \file      MeshCartesian.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot
//! \version   1.0
//! \date      September 06 2018

#include <iostream>
#include <algorithm>
#include <sstream>

#include "MeshCartesian.h"

using namespace std;

//***********************************************************************

MeshCartesian::MeshCartesian(double lX, int numberCellsX, double lY, int numberCellsY, double lZ, int numberCellsZ,
  std::vector<stretchZone> stretchX, std::vector<stretchZone> stretchY, std::vector<stretchZone> stretchZ) :
  m_lX(lX), m_numberCellsXGlobal(numberCellsX),
  m_lY(lY), m_numberCellsYGlobal(numberCellsY),
  m_lZ(lZ), m_numberCellsZGlobal(numberCellsZ)
{
  for (unsigned int i = 0; i < stretchX.size(); i++) { m_stretchX.push_back(stretchX[i]); }
  for (unsigned int i = 0; i < stretchY.size(); i++) { m_stretchY.push_back(stretchY[i]); }
  for (unsigned int i = 0; i < stretchZ.size(); i++) { m_stretchZ.push_back(stretchZ[i]); }

  m_numberCellsCalcul = m_numberCellsXGlobal*m_numberCellsYGlobal*m_numberCellsZGlobal;
  m_geometrie = 0;
  if (numberCellsX != 1) { m_geometrie += 1; }
  if (numberCellsY != 1) { m_geometrie += 1; }
  if (numberCellsZ != 1) { m_geometrie += 1; }
  m_type = REC;
  m_numberCpuX = 1;
  m_numberCpuY = 1;
  m_numberCpuZ = 1;
  m_CpuCoordX = 0;
  m_CpuCoordY = 0;
  m_CpuCoordZ = 0;
  m_offsetX = 0;
  m_offsetY = 0;
  m_offsetZ = 0;
}

//***********************************************************************

MeshCartesian::~MeshCartesian(){
  delete[] m_elements;
  delete[] m_faces;
  for (int j = m_numberBoundCondInit; j < 6; j++) {
    switch (j) {
    case 0:
      delete m_limXm; break;
    case 1:
      delete m_limXp; break;
    case 2:
      delete m_limYm; break;
    case 3:
      delete m_limYp; break;
    case 4:
      delete m_limZm; break;
    case 5:
      delete m_limZp; break;
    default:
      Errors::errorMessage("Probleme de limites dans delete(MeshCartesian)"); break;
    }
  }
}

//***********************************************************************

void MeshCartesian::attributLimites(std::vector<BoundCond*> &boundCond)
{
  m_numberBoundCondInit = boundCond.size();
  for (int i = 0; i < m_numberBoundCondInit; i++) {
    switch (i) {
    case 0:
      m_limXm = boundCond[i]; break;
    case 1:
      m_limXp = boundCond[i]; break;
    case 2:
      m_limYm = boundCond[i]; break;
    case 3:
      m_limYp = boundCond[i]; break;
    case 4:
      m_limZm = boundCond[i]; break;
    case 5:
      m_limZp = boundCond[i]; break;
    default:
      Errors::errorMessage("Probleme de limites dans attributLimites"); break;
    }
  }

  for (int j = m_numberBoundCondInit; j < 6; j++) {
    switch (j) {
    case 0:
      m_limXm = new BoundCondAbs(j+1); break;
    case 1:
      m_limXp = new BoundCondAbs(j+1); break;
    case 2:
      m_limYm = new BoundCondAbs(j+1); break;
    case 3:
      m_limYp = new BoundCondAbs(j+1); break;
    case 4:
      m_limZm = new BoundCondAbs(j+1); break;
    case 5:
      m_limZp = new BoundCondAbs(j+1); break;
    default:
      Errors::errorMessage("Probleme de limites dans attributLimites"); break;
    }
  }
}

//***********************************************************************

void MeshCartesian::recupereIJK(const int &index, int &i, int &j, int &k) const
{
  int reste;
  k = index / (m_numberCellsX*m_numberCellsY);
  reste = index % (m_numberCellsX*m_numberCellsY);
  j = reste / m_numberCellsX;
  i = reste % m_numberCellsX;
}

//***********************************************************************

void MeshCartesian::construitIGlobal(const int &i, const int &j, const int &k, int &index) const
{
  index = 0;
  if (m_numberCellsX != 1) index += i;
  if (m_numberCellsY != 1) index += j*m_numberCellsX;
  if (m_numberCellsZ != 1) index += k*m_numberCellsX*m_numberCellsY;
}

//***********************************************************************

int MeshCartesian::initializeGeometrie(Cell ***cells, CellInterface ***bord, bool pretraitementParallele, string ordreCalcul)
{
  this->meshStretching();
  if (Ncpu == 1)
  {
    this->initializeGeometrieMonoCpu(cells, bord, ordreCalcul);
  }
  else
  {
    this->initializeGeometrieParallele(cells, bord, ordreCalcul);
  }
  return m_geometrie;
}

//***********************************************************************

void MeshCartesian::meshStretching()
{
  //initial value for cell sizes
  m_numberCellsCalcul = m_numberCellsXGlobal*m_numberCellsYGlobal*m_numberCellsZGlobal;
  double dX = m_lX / static_cast<double>(m_numberCellsXGlobal);
  double dY = m_lY / static_cast<double>(m_numberCellsYGlobal);
  double dZ = m_lZ / static_cast<double>(m_numberCellsZGlobal);

  //Stretching on X
  if (m_stretchX.size() == 0) {
    for(int i=0; i<m_numberCellsXGlobal; i++){ m_dXi.push_back(dX); m_posXi.push_back((i + 0.5)*dX);}
  }
  else {
    m_numberCellsXGlobal = 0;
    for (unsigned int z = 0; z < m_stretchX.size(); z++) { 
      m_numberCellsXGlobal += m_stretchX[z].stretching(m_dXi, m_posXi); 
    }
  }
  //Stretching on Y
  if (m_stretchY.size() == 0) {
    for (int i = 0; i<m_numberCellsYGlobal; i++) { m_dYj.push_back(dY); m_posYj.push_back((i+0.5)*dY); }
  }
  else {
    m_numberCellsYGlobal = 0;
    for (unsigned int z = 0; z < m_stretchY.size(); z++) { 
      m_numberCellsYGlobal += m_stretchY[z].stretching(m_dYj, m_posYj); 
    }
  }
  //Stretching on Z
  if (m_stretchZ.size() == 0) {
    for (int i = 0; i<m_numberCellsZGlobal; i++) { m_dZk.push_back(dZ); m_posZk.push_back((i + 0.5)*dZ); }
  }
  else {
    m_numberCellsZGlobal = 0;
    for (unsigned int z = 0; z < m_stretchZ.size(); z++) {
      m_numberCellsZGlobal += m_stretchZ[z].stretching(m_dZk, m_posZk);
    }
  }
  //Updating cells number
  m_numberCellsCalcul = m_numberCellsXGlobal*m_numberCellsYGlobal*m_numberCellsZGlobal;
}

//***********************************************************************

void MeshCartesian::initializeGeometrieMonoCpu(Cell ***cells, CellInterface ***bord, string ordreCalcul)
{
  int ix, iy, iz;
  int compteMaillesParallele(0);

  m_numberCellsX = m_numberCellsXGlobal;
  m_numberCellsY = m_numberCellsYGlobal;
  m_numberCellsZ = m_numberCellsZGlobal;
  m_numberCellsTotal = m_numberCellsCalcul;

  //Declaration du tableau de mailles et d'elements
  //-----------------------------------------------
  m_elements = new ElementCartesian[m_numberCellsCalcul];
  (*cells) = new Cell*[m_numberCellsCalcul];
  for (int i = 0; i < m_numberCellsCalcul; i++) { 
    if(ordreCalcul == "FIRSTORDER") { (*cells)[i] = new Cell; }
    else { (*cells)[i] = new CellO2; }
    (*cells)[i]->setElement(&m_elements[i], i); 
  }

  //Attribution des donnees geometriques cartesiennes aux mailles
  //-------------------------------------------------------------
  Coord tangent, normal, binormal;
  double surface(1.), volume(0.);
  double posX, posY, posZ;
  for (int i = 0; i < m_numberCellsCalcul; i++)
  {
    this->recupereIJK(i, ix, iy, iz);
    volume = m_dXi[ix] * m_dYj[iy] * m_dZk[iz];
    (*cells)[i]->getElement()->setVolume(volume);

    //CFL lenght
    double lCFL(1.e10);
    if (m_numberCellsX != 1) { lCFL = min(lCFL, m_dXi[ix]); }
    if (m_numberCellsY != 1) { lCFL = min(lCFL, m_dYj[iy]); }
    if (m_numberCellsZ != 1) { lCFL = min(lCFL, m_dZk[iz]); }
    if (m_geometrie > 1) lCFL *= 0.6;

    (*cells)[i]->getElement()->setLCFL(lCFL);
    (*cells)[i]->getElement()->setPos(m_posXi[ix], m_posYj[iy], m_posZk[iz]);
    (*cells)[i]->getElement()->setSize(m_dXi[ix], m_dYj[iy], m_dZk[iz]);
  }

  //Attribution des donnees geometriques cartesiennes aux boundaries de mailles
  //----------------------------------------------------------------------
  //Determination du number de faces selon compute 1D/2D/3D
  //m_numberFaces : number total de faces
  //m_numberFacesLimites : number de faces sur les limites
  //numberFacesInternes : number de faces communes a deux mailles
  m_numberFacesTotal = 0;
  int m_numberFacesLimites(0);
  if (m_numberCellsX != 1)
  {
    m_numberFacesTotal += (m_numberCellsX + 1)*m_numberCellsY*m_numberCellsZ;
    m_numberFacesLimites += 2 * m_numberCellsY*m_numberCellsZ;
  }
  if (m_numberCellsY != 1)
  {
    m_numberFacesTotal += (m_numberCellsY + 1)*m_numberCellsX*m_numberCellsZ;
    m_numberFacesLimites += 2 * m_numberCellsX*m_numberCellsZ;
  }
  if (m_numberCellsZ != 1)
  {
    m_numberFacesTotal += (m_numberCellsZ + 1)*m_numberCellsX*m_numberCellsY;
    m_numberFacesLimites += 2 * m_numberCellsX*m_numberCellsY;
  }
  int numberFacesInternes(m_numberFacesTotal - m_numberFacesLimites);
  //Allocation globale
  (*bord) = new CellInterface*[m_numberFacesTotal];
  m_faces = new FaceCartesian[m_numberFacesTotal];

  //Initialization des faces internes
  //*********************************
  int iMailleG, iMailleD, iFace(0), iTemp;
  Cell *BGM=0, *BGP=0, *BDM=0, *BDP=0;
  //Faces selon X
  tangent.setXYZ(0., 1., 0.); normal.setXYZ(1., 0., 0.); binormal.setXYZ(0., 0., 1.);
  for (ix = 0; ix < m_numberCellsX - 1; ix++)
  {
    for (iy = 0; iy < m_numberCellsY; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        if (ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
        else { (*bord)[iFace] = new CellInterfaceO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix + 1, iy, iz, iMailleD);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        surface = m_dYj[iy] * m_dZk[iz];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(0., m_dYj[iy], m_dZk[iz]);
        posX = m_posXi[ix] + 0.5*m_dXi[ix];
        posY = m_posYj[iy];
        posZ = m_posZk[iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        //Preparation des cells pour ordre 2 multislopes
        BGM = 0; BGP = 0; BDM = 0; BDP = 0;
        if (ix != 0) { this->construitIGlobal(ix - 1, iy, iz, iTemp); BGM = (*cells)[iTemp]; }
        if (ix < m_numberCellsX - 1) { this->construitIGlobal(ix + 1, iy, iz, iTemp); BGP = (*cells)[iTemp]; }
        if (ix < m_numberCellsX - 2) { this->construitIGlobal(ix + 2, iy, iz, iTemp); BDM = (*cells)[iTemp]; }
        this->construitIGlobal(ix, iy, iz, iTemp); BDP = (*cells)[iTemp];
        iFace++;
      }
    }
  }
  //Faces selon Y
  tangent.setXYZ(-1., 0., 0.); normal.setXYZ(0., 1., 0.); binormal.setXYZ(0., 0., 1.);
  for (ix = 0; ix < m_numberCellsX; ix++)
  {
    for (iy = 0; iy < m_numberCellsY - 1; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        if(ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
        else { (*bord)[iFace] = new CellInterfaceO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix, iy + 1, iz, iMailleD);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        surface = m_dXi[ix] * m_dZk[iz];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[ix], 0., m_dZk[iz]);
        posX = m_posXi[ix];
        posY = m_posYj[iy] + 0.5*m_dYj[iy];
        posZ = m_posZk[iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        //Preparation des cells pour ordre 2 multislopes
        BGM = 0; BGP = 0; BDM = 0; BDP = 0;
        if (iy != 0) { this->construitIGlobal(ix, iy - 1, iz, iTemp); BGM = (*cells)[iTemp]; }
        if (iy < m_numberCellsY - 1) { this->construitIGlobal(ix, iy + 1, iz, iTemp); BGP = (*cells)[iTemp]; }
        if (iy < m_numberCellsY - 2) { this->construitIGlobal(ix, iy + 2, iz, iTemp); BDM = (*cells)[iTemp]; }
        this->construitIGlobal(ix, iy, iz, iTemp); BDP = (*cells)[iTemp];
        iFace++;
      }
    }
  }
  //Faces selon Z
  tangent.setXYZ(1., 0., 0.); normal.setXYZ(0., 0., 1.); binormal.setXYZ(0., 1., 0.);
  for (ix = 0; ix < m_numberCellsX; ix++)
  {
    for (iy = 0; iy < m_numberCellsY; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ - 1; iz++)
      {
        if(ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
        else { (*bord)[iFace] = new CellInterfaceO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix, iy, iz + 1, iMailleD);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        surface = m_dXi[ix] * m_dYj[iy];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[ix], m_dYj[iy], 0.);
        posX = m_posXi[ix];
        posY = m_posYj[iy];
        posZ = m_posZk[iz] + 0.5*m_dZk[iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        //Preparation des cells pour ordre 2 multislopes
        BGM = 0; BGP = 0; BDM = 0; BDP = 0;
        if (iz != 0) { this->construitIGlobal(ix, iy, iz - 1, iTemp); BGM = (*cells)[iTemp]; }
        if (iz < m_numberCellsZ - 1) { this->construitIGlobal(ix, iy, iz + 1, iTemp); BGP = (*cells)[iTemp]; }
        if (iz < m_numberCellsZ - 2) { this->construitIGlobal(ix, iy, iz + 2, iTemp); BDM = (*cells)[iTemp]; }
        this->construitIGlobal(ix, iy, iz, iTemp); BDP = (*cells)[iTemp];
        iFace++;
      }
    }
  }

  //Initialization des faces limites en X
  //---------------------------------------
  if (m_numberCellsX != 1)
  {
    //Limite X=0
    ix = 0;
    tangent.setXYZ(0., -1., 0.); normal.setXYZ(-1., 0., 0.); binormal.setXYZ(0., 0., 1.);
    for (iy = 0; iy < m_numberCellsY; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        m_limXm->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        surface = m_dYj[iy] * m_dZk[iz];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(0., m_dYj[iy], m_dZk[iz]);
        posX = 0.;
        posY = m_posYj[iy];
        posZ = m_posZk[iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Limite X=lX
    ix = m_numberCellsX - 1;
    tangent.setXYZ(0., 1., 0.); normal.setXYZ(1., 0., 0.); binormal.setXYZ(0., 0., 1.);
    for (iy = 0; iy < m_numberCellsY; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        m_limXp->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        surface = m_dYj[iy] * m_dZk[iz];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(0., m_dYj[iy], m_dZk[iz]);
        posX = m_posXi[ix] + 0.5*m_dXi[ix];
        posY = m_posYj[iy];
        posZ = m_posZk[iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
  //Initialization des faces limites en Y
  //---------------------------------------
  if (m_numberCellsY != 1)
  {
    //Limite Y=0
    iy = 0;
    tangent.setXYZ(1., 0., 0.); normal.setXYZ(0., -1., 0.); binormal.setXYZ(0., 0., 1.);
    for (ix = 0; ix < m_numberCellsX; ix++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        m_limYm->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        surface = m_dXi[ix] * m_dZk[iz];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[ix], 0., m_dZk[iz]);
        posX = m_posXi[ix];
        posY = 0.;
        posZ = m_posZk[iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Limite Y=lY
    iy = m_numberCellsY - 1;
    tangent.setXYZ(-1., 0., 0.); normal.setXYZ(0., 1., 0.); binormal.setXYZ(0., 0., 1.);
    for (ix = 0; ix < m_numberCellsX; ix++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        m_limYp->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        surface = m_dXi[ix] * m_dZk[iz];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[ix], 0., m_dZk[iz]);
        posX = m_posXi[ix];
        posY = m_posYj[iy] + 0.5*m_dYj[iy];
        posZ = m_posZk[iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
  //Initialization des faces limites en Z
  //-------------------------------------
  if (m_numberCellsZ != 1)
  {
    //Limite Z=0
    iz = 0;
    tangent.setXYZ(-1., 0., 0.); normal.setXYZ(0., 0., -1.); binormal.setXYZ(0., 1., 0.);
    for (ix = 0; ix < m_numberCellsX; ix++)
    {
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        m_limZm->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        surface = m_dXi[ix] * m_dYj[iy];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[ix], m_dYj[iy], 0.);
        posX = m_posXi[ix];
        posY = m_posYj[iy];
        posZ = 0.;
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Limite Z=lZ
    iz = m_numberCellsZ - 1;
    tangent.setXYZ(1., 0., 0.); normal.setXYZ(0., 0., 1.); binormal.setXYZ(0., 1., 0.);
    for (ix = 0; ix < m_numberCellsX; ix++)
    {
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        m_limZp->creeLimite(&(*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        iMailleD = iMailleG;
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        surface = m_dXi[ix] * m_dYj[iy];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[ix], m_dYj[iy], 0.);
        posX = m_posXi[ix];
        posY = m_posYj[iy];
        posZ = m_posZk[iz] + 0.5*m_dZk[iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
}

//***********************************************************************

void MeshCartesian::initializeGeometrieParallele(Cell ***cells, CellInterface ***bord, string ordreCalcul)
{
  int ix, iy, iz;
  int compteMaillesParallele(0);

  m_numberCellsX = m_numberCellsXGlobal;
  m_numberCellsY = m_numberCellsYGlobal;
  m_numberCellsZ = m_numberCellsZGlobal;
  
  //Declaration du tableau de mailles
  //---------------------------------
  //Decoupage parallele du domain geometrique
  this->decoupageParallele();
  m_elements = new ElementCartesian[m_numberCellsTotal];
  (*cells) = new Cell*[m_numberCellsTotal];
  for (int i = 0; i < m_numberCellsCalcul; i++) {
	  if (ordreCalcul == "FIRSTORDER") { (*cells)[i] = new Cell; }
	  else { (*cells)[i] = new CellO2; }
    (*cells)[i]->setElement(&m_elements[i], i);
  }
	for (int i = m_numberCellsCalcul; i < m_numberCellsTotal; i++) {
		if (ordreCalcul == "FIRSTORDER") { (*cells)[i] = new Cell; }
		else { (*cells)[i] = new CellO2Ghost; }
		(*cells)[i]->setElement(&m_elements[i], i);
	}
  
  //Geometrical data settings for computational cells
  //-------------------------------------------------
  Coord tangent, normal, binormal;
  double surface(1.), volume(0.);
  double posX, posY, posZ;
  for (int i = 0; i < m_numberCellsCalcul; i++)
  {
    this->recupereIJK(i, ix, iy, iz);
    volume = m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz];
    (*cells)[i]->getElement()->setVolume(volume);

    //CFL lenght
    double lCFL(1.e10);
    if (m_numberCellsX != 1) { lCFL = min(lCFL, m_dXi[ix]); }
    if (m_numberCellsY != 1) { lCFL = min(lCFL, m_dYj[iy]); }
    if (m_numberCellsZ != 1) { lCFL = min(lCFL, m_dZk[iz]); }
    if (m_geometrie > 1) lCFL *= 0.6;

    (*cells)[i]->getElement()->setLCFL(lCFL);
    (*cells)[i]->getElement()->setPos(m_posXi[m_offsetX + ix], m_posYj[m_offsetY + iy], m_posZk[m_offsetZ + iz]);
    (*cells)[i]->getElement()->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
  }

  //Attribution des donnees geometriques cartesiennes aux boundaries de mailles
  //----------------------------------------------------------------------
  //Determination du number de faces selon compute 1D/2D/3D
  //m_numberFaces : number total de faces
  //m_numberFacesLimites : number de faces sur les limites
  //numberFacesInternes : number de faces communes a deux mailles
  m_numberFacesTotal = 0;
  int m_numberFacesLimites(0);
  if (m_numberCellsX != 1)
  {
    m_numberFacesTotal += (m_numberCellsX + 1)*m_numberCellsY*m_numberCellsZ;
    m_numberFacesLimites += 2 * m_numberCellsY*m_numberCellsZ;
  }
  if (m_numberCellsY != 1)
  {
    m_numberFacesTotal += (m_numberCellsY + 1)*m_numberCellsX*m_numberCellsZ;
    m_numberFacesLimites += 2 * m_numberCellsX*m_numberCellsZ;
  }
  if (m_numberCellsZ != 1)
  {
    m_numberFacesTotal += (m_numberCellsZ + 1)*m_numberCellsX*m_numberCellsY;
    m_numberFacesLimites += 2 * m_numberCellsX*m_numberCellsY;
  }
  int numberFacesInternes(m_numberFacesTotal - m_numberFacesLimites);
  //Allocation globale
  (*bord) = new CellInterface*[m_numberFacesTotal];
  m_faces = new FaceCartesian[m_numberFacesTotal];
  
  //Initialization des faces internes
  //*********************************
  int iMailleG, iMailleD, iFace(0);
  //Faces selon X
  tangent.setXYZ(0., 1., 0.); normal.setXYZ(1., 0., 0.); binormal.setXYZ(0., 0., 1.);
  for (ix = 0; ix < m_numberCellsX - 1; ix++)
  {
    for (iy = 0; iy < m_numberCellsY; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        if(ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
        else { (*bord)[iFace] = new CellInterfaceO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix + 1, iy, iz, iMailleD);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        surface = m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(0., m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix] + 0.5*m_dXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
  //Faces selon Y
  tangent.setXYZ(-1., 0., 0.); normal.setXYZ(0., 1., 0.); binormal.setXYZ(0., 0., 1.);
  for (ix = 0; ix < m_numberCellsX; ix++)
  {
    for (iy = 0; iy < m_numberCellsY - 1; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        if(ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
        else { (*bord)[iFace] = new CellInterfaceO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix, iy + 1, iz, iMailleD);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        surface = m_dXi[m_offsetX + ix] * m_dZk[m_offsetZ + iz];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[m_offsetX + ix], 0., m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy] + 0.5*m_dYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
  //Faces selon Z
  tangent.setXYZ(1., 0., 0.); normal.setXYZ(0., 0., 1.); binormal.setXYZ(0., 1., 0.);
  for (ix = 0; ix < m_numberCellsX; ix++)
  {
    for (iy = 0; iy < m_numberCellsY; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ - 1; iz++)
      {
        if(ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
        else { (*bord)[iFace] = new CellInterfaceO2; }
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iMailleG);
        this->construitIGlobal(ix, iy, iz + 1, iMailleD);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        surface = m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy];
        m_faces[iFace].initializeAutres(surface, normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], 0.);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz] + 0.5*m_dZk[m_offsetZ + iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }

  compteMaillesParallele = m_numberCellsCalcul;

  //Ghosts cells and boundaries initialisation in X direction
  //*********************************************************
  if (m_numberCellsX != 1)
  {
    //X=0 boundary
    ix = 0;
    binormal.setXYZ(0., 0., 1.);
    for (iy = 0; iy < m_numberCellsY; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        //A) CPU neighbour limits treatment
        //---------------------------------
        if (m_numberCpuX > 1 && m_CpuCoordX > 0) {
          tangent.setXYZ(0., 1., 0.); normal.setXYZ(1., 0., 0.); //Inversion for neighbour CPU matching
          //right and Left cells catching
					this->construitIGlobal(ix, iy, iz, iMailleD);
          iMailleG = compteMaillesParallele++; //Ghost cell number taken in order
          //setting ghost cell geometry
					(*cells)[iMailleG]->getElement()->setPos((*cells)[iMailleD]->getElement()->getPosition());
          (*cells)[iMailleG]->getElement()->setPosX(m_posXi[m_offsetX - 1]);
          (*cells)[iMailleG]->getElement()->setSize(m_dXi[m_offsetX -1], m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
          (*cells)[iMailleG]->getElement()->setVolume(m_dXi[m_offsetX - 1] * m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = min(lCFL, m_dXi[m_offsetX - 1]); }
          if (m_numberCellsY != 1) { lCFL = min(lCFL, m_dYj[m_offsetY + iy]); }
          if (m_numberCellsZ != 1) { lCFL = min(lCFL, m_dZk[m_offsetZ + iz]); }
          if (m_geometrie > 1) lCFL *= 0.6;
          (*cells)[iMailleG]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
          else { (*bord)[iFace] = new CellInterfaceO2; }     
					(*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        }
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          tangent.setXYZ(0., -1., 0.); normal.setXYZ(-1., 0., 0.); 
          //right and Left cells equals
					this->construitIGlobal(ix, iy, iz, iMailleG);
          iMailleD = iMailleG;
          //setting boundary
          m_limXm->creeLimite(&(*bord)[iFace]);
        }
        //Common settings
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        m_faces[iFace].initializeAutres(m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz], normal, tangent, binormal);
        m_faces[iFace].setSize(0., m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix] - 0.5*m_dXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //X=lX Boundary
    ix = m_numberCellsX - 1;
    tangent.setXYZ(0., 1., 0.); normal.setXYZ(1., 0., 0.); binormal.setXYZ(0., 0., 1.);
    for (iy = 0; iy < m_numberCellsY; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      { 
        this->construitIGlobal(ix, iy, iz, iMailleG);
        //A) CPU neighbour limits treatment
        //---------------------------------
        if (m_numberCpuX > 1 && m_CpuCoordX < m_numberCpuX - 1) {
          iMailleD = compteMaillesParallele++; //Ghost cell number taken in order
          //setting ghost cell geometry
          (*cells)[iMailleD]->getElement()->setPos((*cells)[iMailleG]->getElement()->getPosition());
          (*cells)[iMailleD]->getElement()->setPosX(m_posXi[m_offsetX + ix + 1]);
          (*cells)[iMailleD]->getElement()->setSize(m_dXi[m_offsetX + ix + 1], m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
          (*cells)[iMailleD]->getElement()->setVolume(m_dXi[m_offsetX + ix + 1] * m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = min(lCFL, m_dXi[m_offsetX + ix + 1]); }
          if (m_numberCellsY != 1) { lCFL = min(lCFL, m_dYj[m_offsetY + iy]); }
          if (m_numberCellsZ != 1) { lCFL = min(lCFL, m_dZk[m_offsetZ + iz]); }
          if (m_geometrie > 1) lCFL *= 0.6;
          (*cells)[iMailleD]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
          else { (*bord)[iFace] = new CellInterfaceO2; }
          (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        }
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          m_limXp->creeLimite(&(*bord)[iFace]);
          iMailleD = iMailleG;
        }
        //Common settings
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        m_faces[iFace].initializeAutres(m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz], normal, tangent, binormal);
        m_faces[iFace].setSize(0., m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix] + 0.5*m_dXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }

  //Ghosts cells and boundaries initialisation in Y direction
  //*********************************************************
  if (m_numberCellsY != 1)
  {
    //Y=0 boundary
    iy = 0;
    binormal.setXYZ(0., 0., 1.);
    for (ix = 0; ix < m_numberCellsX; ix++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        //A) CPU neighbour limits treatment
        //---------------------------------
        if (m_numberCpuY > 1 && m_CpuCoordY > 0) {
          tangent.setXYZ(-1., 0., 0.); normal.setXYZ(0., 1., 0.);
          //right and Left cells catching
					this->construitIGlobal(ix, iy, iz, iMailleD);
          iMailleG = compteMaillesParallele++; //Ghost cell number taken in order
          //setting ghost cell geometry
          (*cells)[iMailleG]->getElement()->setPos((*cells)[iMailleD]->getElement()->getPosition());
          (*cells)[iMailleG]->getElement()->setPosY(m_posYj[m_offsetY - 1]);
          (*cells)[iMailleG]->getElement()->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY - 1], m_dZk[m_offsetZ + iz]);
          (*cells)[iMailleG]->getElement()->setVolume(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY - 1] * m_dZk[m_offsetZ + iz]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = min(lCFL, m_dXi[m_offsetX + ix]); }
          if (m_numberCellsY != 1) { lCFL = min(lCFL, m_dYj[m_offsetY - 1]); }
          if (m_numberCellsZ != 1) { lCFL = min(lCFL, m_dZk[m_offsetZ + iz]); }
          if (m_geometrie > 1) lCFL *= 0.6;
          (*cells)[iMailleG]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
          else { (*bord)[iFace] = new CellInterfaceO2; }
          (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        } 
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          tangent.setXYZ(1., 0., 0.); normal.setXYZ(0., -1., 0.);
          //right and Left cells equals
					this->construitIGlobal(ix, iy, iz, iMailleG);
          iMailleD = iMailleG;
          //setting boundary
          m_limYm->creeLimite(&(*bord)[iFace]);
        }
        //Common settings
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        m_faces[iFace].initializeAutres(m_dXi[m_offsetX + ix] * m_dZk[m_offsetZ + iz], normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[m_offsetX + ix], 0., m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy] - 0.5*m_dYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Y=lY Boudary
    iy = m_numberCellsY - 1;
    tangent.setXYZ(-1., 0., 0.); normal.setXYZ(0., 1., 0.); binormal.setXYZ(0., 0., 1.);
    for (ix = 0; ix < m_numberCellsX; ix++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        this->construitIGlobal(ix, iy, iz, iMailleG);
        //A) CPU neighbour limits treatment
        //---------------------------------
        if (m_numberCpuY > 1 && m_CpuCoordY < m_numberCpuY - 1) {
          iMailleD = compteMaillesParallele++; //Ghost cell number taken in order
          //setting ghost cell geometry
          (*cells)[iMailleD]->getElement()->setPos((*cells)[iMailleG]->getElement()->getPosition());
          (*cells)[iMailleD]->getElement()->setPosY(m_posYj[m_offsetY + iy + 1]);
          (*cells)[iMailleD]->getElement()->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy + 1], m_dZk[m_offsetZ + iz]);
          (*cells)[iMailleD]->getElement()->setVolume(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy + 1] * m_dZk[m_offsetZ + iz]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = min(lCFL, m_dXi[m_offsetX + ix]); }
          if (m_numberCellsY != 1) { lCFL = min(lCFL, m_dYj[m_offsetY + iy + 1]); }
          if (m_numberCellsZ != 1) { lCFL = min(lCFL, m_dZk[m_offsetZ + iz]); }
          if (m_geometrie > 1) lCFL *= 0.6;
          (*cells)[iMailleD]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
          else { (*bord)[iFace] = new CellInterfaceO2; }
          (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        }
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          m_limYp->creeLimite(&(*bord)[iFace]);
          iMailleD = iMailleG;
        }
        //Common settings
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        m_faces[iFace].initializeAutres(m_dXi[m_offsetX + ix] * m_dZk[m_offsetZ + iz], normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[m_offsetX + ix], 0., m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy] + 0.5*m_dYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }

  //Ghosts cells and boundaries initialisation in Z direction
  //*********************************************************
  if (m_numberCellsZ != 1)
  {
    //Z=0 Boundary
    iz = 0;
    binormal.setXYZ(0., 1., 0.);
    for (ix = 0; ix < m_numberCellsX; ix++)
    {
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        //A) CPU neighbour limits treatment
        //---------------------------------
        if (m_numberCpuZ > 1 && m_CpuCoordZ > 0) {
          tangent.setXYZ(1., 0., 0.); normal.setXYZ(0., 0., 1.);  //Inversion for neighbour CPU matching
          //right and Left cells catching
          this->construitIGlobal(ix, iy, iz, iMailleD);
          iMailleG = compteMaillesParallele++; //Ghost cell number taken in order
          //setting ghost cell geometry
          (*cells)[iMailleG]->getElement()->setPos((*cells)[iMailleD]->getElement()->getPosition());
          (*cells)[iMailleG]->getElement()->setPosZ(m_posZk[m_offsetZ - 1]);
          (*cells)[iMailleG]->getElement()->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], m_dZk[m_offsetZ - 1]);
          (*cells)[iMailleG]->getElement()->setVolume(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ - 1]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = min(lCFL, m_dXi[m_offsetX + ix]); }
          if (m_numberCellsY != 1) { lCFL = min(lCFL, m_dYj[m_offsetY + iy]); }
          if (m_numberCellsZ != 1) { lCFL = min(lCFL, m_dZk[m_offsetZ - 1]); }
          if (m_geometrie > 1) lCFL *= 0.6;
          (*cells)[iMailleG]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
          else { (*bord)[iFace] = new CellInterfaceO2; }
          (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        }
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          tangent.setXYZ(-1., 0., 0.); normal.setXYZ(0., 0., -1.);
          //right and Left cells equals
          this->construitIGlobal(ix, iy, iz, iMailleG);
          iMailleD = iMailleG;
          //setting boundary
          m_limZm->creeLimite(&(*bord)[iFace]);
        }
        //Common settings
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        m_faces[iFace].initializeAutres(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy], normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], 0.);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz] - 0.5*m_dZk[m_offsetZ + iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
    //Z=lZ boundary
    iz = m_numberCellsZ - 1;
    tangent.setXYZ(1., 0., 0.); normal.setXYZ(0., 0., 1.); binormal.setXYZ(0., 1., 0.);
    for (ix = 0; ix < m_numberCellsX; ix++)
    {
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        this->construitIGlobal(ix, iy, iz, iMailleG);
        //A) CPU neighbour limits treatment
        //---------------------------------
        if (m_numberCpuZ > 1 && m_CpuCoordZ < m_numberCpuZ - 1) {
          iMailleD = compteMaillesParallele++; //Ghost cell number taken in order
          //setting ghost cell geometry
          (*cells)[iMailleD]->getElement()->setPos((*cells)[iMailleG]->getElement()->getPosition());
          (*cells)[iMailleD]->getElement()->setPosZ(m_posZk[m_offsetZ + iz + 1]);
          (*cells)[iMailleD]->getElement()->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz + 1]);
          (*cells)[iMailleD]->getElement()->setVolume(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz + 1]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = min(lCFL, m_dXi[m_offsetX + ix]); }
          if (m_numberCellsY != 1) { lCFL = min(lCFL, m_dYj[m_offsetY + iy]); }
          if (m_numberCellsZ != 1) { lCFL = min(lCFL, m_dZk[m_offsetZ + iz + 1]); }
          if (m_geometrie > 1) lCFL *= 0.6;
          (*cells)[iMailleD]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { (*bord)[iFace] = new CellInterface; }
          else { (*bord)[iFace] = new CellInterfaceO2; }
          (*cells)[iMailleD]->addBoundary((*bord)[iFace]);
        }
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          m_limZp->creeLimite(&(*bord)[iFace]);
          iMailleD = iMailleG;
        }
        //Common settings
        (*cells)[iMailleG]->addBoundary((*bord)[iFace]);
        (*bord)[iFace]->setFace(&m_faces[iFace]);
        (*bord)[iFace]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
        m_faces[iFace].initializeAutres(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy], normal, tangent, binormal);
        m_faces[iFace].setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], 0.);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz] + 0.5*m_dZk[m_offsetZ + iz];
        m_faces[iFace].setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
}

//***********************************************************************

void MeshCartesian::decoupageParallele()
{
  int ix, iy, iz;
  int maille_par_cpu, reste, iMaille, neighbour;
  int numberElements, compteMaillesParallele(0), countElements(0);
  int *elements_rec, *elements_env;
  string whichCpuAmIForNeighbour;

  if (m_geometrie == 1) {
    //1D Cartesian Processor Topology
    //-------------------------------

    m_numberCpuX = Ncpu;
    m_CpuCoordX = rankCpu;
    maille_par_cpu = m_numberCellsXGlobal / Ncpu;
    reste = m_numberCellsXGlobal % Ncpu;

    if (rankCpu < reste){ ++maille_par_cpu; }
    m_numberCellsX = maille_par_cpu;
    m_numberCellsY = m_numberCellsYGlobal;
    m_numberCellsZ = m_numberCellsZGlobal;

    m_offsetX = rankCpu * maille_par_cpu;
    if (rankCpu >= reste) { m_offsetX += reste; }

    //Number de mailles sur ce Cpu;
    m_numberCellsCalcul = m_numberCellsX*m_numberCellsY*m_numberCellsZ;

    //Determination du number d'element a envoyer/recevoir
    numberElements = m_numberCellsY*m_numberCellsZ;

    //Number de mailles total en comptant les mailles necessaires au parallele;
    if (rankCpu > 0 && rankCpu < Ncpu - 1){ m_numberCellsTotal = m_numberCellsCalcul + 2 * numberElements; }
    else{ m_numberCellsTotal = m_numberCellsCalcul + numberElements; }

    ////***************Table de connectivite**********

    //Le number des mailles paralleles est situe en dehors du number de maille "vraies"
    //donc commence a numberCells
    compteMaillesParallele = m_numberCellsCalcul;
    elements_rec = new int[numberElements];
    elements_env = new int[numberElements];
    countElements = 0;
    if (rankCpu > 0)
    {
      neighbour = rankCpu - 1;
      whichCpuAmIForNeighbour = "RIGHT";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      ix = 0;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          elements_rec[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_env[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_env, numberElements);
      parallel.setElementsToReceive(neighbour, elements_rec, numberElements);
    }
    countElements = 0;
    if (rankCpu < Ncpu - 1)
    {
      neighbour = rankCpu + 1;
      whichCpuAmIForNeighbour = "LEFT";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      ix = m_numberCellsX - 1;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          elements_rec[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_env[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_env, numberElements);
      parallel.setElementsToReceive(neighbour, elements_rec, numberElements);
    }
  }

  else if (m_geometrie == 2) {
    //2D Cartesian Processor Topology
    //-------------------------------

    //Mininum # of stencils (5) times order of the method (always considered as second order) gives the minimum number of cells in one direction
    int minCellNumberOneDirection(2);

    //Initial estimate of optimal processor topology
    m_numberCpuX = 1;      //Optimal number of processors in the x-direction
    m_numberCpuY = Ncpu;   //Optimal number of processors in the y-direction
    int error(-1);

    //Benchmarking the quality of this initial guess
    double tempNumberCpuX(static_cast<double>(m_numberCpuX)); //Non-optimal number of processors in the x-direction
    double tempNumberCpuY(static_cast<double>(m_numberCpuY)); //Non-optimal number of processors in the y-direction
    double CpuFactorizationMin;                               //Processor factorization (fct) minimization parameter
    CpuFactorizationMin = 10.*abs(static_cast<double>(m_numberCellsXGlobal) / tempNumberCpuX - static_cast<double>(m_numberCellsYGlobal) / tempNumberCpuY);

    //Optimization of the initial processor topology
    for (int i = 1; i <= Ncpu; i++) {
      if ((Ncpu % i == 0) && (m_numberCellsXGlobal / i >= minCellNumberOneDirection)) {
        tempNumberCpuX = i;
        tempNumberCpuY = Ncpu / i;
        if (CpuFactorizationMin >= abs(m_numberCellsXGlobal / tempNumberCpuX -m_numberCellsYGlobal / tempNumberCpuY) && m_numberCellsYGlobal / tempNumberCpuY >= minCellNumberOneDirection) {
          m_numberCpuX = i;
          m_numberCpuY = Ncpu / i;
          CpuFactorizationMin = abs(m_numberCellsXGlobal / tempNumberCpuX - m_numberCellsYGlobal / tempNumberCpuY);
          error = 0;
        }
      }
    }

    //Verifying that a valid decomposition of the computational domain has been established. If not, the simulation exits.
    if (rankCpu == 0 && error == -1) {
      Errors::errorMessage("Unsupported combination of values of number of CPU, cells in x- and y-directions");
    }

    //Coordinates of the current CPU
    int rest;
    rest = rankCpu % (m_numberCpuX*m_numberCpuY);
    m_CpuCoordY = rest / m_numberCpuX;
    m_CpuCoordX = rest % m_numberCpuX;

    //Determination of the number of cells in each direction of the current CPU
    //Number of remaining cells
    int remainingCellsX(m_numberCellsXGlobal % m_numberCpuX);
    int remainingCellsY(m_numberCellsYGlobal % m_numberCpuY);
    m_numberCellsX = m_numberCellsXGlobal / m_numberCpuX;
    m_numberCellsY = m_numberCellsYGlobal / m_numberCpuY;
    m_numberCellsZ = m_numberCellsZGlobal;
    if (m_CpuCoordX < remainingCellsX) { ++m_numberCellsX; }
    if (m_CpuCoordY < remainingCellsY) { ++m_numberCellsY; }

    //Determination of the offset
    m_offsetX = m_CpuCoordX * m_numberCellsX;
    m_offsetY = m_CpuCoordY * m_numberCellsY;
    if (m_CpuCoordX >= remainingCellsX) { m_offsetX += remainingCellsX; }
    if (m_CpuCoordY >= remainingCellsY) { m_offsetY += remainingCellsY; }

    //Number of cells on the current CPU;
    m_numberCellsCalcul = m_numberCellsX * m_numberCellsY * m_numberCellsZ;

    //Determination of the number of elements to send/receive for each parallel boundary in each direction
    int numberElementsX(m_numberCellsY * m_numberCellsZ);
    int numberElementsY(m_numberCellsX * m_numberCellsZ);

    //Total number of cells taking into account the necessary cells for the parallel (ghost cells)
    m_numberCellsTotal = m_numberCellsCalcul;
    if (m_numberCpuX > 1) {
      if (m_CpuCoordX > 0 && m_CpuCoordX < m_numberCpuX - 1) { m_numberCellsTotal += 2 * numberElementsX; }
      else { m_numberCellsTotal += numberElementsX; }
    }
    if (m_numberCpuY > 1) {
      if (m_CpuCoordY > 0 && m_CpuCoordY < m_numberCpuY - 1) { m_numberCellsTotal += 2 * numberElementsY; }
      else { m_numberCellsTotal += numberElementsY; }
    }

    //Connectivity table
    //The number of the ghost cells is located apart from the number of "true" cells, hence start at m_numberCellsCalcul
    compteMaillesParallele = m_numberCellsCalcul;
    int neighbourCpuCoordX, neighbourCpuCoordY;
    //In the x-direction
    int *elements_recX, *elements_envX;
    elements_recX = new int[numberElementsX];
    elements_envX = new int[numberElementsX];
    countElements = 0;
    if (m_CpuCoordX > 0)
    {
      neighbourCpuCoordX = m_CpuCoordX - 1;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      whichCpuAmIForNeighbour = "RIGHT";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      ix = 0;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          elements_recX[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_envX[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_envX, numberElementsX);
      parallel.setElementsToReceive(neighbour, elements_recX, numberElementsX);
    }
    countElements = 0;
    if (m_CpuCoordX < m_numberCpuX - 1)
    {
      neighbourCpuCoordX = m_CpuCoordX + 1;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      whichCpuAmIForNeighbour = "LEFT";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      ix = m_numberCellsX - 1;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          elements_recX[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_envX[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_envX, numberElementsX);
      parallel.setElementsToReceive(neighbour, elements_recX, numberElementsX);
    }
    //In the y-direction
    int *elements_recY, *elements_envY;
    elements_recY = new int[numberElementsY];
    elements_envY = new int[numberElementsY];
    countElements = 0;
    if (m_CpuCoordY > 0)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY - 1;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      whichCpuAmIForNeighbour = "TOP";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      iy = 0;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          elements_recY[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_envY[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_envY, numberElementsY);
      parallel.setElementsToReceive(neighbour, elements_recY, numberElementsY);
    }
    countElements = 0;
    if (m_CpuCoordY < m_numberCpuY - 1)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY + 1;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      whichCpuAmIForNeighbour = "BOTTOM";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      iy = m_numberCellsY - 1;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          elements_recY[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_envY[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_envY, numberElementsY);
      parallel.setElementsToReceive(neighbour, elements_recY, numberElementsY);
    }
  }

  else {
    //3D Cartesian Processor Topology
    //-------------------------------

    //Mininum # of stencils (5) times order of the method (always considered as second order) gives the minimum number of cells in one direction
    int minCellNumberOneDirection(2);

    //Initial estimate of optimal processor topology
    m_numberCpuX = 1;
    m_numberCpuY = 1;
    m_numberCpuZ = Ncpu;
    int error(-1);

    //Benchmarking the quality of this initial guess
    double tempNumberCpuX(static_cast<double>(m_numberCpuX)); //Non-optimal number of processors in the x-direction
    double tempNumberCpuY(static_cast<double>(m_numberCpuY)); //Non-optimal number of processors in the y-direction
    double tempNumberCpuZ(static_cast<double>(m_numberCpuZ)); //Non-optimal number of processors in the z-direction
    double CpuFactorizationMin;                               //Processor factorization (fct) minimization parameter
    CpuFactorizationMin = 10.*abs(static_cast<double>(m_numberCellsXGlobal) / tempNumberCpuX - static_cast<double>(m_numberCellsYGlobal) / tempNumberCpuY)
      + 10.*abs(static_cast<double>(m_numberCellsYGlobal) / tempNumberCpuY - static_cast<double>(m_numberCellsZGlobal) / tempNumberCpuZ);

    //Optimization of the initial processor topology
    for (int i = 1; i <= Ncpu; i++) {
      if ((Ncpu % i == 0) && (m_numberCellsXGlobal / i >= minCellNumberOneDirection)) {
        for (int j = 1; j <= Ncpu/i; j++) {
          if ((Ncpu/i % j == 0) && (m_numberCellsYGlobal / j >= minCellNumberOneDirection)) {
            tempNumberCpuX = i;
            tempNumberCpuY = j;
            tempNumberCpuZ = Ncpu / (i*j);
            if (CpuFactorizationMin >= abs(m_numberCellsXGlobal / tempNumberCpuX - m_numberCellsYGlobal / tempNumberCpuY) + abs(m_numberCellsYGlobal / tempNumberCpuY - m_numberCellsZGlobal / tempNumberCpuZ)
              && m_numberCellsZGlobal / tempNumberCpuZ >= minCellNumberOneDirection) {
              m_numberCpuX = i;
              m_numberCpuY = j;
              m_numberCpuZ = Ncpu / (i*j);
              CpuFactorizationMin = abs(m_numberCellsXGlobal / tempNumberCpuX - m_numberCellsYGlobal / tempNumberCpuY) + abs(m_numberCellsYGlobal / tempNumberCpuY - m_numberCellsZGlobal / tempNumberCpuZ);
              error = 0;
            }
          }
        }
      }
    }

    //Verifying that a valid decomposition of the computational domain has been established. If not, the simulation exits.
    if (rankCpu == 0 && error == -1) {
      Errors::errorMessage("Unsupported combination of values of number of CPU, cells in x- and y-directions");
    }

    //Coordinates of the current CPU
    int rest;
    m_CpuCoordZ = rankCpu / (m_numberCpuX*m_numberCpuY);
    rest = rankCpu % (m_numberCpuX*m_numberCpuY);
    m_CpuCoordY = rest / m_numberCpuX;
    m_CpuCoordX = rest % m_numberCpuX;

    //Determination of the number of cells in each direction of the current CPU
    //Number of remaining cells
    int remainingCellsX(m_numberCellsXGlobal % m_numberCpuX);
    int remainingCellsY(m_numberCellsYGlobal % m_numberCpuY);
    int remainingCellsZ(m_numberCellsZGlobal % m_numberCpuZ);
    m_numberCellsX = m_numberCellsXGlobal / m_numberCpuX;
    m_numberCellsY = m_numberCellsYGlobal / m_numberCpuY;
    m_numberCellsZ = m_numberCellsZGlobal / m_numberCpuZ;
    if (m_CpuCoordX < remainingCellsX) { ++m_numberCellsX; }
    if (m_CpuCoordY < remainingCellsY) { ++m_numberCellsY; }
    if (m_CpuCoordZ < remainingCellsZ) { ++m_numberCellsZ; }

    //Determination of the offset
    m_offsetX = m_CpuCoordX * m_numberCellsX;
    m_offsetY = m_CpuCoordY * m_numberCellsY;
    m_offsetZ = m_CpuCoordZ * m_numberCellsZ;
    if (m_CpuCoordX >= remainingCellsX) { m_offsetX += remainingCellsX; }
    if (m_CpuCoordY >= remainingCellsY) { m_offsetY += remainingCellsY; }
    if (m_CpuCoordZ >= remainingCellsZ) { m_offsetZ += remainingCellsZ; }

    //Number of cells on the current CPU;
    m_numberCellsCalcul = m_numberCellsX * m_numberCellsY * m_numberCellsZ;

    //Determination of the number of elements to send/receive for each parallel boundary in each direction
    int numberElementsX(m_numberCellsY * m_numberCellsZ);
    int numberElementsY(m_numberCellsX * m_numberCellsZ);
    int numberElementsZ(m_numberCellsX * m_numberCellsY);

    //Total number of cells taking into account the necessary cells for the parallel (ghost cells)
    m_numberCellsTotal = m_numberCellsCalcul;
    if (m_numberCpuX > 1) {
      if (m_CpuCoordX > 0 && m_CpuCoordX < m_numberCpuX - 1) { m_numberCellsTotal += 2 * numberElementsX; }
      else { m_numberCellsTotal += numberElementsX; }
    }
    if (m_numberCpuY > 1) {
      if (m_CpuCoordY > 0 && m_CpuCoordY < m_numberCpuY - 1) { m_numberCellsTotal += 2 * numberElementsY; }
      else { m_numberCellsTotal += numberElementsY; }
    }
    if (m_numberCpuZ > 1) {
      if (m_CpuCoordZ > 0 && m_CpuCoordZ < m_numberCpuZ - 1) { m_numberCellsTotal += 2 * numberElementsZ; }
      else { m_numberCellsTotal += numberElementsZ; }
    }

    //Connectivity table
    //The number of the ghost cells is located apart from the number of "true" cells, hence start at m_numberCellsCalcul
    compteMaillesParallele = m_numberCellsCalcul;
    int neighbourCpuCoordX, neighbourCpuCoordY, neighbourCpuCoordZ;
    //In the x-direction
    int *elements_recX, *elements_envX;
    elements_recX = new int[numberElementsX];
    elements_envX = new int[numberElementsX];
    countElements = 0;
    if (m_CpuCoordX > 0)
    {
      neighbourCpuCoordX = m_CpuCoordX - 1;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbourCpuCoordZ = m_CpuCoordZ;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      whichCpuAmIForNeighbour = "RIGHT";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      ix = 0;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          elements_recX[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_envX[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_envX, numberElementsX);
      parallel.setElementsToReceive(neighbour, elements_recX, numberElementsX);
    }
    countElements = 0;
    if (m_CpuCoordX < m_numberCpuX - 1)
    {
      neighbourCpuCoordX = m_CpuCoordX + 1;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbourCpuCoordZ = m_CpuCoordZ;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      whichCpuAmIForNeighbour = "LEFT";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      ix = m_numberCellsX - 1;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          elements_recX[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_envX[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_envX, numberElementsX);
      parallel.setElementsToReceive(neighbour, elements_recX, numberElementsX);
    }
    //In the y-direction
    int *elements_recY, *elements_envY;
    elements_recY = new int[numberElementsY];
    elements_envY = new int[numberElementsY];
    countElements = 0;
    if (m_CpuCoordY > 0)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY - 1;
      neighbourCpuCoordZ = m_CpuCoordZ;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      whichCpuAmIForNeighbour = "TOP";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      iy = 0;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          elements_recY[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_envY[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_envY, numberElementsY);
      parallel.setElementsToReceive(neighbour, elements_recY, numberElementsY);
    }
    countElements = 0;
    if (m_CpuCoordY < m_numberCpuY - 1)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY + 1;
      neighbourCpuCoordZ = m_CpuCoordZ;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      whichCpuAmIForNeighbour = "BOTTOM";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      iy = m_numberCellsY - 1;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          elements_recY[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_envY[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_envY, numberElementsY);
      parallel.setElementsToReceive(neighbour, elements_recY, numberElementsY);
    }
    //In the z-direction
    int *elements_recZ, *elements_envZ;
    elements_recZ = new int[numberElementsZ];
    elements_envZ = new int[numberElementsZ];
    countElements = 0;
    if (m_CpuCoordZ > 0)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbourCpuCoordZ = m_CpuCoordZ - 1;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      whichCpuAmIForNeighbour = "FRONT";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      iz = 0;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iy = 0; iy < m_numberCellsY; iy++)
        {
          elements_recZ[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_envZ[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_envZ, numberElementsZ);
      parallel.setElementsToReceive(neighbour, elements_recZ, numberElementsZ);
    }
    countElements = 0;
    if (m_CpuCoordZ < m_numberCpuZ - 1)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbourCpuCoordZ = m_CpuCoordZ + 1;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      whichCpuAmIForNeighbour = "BACK";
      parallel.setNeighbour(neighbour, whichCpuAmIForNeighbour);
      iz = m_numberCellsZ - 1;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iy = 0; iy < m_numberCellsY; iy++)
        {
          elements_recZ[countElements] = compteMaillesParallele;
          this->construitIGlobal(ix, iy, iz, iMaille);
          elements_envZ[countElements] = iMaille;
          ++countElements;
          ++compteMaillesParallele;
        }
      }
      parallel.setElementsToSend(neighbour, elements_envZ, numberElementsZ);
      parallel.setElementsToReceive(neighbour, elements_recZ, numberElementsZ);
    }
  }
}

//***********************************************************************

string MeshCartesian::whoAmI() const
{
  return "CARTESIAN";
}

//**************************************************************************
//******************************** PRINTING ********************************
//**************************************************************************

//****************************************************************************

string MeshCartesian::recupereChaineExtent(int localRank, bool global) const
{
  stringstream chaineExtent;

  //Determination of the offset of the current CPU in each direction => determination of the first and last element of the current CPU in the complete domain
  //Number of remaining cells
  int remainingCellsX(m_numberCellsXGlobal % m_numberCpuX);
  int remainingCellsY(m_numberCellsYGlobal % m_numberCpuY);
  int remainingCellsZ(m_numberCellsZGlobal % m_numberCpuZ);
  int offset_x(m_CpuCoordX*m_numberCellsX);
  int offset_y(m_CpuCoordY*m_numberCellsY);
  int offset_z(m_CpuCoordZ*m_numberCellsZ);
  if (m_CpuCoordX >= remainingCellsX) { offset_x += remainingCellsX; }
  if (m_CpuCoordY >= remainingCellsY) { offset_y += remainingCellsY; }
  if (m_CpuCoordZ >= remainingCellsZ) { offset_z += remainingCellsZ; }
  
  if (global) { chaineExtent << " 0 " << m_numberCellsXGlobal << " 0 " << m_numberCellsYGlobal << " 0 " << m_numberCellsZGlobal; }
  else{ chaineExtent << offset_x << " " << offset_x + m_numberCellsX << " " << offset_y << " " << offset_y + m_numberCellsY << " " << offset_z << " " << offset_z + m_numberCellsZ; }
  
  return chaineExtent.str();
}

//****************************************************************************

void MeshCartesian::recupereCoord(vector<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, Axe axe) const
{
  switch (axe) {
  case X:
    for (int i = 0; i < m_numberCellsX; i++)
    {
      jeuDonnees.push_back(m_posXi[m_offsetX + i] - 0.5*m_dXi[m_offsetX + i]);
    }
    jeuDonnees.push_back(m_posXi[m_offsetX + m_numberCellsX - 1] + 0.5*m_dXi[m_offsetX + m_numberCellsX - 1]);
    break;
  case Y:
    for (int j = 0; j < m_numberCellsY; j++)
    {
      jeuDonnees.push_back(m_posYj[m_offsetY + j] - 0.5*m_dYj[m_offsetY + j]);
    }
    jeuDonnees.push_back(m_posYj[m_offsetY + m_numberCellsY - 1] + 0.5*m_dYj[m_offsetY + m_numberCellsY - 1]);
    break;
  case Z:
    for (int k = 0; k < m_numberCellsZ; k++)
    {
      jeuDonnees.push_back(m_posZk[m_offsetZ + k] - 0.5*m_dZk[m_offsetZ + k]);
    }
    jeuDonnees.push_back(m_posZk[m_offsetZ + m_numberCellsZ - 1] + 0.5*m_dZk[m_offsetZ + m_numberCellsZ - 1]);
    break;
  }
}

//****************************************************************************

void MeshCartesian::recupereDonnees(vector<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, const int var, int phase, int lvl) const
{
  jeuDonnees.clear();
  int numCell;
  for (int k = 0; k < m_numberCellsZ; k++) {
    for (int j = 0; j < m_numberCellsY; j++) {
      for (int i = 0; i < m_numberCellsX; i++) {
        construitIGlobal(i, j, k, numCell);
        if (var > 0) { //On veut recuperer les donnees scalars
          if (phase >= 0) { jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnScalar(var)); }     //Donnees de phases
          else if (phase == -1) { jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnScalar(var)); }  //Donnees de mixture
          else if (phase == -2) { jeuDonnees.push_back(cellsLvl[0][numCell]->getTransport(var - 1).getValue()); }
          else if (phase == -3) { jeuDonnees.push_back(cellsLvl[0][numCell]->getXi()); }
          else if (phase == -4) { jeuDonnees.push_back(cellsLvl[0][numCell]->getGradient()); }
          else { Errors::errorMessage("MeshCartesian::recupereDonnees: unknown number of phase: ", phase); }
        }
        else { //On veut recuperer les donnees vectorielles
          if (phase >= 0) { //Donnees de phases
            jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getX());
            jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getY());
            jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getZ());
          }
          else if (phase == -1) {  //Donnees de mixture
            jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getX());
            jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getY());
            jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getZ());
          }
          else { Errors::errorMessage("MeshCartesian::recupereDonnees: unknown number of phase: ", phase); }
        } //Fin vecteur
      } // Fin X
    } //Fin Y
  } //Fin Z
}

//****************************************************************************

void MeshCartesian::setDataSet(std::vector<double> &jeuDonnees, vector<Cell *> *cellsLvl, const int var, int phase, int lvl) const
{
  int iterDataSet(0);
  int numCell;
  Coord vec;
  for (int k = 0; k < m_numberCellsZ; k++) {
    for (int j = 0; j < m_numberCellsY; j++) {
      for (int i = 0; i < m_numberCellsX; i++) {
        construitIGlobal(i, j, k, numCell);
        if (var > 0) { //Scalars data are first set
          if (phase >= 0) { cellsLvl[0][numCell]->getPhase(phase)->setScalar(var, jeuDonnees[iterDataSet++]); } //phases data
          else if (phase == -1) { cellsLvl[0][numCell]->getMixture()->setScalar(var, jeuDonnees[iterDataSet++]); }  //mixture data
          else if (phase == -2) { cellsLvl[0][numCell]->getTransport(var - 1).setValue(jeuDonnees[iterDataSet++]); } //transport data
          else { Errors::errorMessage("MeshCartesian::setDataSet: unknown phase number: ", phase); }
        }
        else { //On veut recuperer les donnees vectorielles
          if (phase >= 0) { //Phases data
            vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet + 1], jeuDonnees[iterDataSet + 2]);
            cellsLvl[0][numCell]->getPhase(phase)->setVector(-var, vec);
            iterDataSet += 3;
          }
          else if (phase == -1) {  //Mixture data
            vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet + 1], jeuDonnees[iterDataSet + 2]);
            cellsLvl[0][numCell]->getMixture()->setVector(-var, vec);
            iterDataSet += 3;
          }
          else { Errors::errorMessage("MeshCartesian::setDataSet: unknown phase number: ", phase); }
        } //Fin vecteur
      } // Fin X
    } // Fin Y
  } // Fin Z
}

//***********************************************************************