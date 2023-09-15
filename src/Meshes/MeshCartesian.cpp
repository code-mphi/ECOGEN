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

#include "MeshCartesian.h"

//***********************************************************************

MeshCartesian::MeshCartesian(double lX, int numberCellsX, double lY, int numberCellsY, double lZ, int numberCellsZ,
  std::vector<stretchZone> stretchX, std::vector<stretchZone> stretchY, std::vector<stretchZone> stretchZ) :
  m_lX(lX), m_lY(lY), m_lZ(lZ),
  m_numberCellsXGlobal(numberCellsX), m_numberCellsYGlobal(numberCellsY), m_numberCellsZGlobal(numberCellsZ)
{
  for (unsigned int i = 0; i < stretchX.size(); i++) { m_stretchX.push_back(stretchX[i]); }
  for (unsigned int i = 0; i < stretchY.size(); i++) { m_stretchY.push_back(stretchY[i]); }
  for (unsigned int i = 0; i < stretchZ.size(); i++) { m_stretchZ.push_back(stretchZ[i]); }

  m_numberCellsCalcul = m_numberCellsXGlobal*m_numberCellsYGlobal*m_numberCellsZGlobal;
  m_problemDimension = 0;
  if (numberCellsX != 1) { m_problemDimension += 1; }
  if (numberCellsY != 1) { m_problemDimension += 1; }
  if (numberCellsZ != 1) { m_problemDimension += 1; }
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

void MeshCartesian::assignLimits(std::vector<BoundCond*>& boundCond)
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
      Errors::errorMessage("Probleme de limites dans assignLimits"); break;
    }
  }

  for (int j = m_numberBoundCondInit; j < 6; j++) {
    switch (j) {
    case 0:
      m_limXm = new BoundCondNonReflecting(j+1); break;
    case 1:
      m_limXp = new BoundCondNonReflecting(j+1); break;
    case 2:
      m_limYm = new BoundCondNonReflecting(j+1); break;
    case 3:
      m_limYp = new BoundCondNonReflecting(j+1); break;
    case 4:
      m_limZm = new BoundCondNonReflecting(j+1); break;
    case 5:
      m_limZp = new BoundCondNonReflecting(j+1); break;
    default:
      Errors::errorMessage("Probleme de limites dans assignLimits"); break;
    }
  }
}

//***********************************************************************

void MeshCartesian::getIJK(const int& index, int& i, int& j, int& k) const
{
  int reste;
  k = index / (m_numberCellsX*m_numberCellsY);
  reste = index % (m_numberCellsX*m_numberCellsY);
  j = reste / m_numberCellsX;
  i = reste % m_numberCellsX;
}

//***********************************************************************

void MeshCartesian::construitIGlobal(const int& i, const int& j, const int& k, int& index) const
{
  index = 0;
  if (m_numberCellsX != 1) index += i;
  if (m_numberCellsY != 1) index += j*m_numberCellsX;
  if (m_numberCellsZ != 1) index += k*m_numberCellsX*m_numberCellsY;
}

//***********************************************************************

int MeshCartesian::initializeGeometrie(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<Cell*>& cellsGhost, TypeMeshContainer<CellInterface*>& cellInterfaces,
  const int& /*restartSimulation*/, bool /*pretraitementParallele*/, std::string ordreCalcul)
{
  this->meshStretching();
  if (Ncpu == 1)
  {
    this->initializeGeometrieMonoCpu(cells, cellInterfaces, ordreCalcul);
  }
  else
  {
    this->initializeGeometrieParallele(cells, cellsGhost, cellInterfaces, ordreCalcul);
  }
  return m_problemDimension;
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

void MeshCartesian::initializeGeometrieMonoCpu(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<CellInterface*>& cellInterfaces, std::string ordreCalcul)
{
  int ix, iy, iz;

  m_numberCellsX = m_numberCellsXGlobal;
  m_numberCellsY = m_numberCellsYGlobal;
  m_numberCellsZ = m_numberCellsZGlobal;
  m_numberCellsTotal = m_numberCellsCalcul;

  //Declaration array of elements and cells
  //---------------------------------------
  for (int i = 0; i < m_numberCellsCalcul; i++) {
    if(ordreCalcul == "FIRSTORDER") { cells.push_back(new Cell); }
    else { cells.push_back(new CellO2Cartesian); }
    m_elements.push_back(new ElementCartesian());
    cells[i]->setElement(m_elements[i], i); 
  }

  //Assign Cartesian geometric data to element
  //------------------------------------------
  Coord tangent, normal, binormal;
  double surface(1.), volume(0.);
  double posX, posY, posZ;
  for (int i = 0; i < m_numberCellsCalcul; i++)
  {
    this->getIJK(i, ix, iy, iz);
    volume = m_dXi[ix] * m_dYj[iy] * m_dZk[iz];
    cells[i]->getElement()->setVolume(volume);

    //CFL lenght
    double lCFL(1.e10);
    if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[ix]); }
    if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[iy]); }
    if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[iz]); }
    if (m_problemDimension > 1) lCFL *= 0.6;

    cells[i]->getElement()->setLCFL(lCFL);
    cells[i]->getElement()->setPos(m_posXi[ix], m_posYj[iy], m_posZk[iz]);
    cells[i]->getElement()->setSize(m_dXi[ix], m_dYj[iy], m_dZk[iz]);
  }

  //Assign Cartesian geometric data to faces
  //----------------------------------------
  //Determination of number of faces for 1D/2D/3D computation
  //m_numberFaces : number total of faces
  //numberFacesLimites : number of faces on limits
  //numberFacesInternes : number of faces common to two elements
  m_numberFacesTotal = 0;
  // int numberFacesLimites(0);
  if (m_numberCellsX != 1)
  {
    m_numberFacesTotal += (m_numberCellsX + 1)*m_numberCellsY*m_numberCellsZ;
    // numberFacesLimites += 2 * m_numberCellsY*m_numberCellsZ;
  }
  if (m_numberCellsY != 1)
  {
    m_numberFacesTotal += (m_numberCellsY + 1)*m_numberCellsX*m_numberCellsZ;
    // numberFacesLimites += 2 * m_numberCellsX*m_numberCellsZ;
  }
  if (m_numberCellsZ != 1)
  {
    m_numberFacesTotal += (m_numberCellsZ + 1)*m_numberCellsX*m_numberCellsY;
    // numberFacesLimites += 2 * m_numberCellsX*m_numberCellsY;
  }
  //int numberFacesInternes(m_numberFacesTotal - numberFacesLimites);

  //Initialization des faces internes
  //*********************************
  int iCellL, iCellR, iFace(0), iTemp;
  //Faces selon X
  tangent.setXYZ(0., 1., 0.); normal.setXYZ(1., 0., 0.); binormal.setXYZ(0., 0., 1.);
  for (ix = 0; ix < m_numberCellsX - 1; ix++)
  {
    for (iy = 0; iy < m_numberCellsY; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
        else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        this->construitIGlobal(ix + 1, iy, iz, iCellR);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dYj[iy] * m_dZk[iz];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(0., m_dYj[iy], m_dZk[iz]);
        posX = m_posXi[ix] + 0.5*m_dXi[ix];
        posY = m_posYj[iy];
        posZ = m_posZk[iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
        //Preparation des cells pour ordre 2 multislopes
        if (ix != 0) { this->construitIGlobal(ix - 1, iy, iz, iTemp); }
        if (ix < m_numberCellsX - 1) { this->construitIGlobal(ix + 1, iy, iz, iTemp); }
        if (ix < m_numberCellsX - 2) { this->construitIGlobal(ix + 2, iy, iz, iTemp); }
        this->construitIGlobal(ix, iy, iz, iTemp);
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
        if(ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
        else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        this->construitIGlobal(ix, iy + 1, iz, iCellR);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dXi[ix] * m_dZk[iz];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[ix], 0., m_dZk[iz]);
        posX = m_posXi[ix];
        posY = m_posYj[iy] + 0.5*m_dYj[iy];
        posZ = m_posZk[iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
        //Preparation des cells pour ordre 2 multislopes
        if (iy != 0) { this->construitIGlobal(ix, iy - 1, iz, iTemp); }
        if (iy < m_numberCellsY - 1) { this->construitIGlobal(ix, iy + 1, iz, iTemp); }
        if (iy < m_numberCellsY - 2) { this->construitIGlobal(ix, iy + 2, iz, iTemp); }
        this->construitIGlobal(ix, iy, iz, iTemp);
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
        if(ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
        else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        this->construitIGlobal(ix, iy, iz + 1, iCellR);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dXi[ix] * m_dYj[iy];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[ix], m_dYj[iy], 0.);
        posX = m_posXi[ix];
        posY = m_posYj[iy];
        posZ = m_posZk[iz] + 0.5*m_dZk[iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
        //Preparation des cells pour ordre 2 multislopes
        if (iz != 0) { this->construitIGlobal(ix, iy, iz - 1, iTemp); }
        if (iz < m_numberCellsZ - 1) { this->construitIGlobal(ix, iy, iz + 1, iTemp); }
        if (iz < m_numberCellsZ - 2) { this->construitIGlobal(ix, iy, iz + 2, iTemp); }
        this->construitIGlobal(ix, iy, iz, iTemp);
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
        m_limXm->createBoundary(cellInterfaces);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        iCellR = iCellL;
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dYj[iy] * m_dZk[iz];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal); //FP//Q// Interet d'avoir deux fonctions initialize ?
        m_faces[iFace]->setSize(0., m_dYj[iy], m_dZk[iz]);
        posX = 0.;
        posY = m_posYj[iy];
        posZ = m_posZk[iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
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
        m_limXp->createBoundary(cellInterfaces);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        iCellR = iCellL;
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dYj[iy] * m_dZk[iz];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(0., m_dYj[iy], m_dZk[iz]);
        posX = m_posXi[ix] + 0.5*m_dXi[ix];
        posY = m_posYj[iy];
        posZ = m_posZk[iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
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
        this->construitIGlobal(ix, iy, iz, iCellL);
        iCellR = iCellL;

        // Hard-coded boundary condition for Blasius test case
        // ---------------------------------------------------
        // Allows to set slip and no-slip wall on same boundary of a Cartesian mesh.
        // Requires to set wall boundary condition in file initialConditions.xml 
        // to define no-slip wall and slipping wall is defined thanks to a
        // symmetry before a given position.
        // To comment if not needed and be carefull when using it.
        // if (cells[iCellL]->getPosition().getX() < 0.1) {
        //   BoundCond* limBuff(new BoundCondSymmetry(4));
        //   limBuff->createBoundary(cellInterfaces);
        // }
        // else {
        //   m_limYm->createBoundary(cellInterfaces);
        // }
        //JC//WARNING This hard-coded boundary condition is not particularly useful anymore
        // since the same effect can be achieved with an unstructured mesh.

        m_limYm->createBoundary(cellInterfaces);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dXi[ix] * m_dZk[iz];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[ix], 0., m_dZk[iz]);
        posX = m_posXi[ix];
        posY = 0.;
        posZ = m_posZk[iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
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
        m_limYp->createBoundary(cellInterfaces);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        iCellR = iCellL;
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dXi[ix] * m_dZk[iz];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[ix], 0., m_dZk[iz]);
        posX = m_posXi[ix];
        posY = m_posYj[iy] + 0.5*m_dYj[iy];
        posZ = m_posZk[iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
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
        m_limZm->createBoundary(cellInterfaces);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        iCellR = iCellL;
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dXi[ix] * m_dYj[iy];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[ix], m_dYj[iy], 0.);
        posX = m_posXi[ix];
        posY = m_posYj[iy];
        posZ = 0.;
        m_faces[iFace]->setPos(posX, posY, posZ);
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
        m_limZp->createBoundary(cellInterfaces);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        iCellR = iCellL;
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dXi[ix] * m_dYj[iy];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[ix], m_dYj[iy], 0.);
        posX = m_posXi[ix];
        posY = m_posYj[iy];
        posZ = m_posZk[iz] + 0.5*m_dZk[iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }
}

//***********************************************************************

void MeshCartesian::initializeGeometrieParallele(TypeMeshContainer<Cell*>& cells, TypeMeshContainer<Cell*>& cellsGhost, TypeMeshContainer<CellInterface*>& cellInterfaces, std::string ordreCalcul)
{
  int ix, iy, iz;
  int countParallelCells(0);

  m_numberCellsX = m_numberCellsXGlobal;
  m_numberCellsY = m_numberCellsYGlobal;
  m_numberCellsZ = m_numberCellsZGlobal;

  //Declaration of the cell array
  //-----------------------------
  //Parallel cutting of the geometric domain
  this->decoupageParallele(ordreCalcul, cells);
  
  //Geometrical data settings for computational cells
  //-------------------------------------------------
  Coord tangent, normal, binormal;
  double surface(1.), volume(0.);
  double posX, posY, posZ;
  for (int i = 0; i < m_numberCellsCalcul; i++)
  {
    this->getIJK(i, ix, iy, iz);
    volume = m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz];
    cells[i]->getElement()->setVolume(volume);

    //CFL lenght
    double lCFL(1.e10);
    if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[m_offsetX + ix]); }
    if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[m_offsetY + iy]); }
    if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[m_offsetZ + iz]); }
    if (m_problemDimension > 1) lCFL *= 0.6;

    cells[i]->getElement()->setLCFL(lCFL);
    cells[i]->getElement()->setPos(m_posXi[m_offsetX + ix], m_posYj[m_offsetY + iy], m_posZk[m_offsetZ + iz]);
    cells[i]->getElement()->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
  }

  //Assign Cartesian geometric data to faces
  //----------------------------------------
  //Determination of number of faces for 1D/2D/3D computation
  //m_numberFaces : number total of faces
  //numberFacesLimites : number of faces on limits
  //numberFacesInternes : number of faces common to two elements
  m_numberFacesTotal = 0;
  // int numberFacesLimites(0);
  if (m_numberCellsX != 1)
  {
    m_numberFacesTotal += (m_numberCellsX + 1)*m_numberCellsY*m_numberCellsZ;
    // numberFacesLimites += 2 * m_numberCellsY*m_numberCellsZ;
  }
  if (m_numberCellsY != 1)
  {
    m_numberFacesTotal += (m_numberCellsY + 1)*m_numberCellsX*m_numberCellsZ;
    // numberFacesLimites += 2 * m_numberCellsX*m_numberCellsZ;
  }
  if (m_numberCellsZ != 1)
  {
    m_numberFacesTotal += (m_numberCellsZ + 1)*m_numberCellsX*m_numberCellsY;
    // numberFacesLimites += 2 * m_numberCellsX*m_numberCellsY;
  }
  //int numberFacesInternes(m_numberFacesTotal - numberFacesLimites);
  
  //Initialization des faces internes
  //*********************************
  int iCellL, iCellR, iFace(0);
  //Faces selon X
  tangent.setXYZ(0., 1., 0.); normal.setXYZ(1., 0., 0.); binormal.setXYZ(0., 0., 1.);
  for (ix = 0; ix < m_numberCellsX - 1; ix++)
  {
    for (iy = 0; iy < m_numberCellsY; iy++)
    {
      for (iz = 0; iz < m_numberCellsZ; iz++)
      {
        if(ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
        else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        this->construitIGlobal(ix + 1, iy, iz, iCellR);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(0., m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix] + 0.5*m_dXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
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
        if(ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
        else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        this->construitIGlobal(ix, iy + 1, iz, iCellR);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dXi[m_offsetX + ix] * m_dZk[m_offsetZ + iz];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[m_offsetX + ix], 0., m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy] + 0.5*m_dYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
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
        if(ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
        else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        this->construitIGlobal(ix, iy, iz, iCellL);
        this->construitIGlobal(ix, iy, iz + 1, iCellR);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        surface = m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy];
        m_faces[iFace]->initializeOthers(surface, normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], 0.);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz] + 0.5*m_dZk[m_offsetZ + iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }

  countParallelCells = m_numberCellsCalcul;

  //Ghosts cells and cell interfaces initialisation in X direction
  //**************************************************************
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
					this->construitIGlobal(ix, iy, iz, iCellR);
          iCellL = countParallelCells++; //Ghost cell number taken in order
          //setting ghost cell geometry
					cells[iCellL]->getElement()->setPos(cells[iCellR]->getElement()->getPosition());
          cells[iCellL]->getElement()->setPosX(m_posXi[m_offsetX - 1]);
          cells[iCellL]->getElement()->setSize(m_dXi[m_offsetX -1], m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
          cells[iCellL]->getElement()->setVolume(m_dXi[m_offsetX - 1] * m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[m_offsetX - 1]); }
          if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[m_offsetY + iy]); }
          if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[m_offsetZ + iz]); }
          if (m_problemDimension > 1) lCFL *= 0.6;
          cells[iCellL]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
          else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }     
					cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        }
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          tangent.setXYZ(0., -1., 0.); normal.setXYZ(-1., 0., 0.); 
          //right and Left cells equals
					this->construitIGlobal(ix, iy, iz, iCellL);
          iCellR = iCellL;
          //setting boundary
          m_limXm->createBoundary(cellInterfaces);
        }
        //Common settings
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        m_faces[iFace]->initializeOthers(m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz], normal, tangent, binormal);
        m_faces[iFace]->setSize(0., m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix] - 0.5*m_dXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
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
        this->construitIGlobal(ix, iy, iz, iCellL);
        //A) CPU neighbour limits treatment
        //---------------------------------
        if (m_numberCpuX > 1 && m_CpuCoordX < m_numberCpuX - 1) {
          iCellR = countParallelCells++; //Ghost cell number taken in order
          //setting ghost cell geometry
          cells[iCellR]->getElement()->setPos(cells[iCellL]->getElement()->getPosition());
          cells[iCellR]->getElement()->setPosX(m_posXi[m_offsetX + ix + 1]);
          cells[iCellR]->getElement()->setSize(m_dXi[m_offsetX + ix + 1], m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
          cells[iCellR]->getElement()->setVolume(m_dXi[m_offsetX + ix + 1] * m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[m_offsetX + ix + 1]); }
          if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[m_offsetY + iy]); }
          if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[m_offsetZ + iz]); }
          if (m_problemDimension > 1) lCFL *= 0.6;
          cells[iCellR]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
          else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
          cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        }
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          m_limXp->createBoundary(cellInterfaces);
          iCellR = iCellL;
        }
        //Common settings
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        m_faces[iFace]->initializeOthers(m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz], normal, tangent, binormal);
        m_faces[iFace]->setSize(0., m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix] + 0.5*m_dXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }

  //Ghosts cells and cell interfaces initialisation in Y direction
  //**************************************************************
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
					this->construitIGlobal(ix, iy, iz, iCellR);
          iCellL = countParallelCells++; //Ghost cell number taken in order
          //setting ghost cell geometry
          cells[iCellL]->getElement()->setPos(cells[iCellR]->getElement()->getPosition());
          cells[iCellL]->getElement()->setPosY(m_posYj[m_offsetY - 1]);
          cells[iCellL]->getElement()->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY - 1], m_dZk[m_offsetZ + iz]);
          cells[iCellL]->getElement()->setVolume(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY - 1] * m_dZk[m_offsetZ + iz]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[m_offsetX + ix]); }
          if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[m_offsetY - 1]); }
          if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[m_offsetZ + iz]); }
          if (m_problemDimension > 1) lCFL *= 0.6;
          cells[iCellL]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
          else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
          cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        } 
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          tangent.setXYZ(1., 0., 0.); normal.setXYZ(0., -1., 0.);
          //right and Left cells equals
          this->construitIGlobal(ix, iy, iz, iCellL);
          iCellR = iCellL;

          // Hard-coded boundary condition for Blasius test case
          // ---------------------------------------------------
          // Allows to set slip and no-slip wall on same boundary of a Cartesian mesh.
          // Requires to set wall boundary condition in file initialConditions.xml 
          // to define no-slip wall and slipping wall is defined thanks to a
          // symmetry before a given position.
          // To comment if not needed and be carefull when using it.
          // if (cells[iCellL]->getPosition().getX() < 0.1) {
          //   BoundCond* limBuff(new BoundCondSymmetry(4));
          //   limBuff->createBoundary(cellInterfaces);
          // }
          // else {
          //   m_limYm->createBoundary(cellInterfaces);
          // }
          //JC//WARNING This hard-coded boundary condition is not particularly useful anymore
          // since the same effect can be achieved with an unstructured mesh.

          //setting boundary
          m_limYm->createBoundary(cellInterfaces);
        }
        //Common settings
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        m_faces[iFace]->initializeOthers(m_dXi[m_offsetX + ix] * m_dZk[m_offsetZ + iz], normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[m_offsetX + ix], 0., m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy] - 0.5*m_dYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
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
        this->construitIGlobal(ix, iy, iz, iCellL);
        //A) CPU neighbour limits treatment
        //---------------------------------
        if (m_numberCpuY > 1 && m_CpuCoordY < m_numberCpuY - 1) {
          iCellR = countParallelCells++; //Ghost cell number taken in order
          //setting ghost cell geometry
          cells[iCellR]->getElement()->setPos(cells[iCellL]->getElement()->getPosition());
          cells[iCellR]->getElement()->setPosY(m_posYj[m_offsetY + iy + 1]);
          cells[iCellR]->getElement()->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy + 1], m_dZk[m_offsetZ + iz]);
          cells[iCellR]->getElement()->setVolume(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy + 1] * m_dZk[m_offsetZ + iz]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[m_offsetX + ix]); }
          if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[m_offsetY + iy + 1]); }
          if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[m_offsetZ + iz]); }
          if (m_problemDimension > 1) lCFL *= 0.6;
          cells[iCellR]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
          else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
          cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        }
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          m_limYp->createBoundary(cellInterfaces);
          iCellR = iCellL;
        }
        //Common settings
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        m_faces[iFace]->initializeOthers(m_dXi[m_offsetX + ix] * m_dZk[m_offsetZ + iz], normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[m_offsetX + ix], 0., m_dZk[m_offsetZ + iz]);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy] + 0.5*m_dYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }

  //Ghosts cells and cell interfaces initialisation in Z direction
  //**************************************************************
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
          this->construitIGlobal(ix, iy, iz, iCellR);
          iCellL = countParallelCells++; //Ghost cell number taken in order
          //setting ghost cell geometry
          cells[iCellL]->getElement()->setPos(cells[iCellR]->getElement()->getPosition());
          cells[iCellL]->getElement()->setPosZ(m_posZk[m_offsetZ - 1]);
          cells[iCellL]->getElement()->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], m_dZk[m_offsetZ - 1]);
          cells[iCellL]->getElement()->setVolume(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ - 1]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[m_offsetX + ix]); }
          if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[m_offsetY + iy]); }
          if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[m_offsetZ - 1]); }
          if (m_problemDimension > 1) lCFL *= 0.6;
          cells[iCellL]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
          else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
          cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        }
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          tangent.setXYZ(-1., 0., 0.); normal.setXYZ(0., 0., -1.);
          //right and Left cells equals
          this->construitIGlobal(ix, iy, iz, iCellL);
          iCellR = iCellL;
          //setting boundary
          m_limZm->createBoundary(cellInterfaces);
        }
        //Common settings
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        m_faces[iFace]->initializeOthers(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy], normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], 0.);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz] - 0.5*m_dZk[m_offsetZ + iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
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
        this->construitIGlobal(ix, iy, iz, iCellL);
        //A) CPU neighbour limits treatment
        //---------------------------------
        if (m_numberCpuZ > 1 && m_CpuCoordZ < m_numberCpuZ - 1) {
          iCellR = countParallelCells++; //Ghost cell number taken in order
          //setting ghost cell geometry
          cells[iCellR]->getElement()->setPos(cells[iCellL]->getElement()->getPosition());
          cells[iCellR]->getElement()->setPosZ(m_posZk[m_offsetZ + iz + 1]);
          cells[iCellR]->getElement()->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], m_dZk[m_offsetZ + iz + 1]);
          cells[iCellR]->getElement()->setVolume(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy] * m_dZk[m_offsetZ + iz + 1]);
          double lCFL(1.e10);
          if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[m_offsetX + ix]); }
          if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[m_offsetY + iy]); }
          if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[m_offsetZ + iz + 1]); }
          if (m_problemDimension > 1) lCFL *= 0.6;
          cells[iCellR]->getElement()->setLCFL(lCFL);
          //setting boundary
          if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
          else { cellInterfaces.push_back(new CellInterfaceO2Cartesian); }
          cells[iCellR]->addCellInterface(cellInterfaces[iFace]);
        }
        //B) Physical boundary condition treatment
        //----------------------------------------
        else {
          m_limZp->createBoundary(cellInterfaces);
          iCellR = iCellL;
        }
        //Common settings
        cells[iCellL]->addCellInterface(cellInterfaces[iFace]);
        m_faces.push_back(new FaceCartesian());
        cellInterfaces[iFace]->setFace(m_faces[iFace]);
        cellInterfaces[iFace]->initialize(cells[iCellL], cells[iCellR]);
        m_faces[iFace]->initializeOthers(m_dXi[m_offsetX + ix] * m_dYj[m_offsetY + iy], normal, tangent, binormal);
        m_faces[iFace]->setSize(m_dXi[m_offsetX + ix], m_dYj[m_offsetY + iy], 0.);
        posX = m_posXi[m_offsetX + ix];
        posY = m_posYj[m_offsetY + iy];
        posZ = m_posZk[m_offsetZ + iz] + 0.5*m_dZk[m_offsetZ + iz];
        m_faces[iFace]->setPos(posX, posY, posZ);
        iFace++;
      }
    }
  }

  //Update of cellsGhost
  cellsGhost.insert(cellsGhost.begin(), cells.begin()+m_numberCellsCalcul, cells.end());
  cells.erase(cells.begin()+m_numberCellsCalcul, cells.end());
}

//***********************************************************************

void MeshCartesian::decoupageParallele(std::string ordreCalcul, TypeMeshContainer<Cell*>& cells)
{
  int ix, iy, iz;
  int cellPerCpu, reste, iCell, neighbour;
  int numberElements, countParallelCells(0);
  
  if (m_problemDimension == 1) {
    //1D Cartesian Processor Topology
    //-------------------------------

    m_numberCpuX = Ncpu;
    m_CpuCoordX = rankCpu;
    cellPerCpu = m_numberCellsXGlobal / Ncpu;
    reste = m_numberCellsXGlobal % Ncpu;

    if (rankCpu < reste){ ++cellPerCpu; }
    m_numberCellsX = cellPerCpu;
    m_numberCellsY = m_numberCellsYGlobal;
    m_numberCellsZ = m_numberCellsZGlobal;

    m_offsetX = rankCpu * cellPerCpu;
    if (rankCpu >= reste) { m_offsetX += reste; }

    //Number of cells on this CPU
    m_numberCellsCalcul = m_numberCellsX*m_numberCellsY*m_numberCellsZ;

    //Determination du number d'element a envoyer/recevoir
    numberElements = m_numberCellsY*m_numberCellsZ;

    //Number of total cells counting the cells necessary for parallel;
    if (rankCpu > 0 && rankCpu < Ncpu - 1){ m_numberCellsTotal = m_numberCellsCalcul + 2 * numberElements; }
    else{ m_numberCellsTotal = m_numberCellsCalcul + numberElements; }


    //Generating cells 
    for (int i = 0; i < m_numberCellsCalcul; i++) {
        if (ordreCalcul == "FIRSTORDER") { cells.push_back(new Cell); }
        else { cells.push_back(new CellO2Cartesian); }
        m_elements.push_back(new ElementCartesian());
        cells[i]->setElement(m_elements[i], i);
    }
    for (int i = m_numberCellsCalcul; i < m_numberCellsTotal; i++) {
        if (ordreCalcul == "FIRSTORDER") { cells.push_back(new CellGhost); }
        else { cells.push_back(new CellO2GhostCartesian); }
        m_elements.push_back(new ElementCartesian());
        cells[i]->setElement(m_elements[i], i);
        cells[i]->pushBackSlope();
    }

    ////***************Connectivity table**********

    //The number of parallel cells is located apart from of the "true" cells
    //so we start at numberCells
    countParallelCells = m_numberCellsCalcul;
    if (rankCpu > 0)
    {
      neighbour = rankCpu - 1;
      parallel.setNeighbour(neighbour);
      ix = 0;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
    if (rankCpu < Ncpu - 1)
    {
      neighbour = rankCpu + 1;
      parallel.setNeighbour(neighbour);
      ix = m_numberCellsX - 1;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
  }

  else if (m_problemDimension == 2) {
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
    CpuFactorizationMin = 10.*std::fabs(static_cast<double>(m_numberCellsXGlobal) / tempNumberCpuX - static_cast<double>(m_numberCellsYGlobal) / tempNumberCpuY);

    //Optimization of the initial processor topology
    for (int i = 1; i <= Ncpu; i++) {
      if ((Ncpu % i == 0) && (m_numberCellsXGlobal / i >= minCellNumberOneDirection)) {
        tempNumberCpuX = i;
        tempNumberCpuY = Ncpu / i;
        if (CpuFactorizationMin >= std::fabs(m_numberCellsXGlobal / tempNumberCpuX -m_numberCellsYGlobal / tempNumberCpuY) && m_numberCellsYGlobal / tempNumberCpuY >= minCellNumberOneDirection) {
          m_numberCpuX = i;
          m_numberCpuY = Ncpu / i;
          CpuFactorizationMin = std::fabs(m_numberCellsXGlobal / tempNumberCpuX - m_numberCellsYGlobal / tempNumberCpuY);
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

    //Generating cells 
    for (int i = 0; i < m_numberCellsCalcul; i++) {
        if (ordreCalcul == "FIRSTORDER") { cells.push_back(new Cell); }
        else { cells.push_back(new CellO2Cartesian); }
        m_elements.push_back(new ElementCartesian());
        cells[i]->setElement(m_elements[i], i);
    }
    for (int i = m_numberCellsCalcul; i < m_numberCellsTotal; i++) {
        if (ordreCalcul == "FIRSTORDER") { cells.push_back(new CellGhost); }
        else { cells.push_back(new CellO2GhostCartesian); }
        m_elements.push_back(new ElementCartesian());
        cells[i]->setElement(m_elements[i], i);
        cells[i]->pushBackSlope();
    }

    //Connectivity table
    //The number of the ghost cells is located apart from the number of "true" cells, hence start at m_numberCellsCalcul
    countParallelCells = m_numberCellsCalcul;
    int neighbourCpuCoordX, neighbourCpuCoordY;
    //In the x-direction
    if (m_CpuCoordX > 0)
    {
      neighbourCpuCoordX = m_CpuCoordX - 1;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      parallel.setNeighbour(neighbour);
      ix = 0;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
    if (m_CpuCoordX < m_numberCpuX - 1)
    {
      neighbourCpuCoordX = m_CpuCoordX + 1;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      parallel.setNeighbour(neighbour);
      ix = m_numberCellsX - 1;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
    //In the y-direction
    if (m_CpuCoordY > 0)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY - 1;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      parallel.setNeighbour(neighbour);
      iy = 0;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
    if (m_CpuCoordY < m_numberCpuY - 1)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY + 1;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      parallel.setNeighbour(neighbour);
      iy = m_numberCellsY - 1;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
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
    CpuFactorizationMin = 10.*std::fabs(static_cast<double>(m_numberCellsXGlobal) / tempNumberCpuX - static_cast<double>(m_numberCellsYGlobal) / tempNumberCpuY)
      + 10.*std::fabs(static_cast<double>(m_numberCellsYGlobal) / tempNumberCpuY - static_cast<double>(m_numberCellsZGlobal) / tempNumberCpuZ);

    //Optimization of the initial processor topology
    for (int i = 1; i <= Ncpu; i++) {
      if ((Ncpu % i == 0) && (m_numberCellsXGlobal / i >= minCellNumberOneDirection)) {
        for (int j = 1; j <= Ncpu/i; j++) {
          if ((Ncpu/i % j == 0) && (m_numberCellsYGlobal / j >= minCellNumberOneDirection)) {
            tempNumberCpuX = i;
            tempNumberCpuY = j;
            tempNumberCpuZ = Ncpu / (i*j);
            if (CpuFactorizationMin >= std::fabs(m_numberCellsXGlobal / tempNumberCpuX - m_numberCellsYGlobal / tempNumberCpuY) + std::fabs(m_numberCellsYGlobal / tempNumberCpuY - m_numberCellsZGlobal / tempNumberCpuZ)
              && m_numberCellsZGlobal / tempNumberCpuZ >= minCellNumberOneDirection) {
              m_numberCpuX = i;
              m_numberCpuY = j;
              m_numberCpuZ = Ncpu / (i*j);
              CpuFactorizationMin = std::fabs(m_numberCellsXGlobal / tempNumberCpuX - m_numberCellsYGlobal / tempNumberCpuY) + std::fabs(m_numberCellsYGlobal / tempNumberCpuY - m_numberCellsZGlobal / tempNumberCpuZ);
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

    //Generating cells 
    for (int i = 0; i < m_numberCellsCalcul; i++) {
        if (ordreCalcul == "FIRSTORDER") { cells.push_back(new Cell); }
        else { cells.push_back(new CellO2Cartesian); }
        m_elements.push_back(new ElementCartesian());
        cells[i]->setElement(m_elements[i], i);
    }
    for (int i = m_numberCellsCalcul; i < m_numberCellsTotal; i++) {
        if (ordreCalcul == "FIRSTORDER") { cells.push_back(new CellGhost); }
        else { cells.push_back(new CellO2GhostCartesian); }
        m_elements.push_back(new ElementCartesian());
        cells[i]->setElement(m_elements[i], i);
        cells[i]->pushBackSlope();
    }

    //Connectivity table
    //The number of the ghost cells is located apart from the number of "true" cells, hence start at m_numberCellsCalcul
    countParallelCells = m_numberCellsCalcul;
    int neighbourCpuCoordX, neighbourCpuCoordY, neighbourCpuCoordZ;
    //In the x-direction
    if (m_CpuCoordX > 0)
    {
      neighbourCpuCoordX = m_CpuCoordX - 1;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbourCpuCoordZ = m_CpuCoordZ;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      parallel.setNeighbour(neighbour);
      ix = 0;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
    if (m_CpuCoordX < m_numberCpuX - 1)
    {
      neighbourCpuCoordX = m_CpuCoordX + 1;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbourCpuCoordZ = m_CpuCoordZ;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      parallel.setNeighbour(neighbour);
      ix = m_numberCellsX - 1;
      for (iy = 0; iy < m_numberCellsY; iy++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
    //In the y-direction
    if (m_CpuCoordY > 0)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY - 1;
      neighbourCpuCoordZ = m_CpuCoordZ;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      parallel.setNeighbour(neighbour);
      iy = 0;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
    if (m_CpuCoordY < m_numberCpuY - 1)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY + 1;
      neighbourCpuCoordZ = m_CpuCoordZ;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      parallel.setNeighbour(neighbour);
      iy = m_numberCellsY - 1;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iz = 0; iz < m_numberCellsZ; iz++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
    //In the z-direction
    if (m_CpuCoordZ > 0)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbourCpuCoordZ = m_CpuCoordZ - 1;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      parallel.setNeighbour(neighbour);
      iz = 0;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iy = 0; iy < m_numberCellsY; iy++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
    if (m_CpuCoordZ < m_numberCpuZ - 1)
    {
      neighbourCpuCoordX = m_CpuCoordX;
      neighbourCpuCoordY = m_CpuCoordY;
      neighbourCpuCoordZ = m_CpuCoordZ + 1;
      neighbour = 0;
      neighbour += neighbourCpuCoordX;
      neighbour += neighbourCpuCoordY * m_numberCpuX;
      neighbour += neighbourCpuCoordZ * m_numberCpuX*m_numberCpuY;
      parallel.setNeighbour(neighbour);
      iz = m_numberCellsZ - 1;
      for (ix = 0; ix < m_numberCellsX; ix++)
      {
        for (iy = 0; iy < m_numberCellsY; iy++)
        {
          this->construitIGlobal(ix, iy, iz, iCell);
          parallel.addElementToSend(neighbour, cells[iCell]);
          parallel.addElementToReceive(neighbour, cells[countParallelCells]);
          parallel.addSlopesToSend(neighbour);
          parallel.addSlopesToReceive(neighbour);
          cells[countParallelCells]->setRankOfNeighborCPU(neighbour);
          ++countParallelCells;
        }
      }
    }
  }
}

//***********************************************************************

std::string MeshCartesian::whoAmI() const
{
  return "CARTESIAN";
}

//***********************************************************************

void MeshCartesian::setImmersedBoundaries(TypeMeshContainer<CellInterface*>* cellInterfacesLvl, std::string ordreCalcul) const
{
  for (unsigned int i = 0; i < cellInterfacesLvl[0].size(); i++) { 
    // N'est pas une CL
    if(cellInterfacesLvl[0][i]->whoAmI()==0){

      //2 murs
      if(cellInterfacesLvl[0][i]->getCellLeft()->getWall() && cellInterfacesLvl[0][i]->getCellRight()->getWall()){
        int b(0);
        for(b=0; b < cellInterfacesLvl[0][i]->getCellRight()->getCellInterfacesSize(); b++){
          if(cellInterfacesLvl[0][i]->getCellRight()->getCellInterface(b) == cellInterfacesLvl[0][i]) break;
        }
        cellInterfacesLvl[0][i]->getCellRight()->deleteInterface(b);
        for(b=0; b < cellInterfacesLvl[0][i]->getCellLeft()->getCellInterfacesSize(); b++){
          if(cellInterfacesLvl[0][i]->getCellLeft()->getCellInterface(b) == cellInterfacesLvl[0][i]) break;
        }
        cellInterfacesLvl[0][i]->getCellLeft()->deleteInterface(b);
        //Traitement particulier pour remove les communications si un des murs apparait en bord de domaine communicant
        if(cellInterfacesLvl[0][i]->getCellLeft()->isCellGhost()){
          parallel.deleteSlopesToSend(cellInterfacesLvl[0][i]->getCellLeft()->getRankOfNeighborCPU());
          parallel.deleteSlopesToReceive(cellInterfacesLvl[0][i]->getCellLeft()->getRankOfNeighborCPU());
          cellInterfacesLvl[0][i]->getCellLeft()->popBackSlope();
        }
        if(cellInterfacesLvl[0][i]->getCellRight()->isCellGhost()){
          parallel.deleteSlopesToSend(cellInterfacesLvl[0][i]->getCellRight()->getRankOfNeighborCPU());
          parallel.deleteSlopesToReceive(cellInterfacesLvl[0][i]->getCellRight()->getRankOfNeighborCPU());
          cellInterfacesLvl[0][i]->getCellRight()->popBackSlope();
        }
        delete cellInterfacesLvl[0][i];
        cellInterfacesLvl[0].erase(cellInterfacesLvl[0].begin()+i);  i--;
      }

      //mur a droite
      else if(cellInterfacesLvl[0][i]->getCellRight()->getWall()){
        Cell* buffCell= cellInterfacesLvl[0][i]->getCellLeft();
        Face* newFace = new FaceCartesian(*(static_cast<FaceCartesian*>(cellInterfacesLvl[0][i]->getFace())));
        int b(0);
        for(b=0; b < cellInterfacesLvl[0][i]->getCellRight()->getCellInterfacesSize(); b++){
          if(cellInterfacesLvl[0][i]->getCellRight()->getCellInterface(b) == cellInterfacesLvl[0][i]) break;
        }
        cellInterfacesLvl[0][i]->getCellRight()->deleteInterface(b);

        for(b=0; b < buffCell->getCellInterfacesSize(); b++){
          if(buffCell->getCellInterface(b) == cellInterfacesLvl[0][i]) break;
        }

        if(cellInterfacesLvl[0][i]->getCellRight()->isCellGhost()){
          parallel.deleteSlopesToSend(cellInterfacesLvl[0][i]->getCellRight()->getRankOfNeighborCPU());
          parallel.deleteSlopesToReceive(cellInterfacesLvl[0][i]->getCellRight()->getRankOfNeighborCPU());
          cellInterfacesLvl[0][i]->getCellRight()->popBackSlope();
        }

        if(cellInterfacesLvl[0][i]->getCellLeft()->isCellGhost()){
          parallel.deleteSlopesToSend(cellInterfacesLvl[0][i]->getCellLeft()->getRankOfNeighborCPU());
          parallel.deleteSlopesToReceive(cellInterfacesLvl[0][i]->getCellLeft()->getRankOfNeighborCPU());
          cellInterfacesLvl[0][i]->getCellLeft()->popBackSlope();
          cellInterfacesLvl[0][i]->getCellLeft()->deleteInterface(b);
          delete cellInterfacesLvl[0][i];
          cellInterfacesLvl[0].erase(cellInterfacesLvl[0].begin()+i);  i--;
        }
        else{
          delete cellInterfacesLvl[0][i];
          if (ordreCalcul == "FIRSTORDER") { cellInterfacesLvl[0][i] = new BoundCondWall(-1); }
          else { cellInterfacesLvl[0][i] = new BoundCondWallO2Cartesian(-1); }
          cellInterfacesLvl[0][i]->setFace(newFace);
          cellInterfacesLvl[0][i]->initialize(buffCell,buffCell);
          buffCell->setCellInterface(b,cellInterfacesLvl[0][i]);
          int allocateSlopeLocal = 0;
          cellInterfacesLvl[0][i]->allocateSlopes(allocateSlopeLocal);
        }
      }

      //mur a gauche
      else if(cellInterfacesLvl[0][i]->getCellLeft()->getWall()){
        Cell* buffCell= cellInterfacesLvl[0][i]->getCellRight();
        Face* newFace = new FaceCartesian(*(static_cast<FaceCartesian*>(cellInterfacesLvl[0][i]->getFace())));
        int b(0);
        for(b=0; b < cellInterfacesLvl[0][i]->getCellLeft()->getCellInterfacesSize(); b++){
          if(cellInterfacesLvl[0][i]->getCellLeft()->getCellInterface(b) == cellInterfacesLvl[0][i]) break;
        }
        cellInterfacesLvl[0][i]->getCellLeft()->deleteInterface(b);
        for(b=0; b < buffCell->getCellInterfacesSize(); b++){
          if(buffCell->getCellInterface(b) == cellInterfacesLvl[0][i]) break;
        }

        if(cellInterfacesLvl[0][i]->getCellLeft()->isCellGhost()){
          parallel.deleteSlopesToSend(cellInterfacesLvl[0][i]->getCellLeft()->getRankOfNeighborCPU());
          parallel.deleteSlopesToReceive(cellInterfacesLvl[0][i]->getCellLeft()->getRankOfNeighborCPU());
          cellInterfacesLvl[0][i]->getCellLeft()->popBackSlope();
        }

        if(cellInterfacesLvl[0][i]->getCellRight()->isCellGhost()){
          parallel.deleteSlopesToSend(cellInterfacesLvl[0][i]->getCellRight()->getRankOfNeighborCPU());
          parallel.deleteSlopesToReceive(cellInterfacesLvl[0][i]->getCellRight()->getRankOfNeighborCPU());
          cellInterfacesLvl[0][i]->getCellRight()->popBackSlope();
          cellInterfacesLvl[0][i]->getCellRight()->deleteInterface(b);
          delete cellInterfacesLvl[0][i];
          cellInterfacesLvl[0].erase(cellInterfacesLvl[0].begin()+i);  i--;
        }
        else{
          delete cellInterfacesLvl[0][i];
          if (ordreCalcul == "FIRSTORDER") { cellInterfacesLvl[0][i] = new BoundCondWall(-1); }
          else { cellInterfacesLvl[0][i] = new BoundCondWallO2Cartesian(-1); }
          cellInterfacesLvl[0][i]->setFace(newFace);
          cellInterfacesLvl[0][i]->getFace()->setNormal(-newFace->getNormal().getX(),
            -newFace->getNormal().getY(), -newFace->getNormal().getZ());
          cellInterfacesLvl[0][i]->getFace()->setTangent(-newFace->getTangent().getX(),
            -newFace->getTangent().getY(), -newFace->getTangent().getZ());
          cellInterfacesLvl[0][i]->initialize(buffCell,buffCell);
          buffCell->setCellInterface(b,cellInterfacesLvl[0][i]);
          int allocateSlopeLocal = 0;
          cellInterfacesLvl[0][i]->allocateSlopes(allocateSlopeLocal);
        }
      }
    }
    
    //mur a gauche d'une CL
    else{
      if(cellInterfacesLvl[0][i]->getCellLeft()->getWall()){
        int b(0);
        for(b=0; b < cellInterfacesLvl[0][i]->getCellLeft()->getCellInterfacesSize(); b++){
          if(cellInterfacesLvl[0][i]->getCellLeft()->getCellInterface(b) == cellInterfacesLvl[0][i]) break;
        }
        cellInterfacesLvl[0][i]->getCellLeft()->deleteInterface(b);
        delete cellInterfacesLvl[0][i];
        cellInterfacesLvl[0].erase(cellInterfacesLvl[0].begin()+i);  i--;
      }
    }
  }
}

//**************************************************************************
//******************************** PRINTING ********************************
//**************************************************************************

//**************************************************************************

std::string MeshCartesian::getStringExtent(bool global) const
{
  std::stringstream chaineExtent;

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

void MeshCartesian::getCoord(std::vector<double>& dataset, Axis axis) const
{
  switch (axis) {
  case X:
    for (int i = 0; i < m_numberCellsX; i++)
    {
      dataset.push_back(m_posXi[m_offsetX + i] - 0.5*m_dXi[m_offsetX + i]);
    }
    dataset.push_back(m_posXi[m_offsetX + m_numberCellsX - 1] + 0.5*m_dXi[m_offsetX + m_numberCellsX - 1]);
    break;
  case Y:
    for (int j = 0; j < m_numberCellsY; j++)
    {
      dataset.push_back(m_posYj[m_offsetY + j] - 0.5*m_dYj[m_offsetY + j]);
    }
    dataset.push_back(m_posYj[m_offsetY + m_numberCellsY - 1] + 0.5*m_dYj[m_offsetY + m_numberCellsY - 1]);
    break;
  case Z:
    for (int k = 0; k < m_numberCellsZ; k++)
    {
      dataset.push_back(m_posZk[m_offsetZ + k] - 0.5*m_dZk[m_offsetZ + k]);
    }
    dataset.push_back(m_posZk[m_offsetZ + m_numberCellsZ - 1] + 0.5*m_dZk[m_offsetZ + m_numberCellsZ - 1]);
    break;
  }
}

//****************************************************************************

void MeshCartesian::getData(TypeMeshContainer<Cell*>* cellsLvl, std::vector<double>& dataset, const int var, int phase) const
{
  dataset.clear();
  int numCell;
  double transport(0.);
  for (int k = 0; k < m_numberCellsZ; k++) {
    for (int j = 0; j < m_numberCellsY; j++) {
      for (int i = 0; i < m_numberCellsX; i++) {
        construitIGlobal(i, j, k, numCell);
        if (var > 0) { //We want to get the scalar data
          if (phase >= 0) { dataset.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnScalar(var)); }     //data de phases
          else if (phase == -1) { dataset.push_back(cellsLvl[0][numCell]->getMixture()->returnScalar(var)); }  //data de mixture
          else if (phase == -2) {
            transport = cellsLvl[0][numCell]->getTransport(var - 1).getValue();
            if (transport < 1.e-20) { transport = 0.; }
            dataset.push_back(transport);
          }
          else if (phase == -3) { dataset.push_back(cellsLvl[0][numCell]->getXi()); }
          else if (phase == -4) { dataset.push_back(cellsLvl[0][numCell]->getDensityGradient()); }
          else if (phase == -7) { // Saturation pressure
            dataset.push_back(cellsLvl[0][numCell]->getPsat());
          }
          else { Errors::errorMessage("MeshCartesian::getData: unknown number of phase: ", phase); }
        }
        else { //We want to get the vector data
          if (phase >= 0) { //data de phases
            dataset.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getX());
            dataset.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getY());
            dataset.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getZ());
          }
          else if (phase == -1) {  //data de mixture
            dataset.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getX());
            dataset.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getY());
            dataset.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getZ());
          }
          else { Errors::errorMessage("MeshCartesian::getData: unknown number of phase: ", phase); }
        } //End vector
      } // End X
    } //End Y
  } //End Z
}

//****************************************************************************

void MeshCartesian::setDataSet(std::vector<double>& dataset, TypeMeshContainer<Cell*>* cellsLvl, const int var, int phase) const
{
  int iterDataSet(0);
  Coord vec;
  for (unsigned int i = 0; i < cellsLvl[0].size(); i++) {
    if (var > 0) { //Scalars data are first set
      if (phase >= 0) { cellsLvl[0][i]->getPhase(phase)->setScalar(var, dataset[iterDataSet++]); } //phases data
      else if (phase == -1) { cellsLvl[0][i]->getMixture()->setScalar(var, dataset[iterDataSet++]); }  //mixture data
      else if (phase == -2) { cellsLvl[0][i]->getTransport(var - 1).setValue(dataset[iterDataSet++]); } //transport data
      else if (phase == -3) { cellsLvl[0][i]->setXi(dataset[iterDataSet++]); } //xi indicator
      else { Errors::errorMessage("MeshCartesian::setDataSet: unknown phase number: ", phase); }
    }
    else { //We want to get the vector data
      if (phase >= 0) { //Phases data
        vec.setXYZ(dataset[iterDataSet], dataset[iterDataSet + 1], dataset[iterDataSet + 2]);
        cellsLvl[0][i]->getPhase(phase)->setVector(-var, vec);
        iterDataSet += 3;
      }
      else if (phase == -1) {  //Mixture data
        vec.setXYZ(dataset[iterDataSet], dataset[iterDataSet + 1], dataset[iterDataSet + 2]);
        cellsLvl[0][i]->getMixture()->setVector(-var, vec);
        iterDataSet += 3;
      }
      else { Errors::errorMessage("MeshCartesian::setDataSet: unknown phase number: ", phase); }
    } //End vector
  }
}

//***********************************************************************
