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

//! \file      MUSGmshV2.h
//! \author    J. Caze
//! \version   1.0
//! \date      November 19 2019

#include "MUSGmshV2.h"

//***********************************************************************

MUSGmshV2::MUSGmshV2(const std::string& meshFile, const std::string& meshExtension, bool switchTags) : MUSGmsh(meshFile, meshExtension), m_switchTags(switchTags)
{}

//***********************************************************************

MUSGmshV2::~MUSGmshV2() {}

//***********************************************************************

void MUSGmshV2::initGeometryMonoCPU(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<CellInterface *> &cellInterfaces, std::string computeOrder)
{
  try {
    // 1) Reading nodes and elements
    // -----------------------------
    std::vector<ElementNS*>* neighborNodes; // Size number of nodes 
    this->readMeshMonoCPU(&neighborNodes);  // Fill m_nodes and m_elements

	// Display mesh information
	this->writeMeshInfoData();

    // CAUTION: Ordering of m_elements is important. Faces first, then cells.

    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << " D) BUILDING GEOMETRY ..." << std::endl;

    // 2) Assignment of cells to their geometric element
    // -------------------------------------------------
    // Cells and boundary counting
    if (m_numberElements3D == 0 && m_numberElements2D == 0) // 1D case
    {
      m_numberCellsCalcul = m_numberElements1D;
      m_numberBoundFaces = m_numberElements0D;
      m_geometrie = 1;
    }
    else if (m_numberElements3D == 0) // 2D case
    {
      m_numberCellsCalcul = m_numberElements2D;
      m_numberBoundFaces = m_numberElements1D;
      m_geometrie = 2;
    }
    else // 3D case
    {
      m_numberCellsCalcul = m_numberElements3D;
      m_numberBoundFaces = m_numberElements2D;
      m_geometrie = 3;
    }
    m_numberCellsTotal = m_numberCellsCalcul;

    // Assignment of non-reflecting boundary for missing limits
    unsigned int nbLimits(0);
    for (int i = 0; i < m_numberBoundFaces; i++) {
      unsigned int appPhys(m_elements[i]->getAppartenancePhysique());
      if (appPhys > nbLimits) { nbLimits = appPhys; }
    }
    for (unsigned int i = m_bound.size(); i < nbLimits; i++) {
      m_bound.push_back(new BoundCondNonReflecting);
    }

    // Assignment of cells to elements and counting inner faces
    m_numberInnerFaces = 0;
    for (int i = 0; i < m_numberCellsCalcul; i++)
    {
      if (computeOrder == "FIRSTORDER") { cells.push_back(new Cell); }
      else { cells.push_back(new CellO2); }
      cells[i]->setElement(m_elements[i + m_numberBoundFaces], i);
      m_numberInnerFaces += m_elements[i + m_numberBoundFaces]->getNumberFaces();
    }
    //    
    // -----------Adding neighbor to cells for 2nd order----------
    std::vector<ElementNS *>* neighbors; // Temporary vector of neighbors
    neighbors = new std::vector<ElementNS *>[m_numberCellsCalcul];
    for (int i = 0; i < m_numberCellsCalcul; i++)
    {
      ElementNS *e(m_elements[i + m_numberBoundFaces]); // Select element
      //std::cout << "el " << i + m_numberBoundFaces << " : " ;
      //std::vector<ElementNS *> neighbors; // Temporary vector of neighbors
      //1) Building vector of neighbors
      for (int n = 0; n < e->getNumberNoeuds(); n++) { // Node loop of element e
        int numNoeud = e->getNumNoeud(n);
        for (unsigned int v = 0; v < neighborNodes[numNoeud].size(); v++) { // Loop neighborNodes of node n
          bool add(true);
          if (neighborNodes[numNoeud][v]->getIndex() == e->getIndex()) add = false;
          //else if(neighborNodes[numNoeud][v]->getIndex() < m_numberBoundFaces) add = false;
          else {
            for (unsigned int vo = 0; vo < neighbors[i].size(); vo++) {
              if (neighborNodes[numNoeud][v]->getIndex() == neighbors[i][vo]->getIndex()) {
                add = false; break;
              }
            }
          }
          if (add) {
            neighbors[i].push_back(neighborNodes[numNoeud][v]);
            //std::cout << neighborNodes[numNoeud][v]->getIndex() << " ";
          }
        }
      }
    }
    // -----------------------------------------------------------

    m_numberInnerFaces -= m_numberBoundFaces; // We remove the boundaries
    m_numberInnerFaces /= 2; // The internal faces are all counted twice => we restore the truth
    m_numberFacesTotal = m_numberInnerFaces + m_numberBoundFaces;

    double volTot(0.);
    for (int i = 0; i < m_numberCellsCalcul; i++)
    {
      volTot += cells[i]->getElement()->getVolume();
    }
    std::cout << "Total volume : " << volTot << std::endl;

    // 3) Building connectivity table
    // -------------------------------
    // Sizing faces array
    m_faces = new FaceNS*[m_numberFacesTotal];
    int **facesBuff; int *sumNodesBuff; // We build a temporary array of faces to speed up the search process
    facesBuff = new int*[m_numberFacesTotal + 1];
    sumNodesBuff = new int[m_numberFacesTotal + 1];
    // Determination of the maximum number of nodes for the faces 
    int sizeFace; // Will be initialize at maximale size
    if (m_numberElements3D != 0)
    {
      if (m_numberQuadrangles != 0) { sizeFace = 4; }
      else if (m_numberTriangles != 0) { sizeFace = 3; }
      else { Errors::errorMessage("Issue in initGeometryMonoCPU for initialization of facesBuff array"); }
    }
    else if (m_numberElements2D != 0) { sizeFace = 2; }
    else { sizeFace = 1; }
    for (int i = 0; i < m_numberFacesTotal + 1; i++)
    {
      facesBuff[i] = new int[sizeFace];
    }

    // Inner faces
    int indexMaxFaces(0);
    clock_t tTemp(clock());
    float t1(0.);
    std::cout << "  1/Building faces ..." << std::endl;
    int frequenceImpressure(std::max((m_numberElements - m_numberBoundFaces) / 10, 1));
    for (int i = m_numberBoundFaces; i < m_numberElements; i++)
    {
      if ((i - m_numberBoundFaces) % frequenceImpressure == 0) { std::cout << "    " << (100 * (i - m_numberBoundFaces) / (m_numberElements - m_numberBoundFaces)) << "% ... " << std::endl; }
      m_elements[i]->construitFaces(m_nodes, m_faces, indexMaxFaces, facesBuff, sumNodesBuff);
    }
    for (int i = 0; i < m_numberFacesTotal + 1; i++) { delete facesBuff[i]; }
    delete[] facesBuff; delete[] sumNodesBuff;
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    std::cout << "    OK in " << t1 << " seconds" << std::endl;

    // Boundaries
    std::cout << "  2/Boundary elements attribution to boundary faces ..." << std::endl;
    tTemp = clock();
    frequenceImpressure = std::max(m_numberBoundFaces / 10, 1);
    for (int i = 0; i < m_numberBoundFaces; i++)
    {
      if (i%frequenceImpressure == 0)
      {
        std::cout << "    " << (100 * i / m_numberBoundFaces) << "% ... " << std::endl; 
      }
      // Assigning the boundary
      m_elements[i]->attributFaceLimite(m_nodes, m_faces, indexMaxFaces);
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    std::cout << "    OK in " << t1 << " seconds" << std::endl;

    // Link Geometry/cellInterfaces of compute
    std::cout << "  3/Linking Geometries -> Physics ..." << std::endl;
    tTemp = clock();
    int iMailleG, iMailleD;
    for (int i = 0; i < m_numberFacesTotal; i++)
    {
      if (m_faces[i]->getEstLimite()) // Physical boundary faces
      {
        int appPhys(m_faces[i]->getElementDroite()->getAppartenancePhysique() - 1); // (appartenance - 1) for array starting at zero
        if (appPhys >= static_cast<int>(m_bound.size()) || appPhys < 0) { Errors::errorMessage("Number of boundary conditions not suited"); }
        m_bound[appPhys]->creeLimite(cellInterfaces);
        cellInterfaces[i]->setFace(m_faces[i]);
        iMailleG = m_faces[i]->getElementGauche()->getIndex() - m_numberBoundFaces;
        iMailleD = iMailleG;
      }
      // Inner faces of domain
      else
      {
        if (computeOrder == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
        else { cellInterfaces.push_back(new CellInterfaceO2); }
        cellInterfaces[i]->setFace(m_faces[i]);
        iMailleG = m_faces[i]->getElementGauche()->getIndex() - m_numberBoundFaces;
        iMailleD = m_faces[i]->getElementDroite()->getIndex() - m_numberBoundFaces;
      }
      cellInterfaces[i]->initialize(cells[iMailleG], cells[iMailleD]);
      cells[iMailleG]->addCellInterface(cellInterfaces[i]);
      cells[iMailleD]->addCellInterface(cellInterfaces[i]);

      //-------------Ici il faut trouver les mailles voisines pour le compute ordre 2 multislopes ------------
      if (computeOrder != "FIRSTORDER") {

        //GAUCHE) Recherche des mailles a gauche
        //**************************************
        //ElementNS *eG(m_faces[i]->getElementGauche()); //Element a gauche
        //int index(eG->getIndex() - m_numberBoundFaces);


        //MeshUnStruct::rechercheElementsArrieres(eG,m_faces[i], cellInterfaces[i], neighbors[index], *cells);
        //MeshUnStruct::rechercheElementsAvants(eG, m_faces[i], cellInterfaces[i], neighbors[index], *cells);
        //
        //
        ////A enlever a terme...
        //double cos(0.),a,b,beta1,beta2;
        //ElementNS *eT;
        //Coord vG = eG->vecteur(m_faces[i]);
        //Coord MB1, MB2, MB;



        //Coord vG = eG->vecteur(m_faces[i]);
        //ElementNS *eT; //Element a tester
        ////1) Recherche du premier element arriere gauche BG1M
        ////---------------------------------------------------
        //double cos(0.);
        //for (unsigned int v = 0; v < neighbors[index].size(); v++) {
        //  eT = neighbors[index][v];
        //  Coord vT = eT->vecteur(eG);
        //  double cosTemp(Coord::cos(vT, vG));
        //  if (cosTemp >= cos) {
        //    cos = cosTemp;
        //    if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG1M a NULL
        //      cellInterfaces[i]->setB(BG1M, 0);
        //    }
        //    else { //Sinon on met a jour avec la maille la plus loin
        //      Cell *c(cells[eT->getNumCellAssociee()]);
        //      cellInterfaces[i]->setB(BG1M, c);
        //    }
        //  }
        //}
        ////2) Recherche du second element arriere gauche BG2M
        ////---------------------------------------------------
        //cos = 0.;
        //if (cellInterfaces[i]->getB(BG1M) != 0) {  //Cas non cellInterface
        //  Element *e1(cellInterfaces[i]->getB(BG1M)->getElement());
        //  Coord v1 = e1->vecteur(eG);
        //  Coord sin1(Coord::sin(v1, vG));
        //  if (std::fabs(sin1.norm()) <= 1e-8) {
        //    cellInterfaces[i]->setB(BG2M, cellInterfaces[i]->getB(BG1M)); // Si le cos est nul, meme cell pour B1 et B2
        //  }
        //  else {
        //    for (unsigned int v = 0; v < neighbors[index].size(); v++) {
        //      eT = neighbors[index][v];
        //      if (eT == e1) continue; // voisin suivant si voisin deja utilise 
        //      Coord vT = eT->vecteur(eG);
        //      Coord sinTemp(Coord::sin(vT, vG));

        //      if (sinTemp.scalar(sin1) <= 0.) {
        //        double cosTemp(Coord::cos(vT, vG));
        //        if (cosTemp >= cos) {
        //          cos = cosTemp;
        //          if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG1M a NULL
        //            cellInterfaces[i]->setB(BG2M, 0);
        //          }
        //          else {  //Sinon on met a jour avec la 2 eme maille la plus loin
        //            Cell *c(cells[eT->getNumCellAssociee()]);
        //            cellInterfaces[i]->setB(BG2M, c);
        //          }
        //        }
        //      } //fin sin*sin <0
        //    }  //fin boucle neighbors
        //  } //fin if 
        //} //Fin if cellInterface

        ////Determination des ponderations arrieres
        //Coord MB1; if (cellInterfaces[i]->getB(BG1M) != 0) { MB1.setFromSubtractedVectors(m_faces[i]->getPos(), cellInterfaces[i]->getB(BG1M)->getPosition()); }
        //Coord MB2; if (cellInterfaces[i]->getB(BG2M) != 0) { MB2.setFromSubtractedVectors(m_faces[i]->getPos(), cellInterfaces[i]->getB(BG2M)->getPosition()); }
        //Coord MB; MB.setFromSubtractedVectors(m_faces[i]->getPos(), eG->getPosition());
        //double a, b, beta1(1.), beta2(0.);
        //a = (MB1.vectoriel(MB)).norm() / MB.norm();
        //b = (MB2.vectoriel(MB)).norm() / MB.norm();
        //if (std::fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
        //cellInterfaces[i]->setBeta(betaG1M, beta1);
        //cellInterfaces[i]->setBeta(betaG2M, beta2);
        ////Calcul de la distance du vertex M (cellInterface de maille) au vertex H (barycenter)
        //if (cellInterfaces[i]->getB(BG1M) != 0) {
        //  cos = Coord::cos(MB, MB1);
        //  double d, e, c;
        //  d = cos*MB1.norm();
        //  Coord B1B2; B1B2 = MB1 - MB2;
        //  e = B1B2.norm()*beta1;
        //  c = sqrt(e*e - a*a);
        //  d += c;
        //  cellInterfaces[i]->setDistanceH(distanceHGM, d);
        //}
        ////cout << eG->getIndex() << " : " << endl;
        ////if (cellInterfaces[i]->getB(BG1M) != 0) cout << "BG1M = " << cellInterfaces[i]->getB(BG1M)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHGM) << endl;
        ////if (cellInterfaces[i]->getB(BG2M) != 0) cout << "BG2M = " << cellInterfaces[i]->getB(BG2M)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHGM) << endl;

        ////3) Recherche du premier element avant gauche BG1P
        ////-------------------------------------------------
        //cos = 0.;
        //for (unsigned int v = 0; v < neighbors[index].size(); v++) {
        //  eT = neighbors[index][v];
        //  Coord vT = eT->vecteur(eG);
        //  double cosTemp(Coord::cos(vT, vG));
        //  if (cosTemp <= cos) {
        //    cos = cosTemp;
        //    if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG1P a NULL
        //      cellInterfaces[i]->setB(BG1P, 0);
        //    }
        //    else { //Sinon on met a jour avec la maille la plus loin
        //      Cell *c(cells[eT->getNumCellAssociee()]);
        //      cellInterfaces[i]->setB(BG1P, c);
        //    }
        //  }
        //}
        ////4) Recherche du second element avant gauche BG2P
        ////------------------------------------------------
        //cos = 0.;
        //if (cellInterfaces[i]->getB(BG1P) != 0) {  //Cas non cellInterface
        //  Element *e1(cellInterfaces[i]->getB(BG1P)->getElement());
        //  Coord v1 = e1->vecteur(eG);
        //  Coord sin1(Coord::sin(v1, vG));
        //  if (std::fabs(sin1.norm()) <= 1e-8) {
        //    cellInterfaces[i]->setB(BG2P, cellInterfaces[i]->getB(BG1P)); // Si le cos est nul, meme cell pour B1 et B2
        //  }
        //  else {
        //    for (unsigned int v = 0; v < neighbors[index].size(); v++) {
        //      eT = neighbors[index][v];
        //      if (eT == e1) continue; // voisin suivant si voisin deja utilise 
        //      Coord vT = eT->vecteur(eG);
        //      Coord sinTemp(Coord::sin(vT, vG));

        //      if (sinTemp.scalar(sin1) <= 0.) {
        //        double cosTemp(Coord::cos(vT, vG));
        //        if (cosTemp <= cos) {
        //          cos = cosTemp;
        //          if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG2P a NULL
        //            cellInterfaces[i]->setB(BG2P, 0);
        //          }
        //          else {  //Sinon on met a jour avec la 2 eme maille la plus loin
        //            Cell *c(cells[eT->getNumCellAssociee()]);
        //            cellInterfaces[i]->setB(BG2P, c);
        //          }
        //        }
        //      } //fin sin*sin <0
        //    }  //fin boucle neighbors
        //  } //fin if 
        //} //Fin if cellInterface

        ////Determination des ponderations avant
        //MB1 = 0.; if (cellInterfaces[i]->getB(BG1P) != 0) { MB1.setFromSubtractedVectors(cellInterfaces[i]->getB(BG1P)->getPosition(), m_faces[i]->getPos()); }
        //MB2 = 0.; if (cellInterfaces[i]->getB(BG2P) != 0) { MB2.setFromSubtractedVectors(cellInterfaces[i]->getB(BG2P)->getPosition(), m_faces[i]->getPos()); }
        //MB.setFromSubtractedVectors(eG->getPosition(), m_faces[i]->getPos());
        //beta1 = 1., beta2 = 0.;
        //a = (MB1.vectoriel(MB)).norm() / MB.norm();
        //b = (MB2.vectoriel(MB)).norm() / MB.norm();
        //if (std::fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
        //cellInterfaces[i]->setBeta(betaG1P, beta1);
        //cellInterfaces[i]->setBeta(betaG2P, beta2);
        ////Calcul de la distance du vertex M (cellInterface de maille) au vertex H (barycenter)
        //if (cellInterfaces[i]->getB(BG1P) != 0) {
        //  cos = Coord::cos(MB, -1.*MB1);
        //  double d, e, c;
        //  d = cos*MB1.norm();
        //  Coord B1B2; B1B2 = MB1 - MB2;
        //  e = B1B2.norm()*beta1;
        //  c = sqrt(e*e - a*a);
        //  d += c;
        //  cellInterfaces[i]->setDistanceH(distanceHGP, d);
        //}
        ////cout << eG->getIndex() << " : " << endl;
        ////if (cellInterfaces[i]->getB(BG1P) != 0) cout << "BG1P = " << cellInterfaces[i]->getB(BG1P)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHGP) << endl;
        ////if (cellInterfaces[i]->getB(BG2P) != 0) cout << "BG2P = " << cellInterfaces[i]->getB(BG2P)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHGP) << endl;


        ////DROITE) Recherche des mailles a droite
        ////**************************************
        //ElementNS *eD(m_faces[i]->getElementDroite()); //Element a droite
        //index = eD->getIndex() - m_numberBoundFaces;
        //if (index >= 0) {
        //  Coord vD = eD->vecteur(m_faces[i]);
        //  //1) Recherche du premier element arriere droite BD1M
        //  //---------------------------------------------------
        //  cos = 0.;
        //  for (unsigned int v = 0; v < neighbors[index].size(); v++) {
        //    eT = neighbors[index][v];
        //    Coord vT = eT->vecteur(eD);
        //    double cosTemp(Coord::cos(vT, vD));
        //    if (cosTemp >= cos) {
        //      cos = cosTemp;
        //      if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG1M a NULL
        //        cellInterfaces[i]->setB(BD1M, 0);
        //      }
        //      else { //Sinon on met a jour avec la maille la plus loin
        //        Cell *c(cells[eT->getNumCellAssociee()]);
        //        cellInterfaces[i]->setB(BD1M, c);
        //      }
        //    }
        //  }
        //  //2) Recherche du second element arriere droite BD2M
        //  //---------------------------------------------------
        //  cos = 0.;
        //  if (cellInterfaces[i]->getB(BD1M) != 0) {  //Cas non cellInterface
        //    Element *e1(cellInterfaces[i]->getB(BD1M)->getElement());
        //    Coord v1 = e1->vecteur(eD);
        //    Coord sin1(Coord::sin(v1, vD));
        //    if (std::fabs(sin1.norm()) <= 1e-8) {
        //      cellInterfaces[i]->setB(BD2M, cellInterfaces[i]->getB(BD1M)); // Si le cos est nul, meme cell pour B1 et B2
        //    }
        //    else {
        //      for (unsigned int v = 0; v < neighbors[index].size(); v++) {
        //        eT = neighbors[index][v];
        //        if (eT == e1) continue; // voisin suivant si voisin deja utilise 
        //        Coord vT = eT->vecteur(eD);
        //        Coord sinTemp(Coord::sin(vT, vD));

        //        if (sinTemp.scalar(sin1) <= 0.) {
        //          double cosTemp(Coord::cos(vT, vD));
        //          if (cosTemp >= cos) {
        //            cos = cosTemp;
        //            if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG1M a NULL
        //              cellInterfaces[i]->setB(BD2M, 0);
        //            }
        //            else {  //Sinon on met a jour avec la 2 eme maille la plus loin
        //              Cell *c(cells[eT->getNumCellAssociee()]);
        //              cellInterfaces[i]->setB(BD2M, c);
        //            }
        //          }
        //        } //fin sin*sin <0
        //      }  //fin boucle neighbors
        //    } //fin if 
        //  } //Fin if cellInterface

        //  //Determination des ponderations arrieres
        //  MB1 = 0.; if (cellInterfaces[i]->getB(BD1M) != 0) { MB1.setFromSubtractedVectors(m_faces[i]->getPos(), cellInterfaces[i]->getB(BD1M)->getPosition()); }
        //  MB2 = 0.; if (cellInterfaces[i]->getB(BD2M) != 0) { MB2.setFromSubtractedVectors(m_faces[i]->getPos(), cellInterfaces[i]->getB(BD2M)->getPosition()); }
        //  MB.setFromSubtractedVectors(m_faces[i]->getPos(), eD->getPosition());
        //  beta1 = 1.; beta2 = 0.;
        //  a = (MB1.vectoriel(MB)).norm() / MB.norm();
        //  b = (MB2.vectoriel(MB)).norm() / MB.norm();
        //  if (std::fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
        //  cellInterfaces[i]->setBeta(betaD1M, beta1);
        //  cellInterfaces[i]->setBeta(betaD2M, beta2);
        //  //Calcul de la distance du vertex M (cellInterface de maille) au vertex H (barycenter)
        //  if (cellInterfaces[i]->getB(BD1M) != 0) {
        //    cos = Coord::cos(MB, MB1);
        //    double d, e, c;
        //    d = cos*MB1.norm();
        //    Coord B1B2; B1B2 = MB1 - MB2;
        //    e = B1B2.norm()*beta1;
        //    c = sqrt(e*e - a*a);
        //    d += c;
        //    cellInterfaces[i]->setDistanceH(distanceHDM, d);
        //  }
        //  //cout << eD->getIndex() << " : " << endl;
        //  //if (cellInterfaces[i]->getB(BD1M) != 0) cout << "BD1M = " << cellInterfaces[i]->getB(BD1M)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHDM) << endl;
        //  //if (cellInterfaces[i]->getB(BD2M) != 0) cout << "BD2M = " << cellInterfaces[i]->getB(BD2M)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHDM) << endl;

        //  //3) Recherche du premier element avant droite BD1P
        //  //-------------------------------------------------
        //  cos = 0.;
        //  for (unsigned int v = 0; v < neighbors[index].size(); v++) {
        //    eT = neighbors[index][v];
        //    Coord vT = eT->vecteur(eD);
        //    double cosTemp(Coord::cos(vT, vD));
        //    if (cosTemp <= cos) {
        //      cos = cosTemp;
        //      if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG1M a NULL
        //        cellInterfaces[i]->setB(BD1P, 0);
        //      }
        //      else { //Sinon on met a jour avec la maille la plus loin
        //        Cell *c(cells[eT->getNumCellAssociee()]);
        //        cellInterfaces[i]->setB(BD1P, c);
        //      }
        //    }
        //  }
        //  //4) Recherche du second element avant droite BD2P
        //  //------------------------------------------------
        //  cos = 0.;
        //  if (cellInterfaces[i]->getB(BD1P) != 0) {  //Cas non cellInterface
        //    Element *e1(cellInterfaces[i]->getB(BD1P)->getElement());
        //    Coord v1 = e1->vecteur(eD);
        //    Coord sin1(Coord::sin(v1, vD));
        //    if (std::fabs(sin1.norm()) <= 1e-8) {
        //      cellInterfaces[i]->setB(BD2P, cellInterfaces[i]->getB(BD1P)); // Si le cos est nul, meme cell pour B1 et B2
        //    }
        //    else {
        //      for (unsigned int v = 0; v < neighbors[index].size(); v++) {
        //        eT = neighbors[index][v];
        //        if (eT == e1) continue; // voisin suivant si voisin deja utilise 
        //        Coord vT = eT->vecteur(eD);
        //        Coord sinTemp(Coord::sin(vT, vD));

        //        if (sinTemp.scalar(sin1) <= 0.) {
        //          double cosTemp(Coord::cos(vT, vD));
        //          if (cosTemp <= cos) {
        //            cos = cosTemp;
        //            if (eT->getIndex() < m_numberBoundFaces) { //Si c est une limite on remet BG1M a NULL
        //              cellInterfaces[i]->setB(BD2P, 0);
        //            }
        //            else {  //Sinon on met a jour avec la 2 eme maille la plus loin
        //              Cell *c(cells[eT->getNumCellAssociee()]);
        //              cellInterfaces[i]->setB(BD2P, c);
        //            }
        //          }
        //        } //fin sin*sin <0
        //      }  //fin boucle neighbors
        //    } //fin if 
        //  } //Fin if cellInterface
        //} //Fin si limite

        ////Determination des ponderations avant
        //MB1 = 0.; if (cellInterfaces[i]->getB(BD1P) != 0) { MB1.setFromSubtractedVectors(cellInterfaces[i]->getB(BD1P)->getPosition(), m_faces[i]->getPos()); }
        //MB2 = 0.; if (cellInterfaces[i]->getB(BD2P) != 0) { MB2.setFromSubtractedVectors(cellInterfaces[i]->getB(BD2P)->getPosition(), m_faces[i]->getPos()); }
        //MB.setFromSubtractedVectors(eD->getPosition(), m_faces[i]->getPos());
        //beta1 = 1., beta2 = 0.;
        //a = (MB1.vectoriel(MB)).norm() / MB.norm();
        //b = (MB2.vectoriel(MB)).norm() / MB.norm();
        //if (std::fabs(a + b) > 1e-8) { beta1 = b / (a + b); beta2 = 1. - beta1; }
        //cellInterfaces[i]->setBeta(betaD1P, beta1);
        //cellInterfaces[i]->setBeta(betaD2P, beta2);
        ////Calcul de la distance du vertex M (cellInterface de maille) au vertex H (barycenter)
        //if (cellInterfaces[i]->getB(BD1P) != 0) {
        //  cos = Coord::cos(MB, -1.*MB1);
        //  double d, e, c;
        //  d = cos*MB1.norm();
        //  Coord B1B2; B1B2 = MB1 - MB2;
        //  e = B1B2.norm()*beta1;
        //  c = sqrt(e*e - a*a);
        //  d += c;
        //  cellInterfaces[i]->setDistanceH(distanceHDP, d);
        //}
        //cout << eD->getIndex() << " : " << endl;
        //if (cellInterfaces[i]->getB(BD1P) != 0) cout << "BD1P = " << cellInterfaces[i]->getB(BD1P)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHDP) << endl;
        //if (cellInterfaces[i]->getB(BD2P) != 0) cout << "BD2P = " << cellInterfaces[i]->getB(BD2P)->getElement()->getIndex() << " " << cellInterfaces[i]->getDistanceH(distanceHDP) << endl;



        //cout << "******" << endl;
        //cout << eG->getIndex() << " " << eD->getIndex() << " : " << endl;
        //if (cellInterfaces[i]->getB(BG1M) != 0) cout << "BG1M = " << cellInterfaces[i]->getB(BG1M)->getElement()->getIndex() << endl;
        //if (cellInterfaces[i]->getB(BG2M) != 0) cout << "BG2M = " << cellInterfaces[i]->getB(BG2M)->getElement()->getIndex() << endl;
        //if (cellInterfaces[i]->getB(BG1P) != 0) cout << "BG1P = " << cellInterfaces[i]->getB(BG1P)->getElement()->getIndex() << endl;
        //if (cellInterfaces[i]->getB(BG2P) != 0) cout << "BG2P = " << cellInterfaces[i]->getB(BG2P)->getElement()->getIndex() << endl;
        //if (cellInterfaces[i]->getB(BD1M) != 0) cout << "BD1M = " << cellInterfaces[i]->getB(BD1M)->getElement()->getIndex() << endl;
        //if (cellInterfaces[i]->getB(BD2M) != 0) cout << "BD2M = " << cellInterfaces[i]->getB(BD2M)->getElement()->getIndex() << endl;
        //if (cellInterfaces[i]->getB(BD1P) != 0) cout << "BD1P = " << cellInterfaces[i]->getB(BD1P)->getElement()->getIndex() << endl;
        //if (cellInterfaces[i]->getB(BD2P) != 0) cout << "BD2P = " << cellInterfaces[i]->getB(BD2P)->getElement()->getIndex() << endl;

      } //Fin preparation neighbors ordre 2
      //---------------------------------------------------------------------------------------

    } // End face
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    std::cout << "    OK in " << t1 << " seconds" << std::endl;
    std::cout << "... BUILDING GEOMETRY COMPLETE " << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    delete[] neighbors;
    delete[] neighborNodes;
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MUSGmshV2::readMeshMonoCPU(std::vector<ElementNS*>** neighborNodes)
{
  // Re-open mesh file from now on 
  try {
    // 1) Opening mesh file at Gmsh format V2
    // --------------------------------------
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << " C) READING MESH FILE " + m_meshFile + " IN PROGRESS ..." << std::endl;
    std::ifstream meshFile(m_meshFile.c_str(), std::ios::in);
    if (!meshFile){ throw ErrorECOGEN("mesh file not found : " + m_meshFile, __FILE__, __LINE__); }
    std::string currentLine;

    // Skiping Gmsh version + go to nodes
    // -----------------------------------
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);

    // 2) Filling m_nodes array
    // ------------------------
    std::cout << "  1/Mesh nodes reading ...";
    meshFile >> m_numberNodes;
    meshFile.ignore(1, '\n');
    m_nodes = new Coord[m_numberNodes];
    *neighborNodes = new std::vector<ElementNS*>[m_numberNodes];
    int useless(0); double x, y, z;
    for (int i = 0; i < m_numberNodes; i++)
    {
      meshFile >> useless >> x >> y >> z;
      m_nodes[i].setXYZ(x, y, z);
    }
    meshFile.ignore(1, '\n');
    getline(meshFile, currentLine); // Skip keyword $EndNodes
    getline(meshFile, currentLine); // Skip keyword $Elements
    std::cout << "OK" << std::endl;

    // 3) 1D/2D/3D elements are stored in m_elements array / counting
    // --------------------------------------------------------------
	std::cout << "  2/0D/1D/2D/3D elements reading ..." << std::endl;
    meshFile >> m_numberElements;
    meshFile.ignore(1, '\n');
    // Allocation array of elements
    m_elements = new ElementNS*[m_numberElements];
    // Reading elements and assigning geometric properties 
    m_numberElements1D = 0, m_numberElements2D = 0, m_numberElements3D = 0;
    int posBeginElement = static_cast<int>(meshFile.tellg()); // id of the position of elements in file for quick return
    int nodeG;
    for (int i = 0; i < m_numberElements; i++) {
      this->readElement(m_nodes, meshFile, &m_elements[i]);
      if (m_elements[i]->getTypeGmsh() == 15) { m_numberElements0D++; }
      else if (m_elements[i]->getTypeGmsh() == 1) { m_numberElements1D++; }
      else if (m_elements[i]->getTypeGmsh() <= 3) { m_numberElements2D++; m_totalSurface += m_elements[i]->getVolume(); }
      else if (m_elements[i]->getTypeGmsh() <= 7) { m_numberElements3D++; m_totalVolume += m_elements[i]->getVolume(); }
      else { throw ErrorECOGEN("Element type of .msh not handled by ECOGEN", __FILE__, __LINE__); }
      // Assignment element i neighbor for concerned nodes (2nd order muiltislopes)
      for (int n = 0; n < m_elements[i]->getNumberNoeuds(); n++) {
        nodeG = m_elements[i]->getNumNoeud(n);
        (*neighborNodes)[nodeG].push_back(m_elements[i]);
      }
    }
    m_numberInnerElements = m_numberElements;
    
    //for (int n = 0; n < m_numberNodes; n++) {
    //  for (int e = 0; e < (*neighborNodes)[n].size(); e++) {
    //    std::cout << (*neighborNodes)[n][e]->getIndex() << " ";
    //  }
    //  std::cout << std::endl;
    //}
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MUSGmshV2::initGeometryParallel(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost, TypeMeshContainer<CellInterface *> &cellInterfaces, std::string computeOrder)
{
  //KS//FP//DEV// Distinction between cells and cellsGhost not done here. To do in the future.
  clock_t totalTime(clock());
  try {
    // 1) Reading nodes and elements
    // -----------------------------
    this->readMeshParallel(); // Filling m_nodes and m_elements
    if (rankCpu == 0)
    {
      std::cout << "------------------------------------------------------" << std::endl;
      std::cout << " D) BUILDING GEOMETRY ..." << std::endl;
    }

    // 2) Assignment of cells to their geometric element
    // -------------------------------------------------
    // Cells and boundary counting
    if (m_numberElements3D == 0 && m_numberElements2D == 0) // 1D case
    {
      m_numberCellsCalcul = m_numberElements1D;
      m_numberBoundFaces = m_numberElements0D;
      m_geometrie = 1;
    }
    else if (m_numberElements3D == 0) // 2D case
    {
      m_numberCellsCalcul = m_numberElements2D;
      m_numberBoundFaces = m_numberElements1D;
      m_geometrie = 2;
    }
    else // 3D case
    {
      m_numberCellsCalcul = m_numberElements3D;
      m_numberBoundFaces = m_numberElements2D;
      m_geometrie = 3;
    }

    // Sizing array of cells
    m_numberCellsTotal = m_numberCellsCalcul + m_numberGhostElements;

    // Assignment of non-reflecting boundary for missing limits
    unsigned int nbLimits(0);
    for (int i = 0; i < m_numberBoundFaces; i++) {
      unsigned int appPhys(m_elements[i]->getAppartenancePhysique());
      if (appPhys > nbLimits) { nbLimits = appPhys; }
    }
    for (unsigned int i = m_bound.size(); i < nbLimits; i++) {
      m_bound.push_back(new BoundCondNonReflecting);
    }

    // Assignment of cells to elements and counting inner faces
    m_numberInnerFaces = 0;
    for (int i = 0; i < m_numberCellsTotal; i++)
    {
      if (computeOrder == "FIRSTORDER") { cells.push_back(new Cell); }
      else { cells.push_back(new CellO2); }
      cells[i]->setElement(m_elements[i + m_numberBoundFaces], i);
      if (i < m_numberCellsCalcul) { m_numberInnerFaces += m_elements[i + m_numberBoundFaces]->getNumberFaces(); }
    }
    m_numberInnerFaces -= m_numberBoundFaces + m_numberFacesParallel; // We remove the boundaries and communicating faces
    m_numberInnerFaces /= 2; // The internal faces are all counted twice => we restore the truth
    m_numberFacesTotal = m_numberInnerFaces + m_numberBoundFaces + m_numberFacesParallel;

    // 3) Building internal connectivity table
    // ---------------------------------------
    
    // Sizing faces array
    m_faces = new FaceNS*[m_numberFacesTotal];
    int **facesBuff; int *sumNodesBuff; // We build a temporary array of faces to speed up the search process
    facesBuff = new int*[m_numberFacesTotal + 1];
    sumNodesBuff = new int[m_numberFacesTotal + 1];
    // Determination of the maximum number of nodes for the faces 
    int sizeFace; // Will be initialize at maximale size
    if (m_numberElements3D != 0)
    {
      if (m_numberQuadrangles != 0) { sizeFace = 4; }
      else if (m_numberTriangles != 0) { sizeFace = 3; }
      else { Errors::errorMessage("Issue in initGeometryMonoCPU for initialization of facesBuff array"); }
    }
    else if (m_numberElements2D != 0) { sizeFace = 2; }
    else { sizeFace = 1; } 
    for (int i = 0; i < m_numberFacesTotal + 1; i++) // +1 is used for the search for the existence of faces
    {
      facesBuff[i] = new int[sizeFace]; // Unknown on the number of points of a face (4 seems to be the max) 
    }

    // Inner faces
    // -----------
    int indexMaxFaces(0);
    MPI_Barrier(MPI_COMM_WORLD);
    int frequenceImpressure;
    clock_t tTemp(clock());
    float t1(0.);
    if (rankCpu == 0)
    {
      std::cout << "  1/Building faces ..." << std::endl;
      frequenceImpressure = std::max((m_numberInnerElements - m_numberBoundFaces) / 10, 1);
    }
    for (int i = m_numberBoundFaces; i < m_numberInnerElements; i++)
    {
      if (rankCpu == 0 && (i - m_numberBoundFaces) % frequenceImpressure == 0) { std::cout << "    " << (100 * (i - m_numberBoundFaces) / (m_numberInnerElements - m_numberBoundFaces)) << "% ... " << std::endl; }
      // Building
      m_elements[i]->construitFaces(m_nodes, m_faces, indexMaxFaces, facesBuff, sumNodesBuff);
    }
    for (int i = 0; i < m_numberFacesTotal + 1; i++) { delete facesBuff[i]; }
    delete[] facesBuff; delete[] sumNodesBuff;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      std::cout << "    OK in " << t1 << " seconds" << std::endl;
    }

    // Boundaries
    // ----------
    tTemp = clock();
    if (rankCpu == 0)
    {
      std::cout << "  2/Boundary elements attribution to boundary faces ..." << std::endl;
      frequenceImpressure = std::max(m_numberBoundFaces / 10, 1);
    }
    for (int i = 0; i < m_numberBoundFaces; i++)
    {
      if (rankCpu == 0 && i%frequenceImpressure == 0) 
      {
         std::cout << "    " << (100 * i / m_numberBoundFaces) << "% ... " << std::endl; 
      }
      // Assigning the boundary
      m_elements[i]->attributFaceLimite(m_nodes, m_faces, indexMaxFaces);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      std::cout << "    OK in " << t1 << " seconds" << std::endl;
    }

    // Communications
    // --------------
    tTemp = clock();
    if (rankCpu == 0)
    {
      std::cout << "  3/Ghost cells attribution to communicating faces ..." << std::endl;
      frequenceImpressure = std::max((m_numberElements - m_numberInnerElements) / 10, 1);
    }
    for (int i = m_numberInnerElements; i < m_numberElements; i++)
    {
      if (rankCpu == 0 && (i - m_numberInnerElements) % frequenceImpressure == 0) { std::cout << "    " << (100 * (i - m_numberInnerElements) / (m_numberElements - m_numberInnerElements)) << "% ... " << std::endl; }
      // Assigning missing boundary
      m_elements[i]->attributFaceCommunicante(m_nodes, m_faces, indexMaxFaces, m_numberInnerNodes);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      std::cout << "    OK in " << t1 << " seconds" << std::endl;
    }

    // Link Geometry/cellInterfaces of compute
    // ---------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    tTemp = clock();
    if (rankCpu == 0)
    {
      std::cout << "  4/Linking Geometries -> Physics ..." << std::endl;
      frequenceImpressure = std::max(m_numberFacesTotal / 10, 1);
    }
    int iMailleG, iMailleD;

    for (int i = 0; i < m_numberFacesTotal; i++)
    {
      if (rankCpu == 0 && i%frequenceImpressure == 0) { std::cout << "    " << 100 * i / m_numberFacesTotal << "% ... " << std::endl; }
      // Boundary faces (physical or communication boundary)
      if (m_faces[i]->getEstLimite())
      {
        // Communication boundary
        if (m_faces[i]->getEstComm())
        {
          if (computeOrder == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
          else { cellInterfaces.push_back(new CellInterfaceO2); }
          cellInterfaces[i]->setFace(m_faces[i]);
          iMailleG = m_faces[i]->getElementGauche()->getNumCellAssociee();
          iMailleD = m_faces[i]->getElementDroite()->getNumCellAssociee();
        }
        // Physical boundary
        else
        {
          int appPhys(m_faces[i]->getElementDroite()->getAppartenancePhysique() - 1); // (appartenance - 1) for array starting at zero
          if (appPhys >= static_cast<int>(m_bound.size()) || appPhys < 0) { throw ErrorECOGEN("Number of boundary conditions not suited", __FILE__, __LINE__); }
          m_bound[appPhys]->creeLimite(cellInterfaces);
          cellInterfaces[i]->setFace(m_faces[i]);
          iMailleG = m_faces[i]->getElementGauche()->getNumCellAssociee();
          iMailleD = iMailleG;
        }
      }
      // Inner faces of domain
      else
      {
        if (computeOrder == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
        else { cellInterfaces.push_back(new CellInterfaceO2); }
        cellInterfaces[i]->setFace(m_faces[i]);
        iMailleG = m_faces[i]->getElementGauche()->getNumCellAssociee();
        iMailleD = m_faces[i]->getElementDroite()->getNumCellAssociee();
      }
      cellInterfaces[i]->initialize(cells[iMailleG], cells[iMailleD]);
      cells[iMailleG]->addCellInterface(cellInterfaces[i]);
      cells[iMailleD]->addCellInterface(cellInterfaces[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      std::cout << "    OK in " << t1 << " seconds" << std::endl;
    }
    
    // 4) Building connectivity table for parallel CPUs
    // ------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    tTemp = clock();
    if (rankCpu == 0)
    {
      std::cout << "  5/Building connectivity tables for CPUs ..." << std::endl;
      frequenceImpressure = std::max((m_numberElements - m_numberBoundFaces) / 10, 1);
    }

    std::vector<int> numberElementsToSend(Ncpu);
    for (int i = 0; i < Ncpu; i++) { numberElementsToSend[i] = 0; }
    std::vector< std::vector<int> > elementsToSend(Ncpu);
    std::vector<int> numberElementsToReceive(Ncpu);
    for (int i = 0; i < Ncpu; i++) { numberElementsToReceive[i] = 0; }
    std::vector< std::vector<int> > elementsToReceive(Ncpu);

    for (int i = m_numberBoundFaces; i < m_numberElements; i++)
    {
      int numberOtherCPUs = m_elements[i]->getNumberAutresCPU();
      if (numberOtherCPUs != 0) // Elements that communicate
      {
        if (m_elements[i]->getCPU() == rankCpu) // The element belongs to the CPU
        {
          for (int p = 0; p < numberOtherCPUs; p++)
          {
            int numCPU = m_elements[i]->getAutreCPU(p);
            numberElementsToSend[numCPU]++;
            elementsToSend[numCPU].push_back(i);
          }
        }
        else // The element belongs to another CPU
        {
          int numCPU = m_elements[i]->getCPU();
          numberElementsToReceive[numCPU]++;
          elementsToReceive[numCPU].push_back(i);
        }
      }
    }
    for (int v = 0; v < Ncpu; v++)
    {
      if (numberElementsToSend[v] != 0) parallel.setNeighbour(v);
      //parallel.addElementToSend(v, buffer, numberElementsAEnvoyer[v]);
      for(int i=0;i<numberElementsToSend[v];++i)
      {
          const auto buffer = elementsToSend[v][i] - m_numberBoundFaces;
          parallel.addElementToSend(v,cells[buffer]);
      }
      for (int i = 0; i < numberElementsToReceive[v]; i++) { 
          const auto buffer= elementsToReceive[v][i] - m_numberBoundFaces;
          parallel.addElementToReceive(v, cells[buffer]);
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      std::cout << "    OK in " << t1 << " seconds" << std::endl;
    }

    if (rankCpu == 0)
    {
      std::cout << "... BUILDING GEOMETRY COMPLETE ";
      totalTime = clock() - totalTime; t1 = static_cast<float>(totalTime) / CLOCKS_PER_SEC;
      std::cout << " Total time of geometrical building : " << t1 << " seconds" << std::endl;
      std::cout << "------------------------------------------------------" << std::endl;
    }

    // Update of cellsGhost
    cellsGhost.insert(cellsGhost.begin(), cells.begin()+m_numberCellsCalcul, cells.end());
    cells.erase(cells.begin()+m_numberCellsCalcul, cells.end());
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MUSGmshV2::readMeshParallel()
{
  // Open mesh file from now on 
  try{
    // 1) Opening mesh file at Gmsh format V2
    // --------------------------------------
    if (rankCpu == 0)
    {
      std::cout << "------------------------------------------------------" << std::endl;
      std::cout << " C) READING MESH FILE " + m_nameMesh + "_CPUX.msh" + " IN PROGRESS ..." << std::endl;
    }
    std::stringstream flux;
    flux << rankCpu;
    m_meshFile = "./libMeshes/" + m_nameMesh + "_CPU" + flux.str() + ".msh";
    std::ifstream meshFile(m_meshFile.c_str(), std::ios::in);
    if (!meshFile) { throw ErrorECOGEN("file mesh absent :" + m_meshFile, __FILE__, __LINE__); }
    std::string currentLine;

    // Skiping Gmsh version + go to nodes
    // ----------------------------------
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);

    // 2) Filling m_nodes array
    // ------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0) { std::cout << "  1/Reading mesh nodes ..."; }
    meshFile >> m_numberNodes;
    meshFile.ignore(1, '\n');
    m_nodes = new Coord[m_numberNodes];
    int useless(0); double x, y, z;
    for (int i = 0; i < m_numberNodes; i++)
    {
      meshFile >> useless >> x >> y >> z;
      m_nodes[i].setXYZ(x, y, z);
    }
    meshFile.ignore(1, '\n');
    getline(meshFile, currentLine); // Skip keyword $EndNodes
    getline(meshFile, currentLine); // Skip keyword $Elements
    if (rankCpu == 0) { std::cout << "OK" << std::endl; }

    // 3) 1D/2D/3D elements are stored in m_elements array / counting
    // --------------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0) { std::cout << "  2/Reading internal 1D/2D/3D elements ..."; }
    meshFile >> m_numberElements;
    meshFile.ignore(1, '\n');
    // Allocation array of elements
    m_elements = new ElementNS*[m_numberElements];
    // Reading elements and assigning geometric properties 
    m_numberElements1D = 0, m_numberElements2D = 0, m_numberElements3D = 0;
    int posBeginElement = static_cast<int>(meshFile.tellg()); // id of the position of elements in file for quick return
    for (int i = 0; i < m_numberElements; i++)
    {
      this->readElement(m_nodes, meshFile, &m_elements[i]);
      if (m_elements[i]->getCPU() == rankCpu)
      {
        if (m_elements[i]->getTypeGmsh() == 1) { m_numberElements1D++; }
        else if (m_elements[i]->getTypeGmsh() <= 3) { m_numberElements2D++; m_totalSurface += m_elements[i]->getVolume(); }
        else if (m_elements[i]->getTypeGmsh() <= 7) { m_numberElements3D++; m_totalVolume += m_elements[i]->getVolume(); }
        else { throw ErrorECOGEN("Element type of .msh not handled by ECOGEN", __FILE__, __LINE__); }
      }
      else { m_numberGhostElements++; }
    }
    meshFile.ignore(1, '\n');
    getline(meshFile, currentLine);
    // Reading informations outside Gmsh
    getline(meshFile, currentLine);
    meshFile >> m_numberFacesParallel;
    meshFile.ignore(1, '\n');
    getline(meshFile, currentLine);
    meshFile >> m_numberInnerNodes;
    meshFile.close();
    // Calculation of the number of elements belonging to the CPU
    m_numberInnerElements = m_numberElements - m_numberGhostElements;
    if (rankCpu == 0) { std::cout << "OK" << std::endl; }
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MUSGmshV2::preProcessMeshFileForParallel()
{
  ElementNS **elementsGlobal;
  Coord *nodesGlobal;
  int frequenceImpressure(0);
  int numberNodesGlobal(0), numberElementsGlobal(0);
  int numberElements0D(0), numberElements1D(0), numberElements2D(0), numberElements3D(0);

  try {
    // Opening mesh file
    // -----------------
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << " C0) MESH FILE PRETRAITEMENT " + m_meshFile + " IN PROGRESS ..." << std::endl;
    clock_t totalTime(clock());
    std::ifstream meshFile(m_meshFile.c_str(), std::ios::in);
    if (!meshFile) { throw ErrorECOGEN("mesh file not found :" + m_meshFile, __FILE__, __LINE__); }
    std::string currentLine;

    // Skiping Gmsh version + go to nodes
    // ----------------------------------
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);

    // 1) Storing the vertex grid in a temporary array of nodes
    // --------------------------------------------------------
    std::cout << "  1/Mesh nodes reading ...";
    meshFile >> numberNodesGlobal;
    meshFile.ignore(1, '\n');
    nodesGlobal = new Coord[numberNodesGlobal];
    int useless(0); double x, y, z;
    for (int i = 0; i < numberNodesGlobal; i++)
    {
      meshFile >> useless >> x >> y >> z;
      nodesGlobal[i].setXYZ(x, y, z);
    }
    meshFile.ignore(1, '\n');
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);
    std::cout << "OK" << std::endl;

    // 2) 1D/2D/3D elements are stored in temporary elements array / counting
    // ----------------------------------------------------------------------
    std::cout << "  2/1D/2D/3D elements reading ...";
    meshFile >> numberElementsGlobal;
    meshFile.ignore(1, '\n');
    // Allocation array of elements
    elementsGlobal = new ElementNS*[numberElementsGlobal];
    // Reading elements and assigning geometric properties 
    int posBeginElement = static_cast<int>(meshFile.tellg()); // id of the position of elements in file for quick return
    for (int i = 0; i < numberElementsGlobal; i++)
    {
      this->readElement(nodesGlobal, meshFile, &elementsGlobal[i]);
      // Counting elements
      if (elementsGlobal[i]->getTypeGmsh() == 15) { numberElements0D++; }
      else if (elementsGlobal[i]->getTypeGmsh() == 1) { numberElements1D++; }
      else if (elementsGlobal[i]->getTypeGmsh() <= 3) { numberElements2D++; }
      else if (elementsGlobal[i]->getTypeGmsh() <= 7) { numberElements3D++; }
      else { Errors::errorMessage("Element type of .msh not handled by ECOGEN"); }
    }
    meshFile.close();

    int numberElementsBoundary;
    // Counting mesh and boundaries
    if (numberElements3D == 0 && numberElements2D == 0) // 1D case
    {
      numberElementsBoundary = numberElements0D;
    }
    else if (numberElements3D == 0) // 2D case
    {
      numberElementsBoundary = numberElements1D;
    }
    else // 3D case
    {
      numberElementsBoundary = numberElements2D;
    }

    std::cout << "OK" << std::endl;
    std::cout << "  -----------------------------------" << std::endl;
    std::cout << "    MESH GENERAL INFORMATIONS :" << std::endl;
    std::cout << "  -----------------------------------" << std::endl;
    std::cout << "    mesh nodes number : " << numberNodesGlobal << std::endl;
    std::cout << "    elements number : " << numberElementsGlobal << std::endl;

    // 3) Tracking and allocation by CPU
    // ---------------------------------
    std::cout << "  3/Attributing nodes to CPUs ..." << std::endl;
    clock_t tTemp(clock()); float t1(0);
    frequenceImpressure = std::max(numberElementsGlobal / 10, 1);
    std::vector< std::vector<int> > nodesCPU(Ncpu);
    std::vector< std::vector<int> > elementsCPU(Ncpu);
    int numCPU; int nodeCurrent; bool nodeExist;
    int numCPUMax(0); // To check if mesh suited or not
    // Seeking nodes of the inner mesh (without ghost) for each CPU
    for (int i = 0; i < numberElementsGlobal; i++)
    {
      if (i%frequenceImpressure == 0) { std::cout << "    " << 100 * i / numberElementsGlobal << "% ... " << std::endl; }
      numCPU = elementsGlobal[i]->getCPU();
      if (numCPU > numCPUMax) numCPUMax = numCPU;
      elementsCPU[numCPU].push_back(i); // Filling number of CPU-specific elements
      for (int n = 0; n < elementsGlobal[i]->getNumberNoeuds(); n++)
      {
        nodeCurrent = elementsGlobal[i]->getNumNoeud(n);
        nodeExist = false;
        for (unsigned int j = 0; j < nodesCPU[numCPU].size(); j++)
        {
          if (nodesCPU[numCPU][j] == nodeCurrent) { nodeExist = true; break; }
        }
        if (!nodeExist) { nodesCPU[numCPU].push_back(nodeCurrent); }
      }
    }
    int *numberInnerNodes = new int[Ncpu];
    for (int p = 0; p < Ncpu; p++)
    {
      numberInnerNodes[p] = nodesCPU[p].size();
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    std::cout << "    OK in " << t1 << " seconds" << std::endl;
    // Checking if mesh matches with CPUs used
    if (numCPUMax != Ncpu - 1) throw ErrorECOGEN("mesh file .msh not suited with the number of CPUs used - Generate the mesh and restard the test case", __FILE__, __LINE__);
    
    // 4) Creating array of faces
    // --------------------------
    // Test estimation of the number of faces (rough)
    std::cout << "  4/Creating faces array ..." << std::endl;
    tTemp = clock();
    frequenceImpressure = std::max((numberElementsGlobal - numberElementsBoundary) / 10, 1);
    int *numberFacesBuff = new int[Ncpu];
    int *iMaxFaces = new int[Ncpu];
    for (int p = 0; p < Ncpu; p++)
    {
      iMaxFaces[p] = 0;
      numberFacesBuff[p] = 0;
      for (unsigned int i = 0; i < elementsCPU[numCPU].size(); i++)
      {
        numberFacesBuff[p] += elementsGlobal[elementsCPU[numCPU][i]]->getNumberFaces();
      }
    }
    int ***facesBuff2 = new int**[Ncpu];
    int **sumNodesBuff2 = new int*[Ncpu];
    for (int p = 0; p < Ncpu; p++)
    {
      sumNodesBuff2[p] = new int[numberFacesBuff[p]];
      facesBuff2[p] = new int*[numberFacesBuff[p]];
      for (int i = 0; i < numberFacesBuff[p]; i++)
      {
        facesBuff2[p][i] = new int[4];
      }
    }
    for (int i = numberElementsBoundary; i < numberElementsGlobal; i++)
    {
      if ((i - numberElementsBoundary) % frequenceImpressure == 0) { std::cout << "    " << (100 * (i - numberElementsBoundary) / (numberElementsGlobal - numberElementsBoundary)) << "% ... " << std::endl; }
      int numCPU(elementsGlobal[i]->getCPU());
      elementsGlobal[i]->construitFacesSimplifie(iMaxFaces[numCPU], facesBuff2[numCPU], sumNodesBuff2[numCPU]);
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    std::cout << "    OK in " << t1 << " seconds" << std::endl;

    // 5) Seeking ghost elements which have a shared face with an inner element (communicators)
    // ----------------------------------------------------------------------------------------
    std::cout << "  5/Looking for ghosts elements ..." << std::endl;
    tTemp = clock();
    frequenceImpressure = std::max(numberElementsGlobal / 10, 1);
    int *numberFacesCommunicatingCPU = new int[Ncpu];
    for (int i = 0; i < Ncpu; i++) { numberFacesCommunicatingCPU[i] = 0; }
    for (int i = 0; i < numberElementsGlobal; i++)
    {
      if (i%frequenceImpressure == 0) { std::cout << "    " << (100 * i / numberElementsGlobal) << "% ... " << std::endl; }
      if (elementsGlobal[i]->getNumberAutresCPU() != 0)
      {
        std::vector<int> CPUAEnlever;
        for (int p = 0; p < elementsGlobal[i]->getNumberAutresCPU(); p++)
        {
          int numCPU(elementsGlobal[i]->getAutreCPU(p));
          // Checking if the element is communicating from a face
          // ****************************************************
          int numberNoeuds(elementsGlobal[i]->getNumberNoeuds());
          bool *isNoeudInterne = new bool[numberNoeuds];
          for (int n = 0; n < numberNoeuds; n++)
          {
            isNoeudInterne[n] = false;
            nodeCurrent = elementsGlobal[i]->getNumNoeud(n);
            for (int j = 0; j < numberInnerNodes[numCPU]; j++)
            {
              if (nodesCPU[numCPU][j] == nodeCurrent)
              {
                isNoeudInterne[n] = true;
                break;
              }
            }
          }
          // Determination of the number of communicating faces 
          // **************************************************
          //int numberFacesCommunicating(elementsGlobal[i]->compteFaceCommunicante(facesTemp[numCPU],sommeNoeudsTemp[numCPU]));
          int numberFacesCommunicating(elementsGlobal[i]->compteFaceCommunicante(iMaxFaces[numCPU], facesBuff2[numCPU], sumNodesBuff2[numCPU]));
          if (numberFacesCommunicating > 0)
          {
            // The element is communicating, we add it as well as his nodes
            elementsCPU[numCPU].push_back(i);
            for (int n = 0; n < numberNoeuds; n++)
            {
              if (!isNoeudInterne[n]) { nodesCPU[numCPU].push_back(elementsGlobal[i]->getNumNoeud(n)); }
            }
            numberFacesCommunicatingCPU[numCPU] += numberFacesCommunicating;
          }
          else
          {
            CPUAEnlever.push_back(p);
          }
          delete[] isNoeudInterne;

        } // End CPU
        elementsGlobal[i]->enleveCPUAutres(CPUAEnlever);
      }
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    std::cout << "    OK in " << t1 << " seconds" << std::endl;

    for (int p = 0; p < Ncpu; p++)
    {
      for (int i = 0; i < numberFacesBuff[p]; i++) { delete facesBuff2[p][i]; }
      delete sumNodesBuff2[p]; delete facesBuff2[p];
    }
    delete[]facesBuff2; delete[]sumNodesBuff2;
    delete[]numberFacesBuff; delete[]iMaxFaces;

    // 6) Writing mesh file for each CPU
    // ------------------------------------------------
    std::cout << "  6/Printing mesh files for each of " << Ncpu << " CPU ..." << std::endl;
    tTemp = clock();
    for (int p = 0; p < Ncpu; p++)
    {
      std::stringstream flux;
      flux << p;
      std::string meshFileCPU("./libMeshes/" + m_nameMesh + "_CPU" + flux.str() + ".msh");
      std::cout << "    print '" << meshFileCPU.c_str() << "' in progress ...' " << std::endl;
      std::ofstream fileStream;
      fileStream.open(meshFileCPU.c_str());

      fileStream << "$MeshFormat" << std::endl;
      fileStream << "2.2 0 8" << std::endl;
      fileStream << "$EndMeshFormat" << std::endl;
      fileStream << "$Nodes" << std::endl;
      //Reordonne les noeuds
      //sort(nodesCPU[p].begin(), nodesCPU[p].end());
      fileStream << nodesCPU[p].size() << std::endl;
      for (unsigned int i = 0; i < nodesCPU[p].size(); i++)
      {
        fileStream << i + 1 << " " << nodesGlobal[nodesCPU[p][i]].getX() << " " << nodesGlobal[nodesCPU[p][i]].getY() << " " << nodesGlobal[nodesCPU[p][i]].getZ() << std::endl;
      }
      fileStream << "$EndNodes" << std::endl;
      fileStream << "$Elements" << std::endl;
      fileStream << elementsCPU[p].size() << std::endl;
      for (unsigned int i = 0; i < elementsCPU[p].size(); i++)
      {
        int e(elementsCPU[p][i]);
        fileStream << i + 1 << " " << elementsGlobal[e]->getTypeGmsh();
        fileStream << " " << 2 + 1 + 1 + elementsGlobal[e]->getNumberAutresCPU();
        fileStream << " " << elementsGlobal[e]->getAppartenancePhysique();
        fileStream << " " << elementsGlobal[e]->getAppartenanceGeometrique();
        fileStream << " " << elementsGlobal[e]->getNumberAutresCPU() + 1;
        fileStream << " " << elementsGlobal[e]->getCPU() + 1;
        for (int cpuAutre = 0; cpuAutre < elementsGlobal[e]->getNumberAutresCPU(); cpuAutre++)
        {
          fileStream << " " << -(elementsGlobal[e]->getAutreCPU(cpuAutre) + 1);
        }
        for (int n = 0; n < elementsGlobal[e]->getNumberNoeuds(); n++)
        {
          // Renumbering locale
          nodeCurrent = elementsGlobal[e]->getNumNoeud(n);
          for (int j = 0; j < static_cast<int>(nodesCPU[p].size()); j++)
          {
            if (nodesCPU[p][j] == nodeCurrent)
            {
              fileStream << " " << j + 1; break;
            }
          }
        }
        fileStream << std::endl;
      }
      fileStream << "$EndElements" << std::endl;
      // Usefull aditionnal data
      fileStream << "Information not read by Gmsh : number of communicating faces" << std::endl;
      fileStream << numberFacesCommunicatingCPU[p] << std::endl;
      fileStream << "Information not read by Gmsh : number of inner nodes (except ghosts)" << std::endl;
      fileStream << numberInnerNodes[p] << std::endl;
      fileStream.close();
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    std::cout << "    OK in " << t1 << " seconds" << std::endl;

    std::cout << "... MESH FILE PRE-PROCESSING COMPLETE ";
    totalTime = clock() - totalTime; t1 = static_cast<float>(totalTime) / CLOCKS_PER_SEC;
    std::cout << " Total time of pre-processing : " << t1 << " seconds" << std::endl;

    // Deallocations
    for (int i = 0; i < numberElementsGlobal; i++)
    {
      delete elementsGlobal[i];
    }
    delete[] elementsGlobal;
    delete[] nodesGlobal;
    delete[] numberInnerNodes;
    delete[] numberFacesCommunicatingCPU;

  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MUSGmshV2::readElement(const Coord *nodesTable, std::ifstream &meshFile, ElementNS **element)
{
  int numberElement,numberTags,typeElement,numberPhysicEntity,numberGeometricEntity;
  meshFile >> numberElement >> typeElement;
  // Reading tags
  meshFile >> numberTags;
  meshFile >> numberPhysicEntity;    // Physic belonging
  meshFile >> numberGeometricEntity; // Geometric belonging

  // 1) Assignation of the number of vertices according to element type
  // ------------------------------------------------------------------
  switch (typeElement)
  {
    case 1: // Segment (two points)
      *element = new ElementSegment;
      m_numberSegments++;
      break;
    case 2: // Triangle (three points)
      *element = new ElementTriangle;
      m_numberTriangles++;
      break;
    case 3: // Quadrangle (four points)
      *element = new ElementQuadrangle;
      m_numberQuadrangles++;
      break;
    case 4: // Tetrahedron (four points)
      *element = new ElementTetrahedron;
      m_numberTetrahedrons++;
      break;
    //case 7: // Quadrangular pyramid (five points) // This element seems to not work with Gmsh, volumes of elements seem to cause this issue 
    //  *element = new ElementPyramid; 
    //  m_numberPyramids++;
    //  break;
    case 15: // Point (a vertex)
      *element = new ElementPoint;
      m_numberPoints++;
      break;
    case 5: // Hexahedron (eight points)
      *element = new ElementHexahedron;
      m_numberHexahedrons++;
      break;
    case 6: // Prism (six points)
      *element = new ElementPrism;
      m_numberHexahedrons++;
      break;
    default:
      Errors::errorMessage("Element type of file .msh inknown for ECOGEN");
      break;
  }
 
  // 2) Specific tags for parallel mesh
  // -----------------------------------
  int numberCPU(0);
  if (numberTags > 2)
  { 
    meshFile >> numberCPU; // Number of mesh partitions to which the elment belongs
    int *numCPU = new int[numberCPU];
    for (int tag = 0; tag < numberCPU; tag++){ meshFile >> numCPU[tag]; }
    (*element)->setAppartenanceCPU(numCPU, numberCPU);
    delete[] numCPU;
  }

  // 3) Building the element and its properties
  // ------------------------------------------
  int currentNode;
  int *numNode = new int[(*element)->getNumberNoeuds()];
  Coord *node = new Coord[(*element)->getNumberNoeuds()];
  for (int i = 0; i < (*element)->getNumberNoeuds(); i++)
  {
    meshFile >> currentNode;
    numNode[i] = currentNode - 1; // Offset because array start at 0
    node[i] = nodesTable[currentNode - 1];
  }
  int indexElement(numberElement - 1);

  if (m_switchTags != false) { // Option to switch between physical/geometrical tag when Gmsh converts .mesh to .msh
	  int buff;
	  buff = numberPhysicEntity;
	  numberPhysicEntity = numberGeometricEntity;
	  numberGeometricEntity = buff;
  }

  (*element)->construitElement(numNode, node, numberPhysicEntity, numberGeometricEntity, indexElement);

  delete[] node;
  delete[] numNode;
}

