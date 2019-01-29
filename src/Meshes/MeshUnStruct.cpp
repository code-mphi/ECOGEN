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

//! \file      MeshUnStruct.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "MeshUnStruct.h"
#include "../Errors.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

MeshUnStruct::MeshUnStruct(const string &fichierMesh) :
  Mesh(),
  m_fichierMesh(fichierMesh),
  m_nameMesh(fichierMesh),
  m_numberNoeuds(0),
  m_numberNoeudsInternes(0),
  m_numberElementsInternes(0),
  m_numberElementsFantomes(0),
  m_numberElementsCommunicants(0),
  m_numberCellsFantomes(0),
  m_numberFacesParallele(0),
  m_numberElements0D(0),
  m_numberElements1D(0),
  m_numberElements2D(0),
  m_numberElements3D(0),
  m_numberSegments(0),
  m_numberTriangles(0),
  m_numberQuadrangles(0),
  m_numberTetrahedrons(0),
  m_numberPyramids(0),
  m_numberPoints(0),
  m_numberHexahedrons(0),
  m_totalSurface(0.),
  m_totalVolume(0.)
{
  m_nameMesh = m_fichierMesh;
  m_nameMesh.resize(m_nameMesh.size() - 4); //On enleve l extension
  m_type = UNS;
}

//***********************************************************************

MeshUnStruct::~MeshUnStruct(){
  for (int f = 0; f < m_numberFacesTotal; f++) { delete m_faces[f]; }
  delete[] m_faces;
  for (int e = 0; e < m_numberElements; e++) { delete m_elements[e]; }
  delete[] m_elements;
  delete[] m_noeuds;
  for (unsigned int l = 0; l < m_lim.size(); l++) { delete m_lim[l]; }
}

//***********************************************************************

void MeshUnStruct::attributLimites(std::vector<BoundCond*> &boundCond)
{
  //Recherche number physique de limite le plus grand
  int maxNumLim(0);
  for (unsigned int i = 0; i < boundCond.size(); i++) {
    if (boundCond[i]->getNumPhys() > maxNumLim) maxNumLim = boundCond[i]->getNumPhys();
  }
  //Attribution des limites dans le tableau m_lim dans l'ordre
  int limite(1);
  int limiteTrouvee(0);
  while (limite <= maxNumLim) {
    for (unsigned int i = 0; i < boundCond.size(); i++) {
      if (boundCond[i]->getNumPhys() == limite) {
        m_lim.push_back(boundCond[i]);
        limiteTrouvee = 1; break;
      }
    }
    if (!limiteTrouvee) m_lim.push_back(new BoundCondAbs(limite));
    limite++; limiteTrouvee = 0;
  }
}

//***********************************************************************

int MeshUnStruct::initializeGeometrie(Cell ***cells, CellInterface ***bord, bool pretraitementParallele, string ordreCalcul)
{
  try {
    if (Ncpu == 1) { this->initializeGeometrieMonoCPU(cells, bord, ordreCalcul); }
    else {
      //Pretraitement du file de mesh par le CPU 0
      if (pretraitementParallele) {
        if (rankCpu == 0) { this->pretraitementFichierMeshGmsh(); }
        MPI_Barrier(MPI_COMM_WORLD);
      }
      this->initializeGeometrieParallele(cells, bord, ordreCalcul);
    }
    return m_geometrie;
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MeshUnStruct::initializeGeometrieMonoCPU(Cell ***cells, CellInterface ***bord, string ordreCalcul)
{
  try {
    //1) Lecture noeuds et elements
    //-----------------------------
    vector<ElementNS*>* voisinsNoeuds; //Dimensions : nb_noeuds
    this->lectureGeometrieGmsh(&voisinsNoeuds);  //Remplissage m_noeuds et m_elements

    //CAUTION: Ordering of m_elements is important. Faces first, then cells.

    cout << "------------------------------------------------------" << endl;
    cout << " B) BUILDING GEOMETRY ..." << endl;

    //2) Attribution des cells a leur element geometrique
    //------------------------------------------------------
    //Comptage mailles et limites
    if (m_numberElements3D == 0 && m_numberElements2D == 0) //Cas 1D
    {
      m_numberCellsCalcul = m_numberElements1D;
      m_numberFacesLimites = m_numberElements0D;
      m_geometrie = 1;
    }
    else if (m_numberElements3D == 0) //Cas 2D
    {
      m_numberCellsCalcul = m_numberElements2D;
      m_numberFacesLimites = m_numberElements1D;
      m_geometrie = 2;
    }
    else //Cas 3D
    {
      m_numberCellsCalcul = m_numberElements3D;
      m_numberFacesLimites = m_numberElements2D;
      m_geometrie = 3;
    }
    m_numberCellsTotal = m_numberCellsCalcul;

    //Dimensionnement du tableau de cells
    (*cells) = new Cell*[m_numberCellsCalcul];

    //Attribution Absorption aux limites manquantes
    unsigned int nbLimites(0);
    for (int i = 0; i < m_numberFacesLimites; i++) {
      unsigned int appPhys(m_elements[i]->getAppartenancePhysique());
      if (appPhys > nbLimites) { nbLimites = appPhys; }
    }
    for (unsigned int i = m_lim.size(); i < nbLimites; i++) {
      m_lim.push_back(new BoundCondAbs);
    }

    //Attribution des cells aux elements et comptage faces internes
    m_numberFacesInternes = 0;
    for (int i = 0; i < m_numberCellsCalcul; i++)
    {
      if (ordreCalcul == "FIRSTORDER") { (*cells)[i] = new Cell; }
      else { (*cells)[i] = new CellO2; }
      (*cells)[i]->setElement(m_elements[i + m_numberFacesLimites], i);
      m_numberFacesInternes += m_elements[i + m_numberFacesLimites]->getNumberFaces();
    }

    //-----------Ajout construction voisins des mailles pour ordre 2----------
    vector<ElementNS *>* voisins; //vecteur temporaire des voisins
    voisins = new vector<ElementNS *>[m_numberCellsCalcul];
    for (int i = 0; i < m_numberCellsCalcul; i++)
    {
      ElementNS *e(m_elements[i + m_numberFacesLimites]); //choix element
      //cout << "el " << i + m_numberFacesLimites << " : " ;
      //vector<ElementNS *> voisins; //vecteur temporaire des voisins
      //1) Construction du vecteur de voisins
      for (int n = 0; n < e->getNumberNoeuds(); n++) { //boucle noeud de element e
        int numNoeud = e->getNumNoeud(n);
        for (unsigned int v = 0; v < voisinsNoeuds[numNoeud].size(); v++) { //boucle voisin du noeud n
          bool ajoute(true);
          if (voisinsNoeuds[numNoeud][v]->getIndex() == e->getIndex()) ajoute = false;
          //else if(voisinsNoeuds[numNoeud][v]->getIndex() < m_numberFacesLimites) ajoute = false;
          else {
            for (unsigned int vo = 0; vo < voisins[i].size(); vo++) {
              if (voisinsNoeuds[numNoeud][v]->getIndex() == voisins[i][vo]->getIndex()) {
                ajoute = false; break;
              }
            }
          }
          if (ajoute) {
            voisins[i].push_back(voisinsNoeuds[numNoeud][v]);
            //cout << voisinsNoeuds[numNoeud][v]->getIndex() << " ";
          }
        }
      }
      //cout << endl;
    }
    //------------------------------------------------------------------------

    m_numberFacesInternes -= m_numberFacesLimites; //On enleve les limites
    m_numberFacesInternes /= 2; //Les faces internes sont toute comptees 2 fois => on retabli la verite !
    m_numberFacesTotal = m_numberFacesInternes + m_numberFacesLimites;

    double volTot(0.);
    for (int i = 0; i < m_numberCellsCalcul; i++)
    {
      volTot += (*cells)[i]->getElement()->getVolume();
    }
    cout << "Total volume : " << volTot << endl;

    //3) Construction de la table de connectivite
    //-------------------------------------------
    //Dimensionnement du tableau de faces
    (*bord) = new CellInterface*[m_numberFacesTotal];
    m_faces = new FaceNS*[m_numberFacesTotal];
    int **facesTemp; int *sommeNoeudsTemp; //On cre un tableau temporaire de faces pour accelerer la recherche d'existance
    facesTemp = new int*[m_numberFacesTotal + 1];
    sommeNoeudsTemp = new int[m_numberFacesTotal + 1];
    //Determination du number de noeuds max pour les faces
    int tailleFace; //Sera initialize a la taille maximale
    if (m_numberElements3D != 0)
    {
      if (m_numberQuadrangles != 0) { tailleFace = 4; }
      else if (m_numberTriangles != 0) { tailleFace = 3; }
      else { Errors::errorMessage("Probleme dans initializeGeometrieMonoCPU pour initialization du tableau facesTemp"); }
    }
    else if (m_numberElements2D != 0) { tailleFace = 2; }
    else { tailleFace = 1; }
    for (int i = 0; i < m_numberFacesTotal + 1; i++)
    {
      facesTemp[i] = new int[tailleFace];
    }

    //Faces internes
    int indexMaxFaces(0);
    clock_t tTemp(clock());
    float t1(0.);
    cout << "  1/Building faces ..." << endl;
    int frequenceImpressure(max((m_numberElements - m_numberFacesLimites) / 10, 1));
    for (int i = m_numberFacesLimites; i < m_numberElements; i++)
    {
      if ((i - m_numberFacesLimites) % frequenceImpressure == 0) { cout << "    " << (100 * (i - m_numberFacesLimites) / (m_numberElements - m_numberFacesLimites)) << "% ... " << endl; }
      m_elements[i]->construitFaces(m_noeuds, m_faces, indexMaxFaces, facesTemp, sommeNoeudsTemp);
    }
    for (int i = 0; i < m_numberFacesTotal + 1; i++) { delete facesTemp[i]; }
    delete[] facesTemp; delete[] sommeNoeudsTemp;
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK in " << t1 << " seconds" << endl;

    //Limites
    cout << "  2/Boundary elements attribution to boundary faces ..." << endl;
    tTemp = clock();
    frequenceImpressure = max(m_numberFacesLimites / 10, 1);
    for (int i = 0; i < m_numberFacesLimites; i++)
    {
      if (i%frequenceImpressure == 0)
      {
        cout << "    " << (100 * i / m_numberFacesLimites) << "% ... " << endl;// << " -> "; 
      }
      //Attribution de la limite
      m_elements[i]->attributFaceLimite(m_noeuds, m_faces, indexMaxFaces);
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK in " << t1 << " seconds" << endl;

    //Liaison Geometrie/Bords de compute
    cout << "  3/Linking Geometries -> Physics ..." << endl;
    tTemp = clock();
    int iMailleG, iMailleD;
    for (int i = 0; i < m_numberFacesTotal; i++)
    {
      if (m_faces[i]->getEstLimite())
      {
        int appPhys(m_faces[i]->getElementDroite()->getAppartenancePhysique() - 1); //appartenance - 1 pour tableau commencant a zero
        if (appPhys >= static_cast<int>(m_lim.size()) || appPhys < 0) { Errors::errorMessage("Number de conditions aux limites non adapte"); }
        m_lim[appPhys]->creeLimite(&(*bord)[i]);
        (*bord)[i]->setFace(m_faces[i]);
        iMailleG = m_faces[i]->getElementGauche()->getIndex() - m_numberFacesLimites;
        iMailleD = iMailleG;
      }
      else
      {
        if (ordreCalcul == "FIRSTORDER") { (*bord)[i] = new CellInterface; }
        else { (*bord)[i] = new CellInterfaceO2; }
        (*bord)[i]->setFace(m_faces[i]);
        iMailleG = m_faces[i]->getElementGauche()->getIndex() - m_numberFacesLimites;
        iMailleD = m_faces[i]->getElementDroite()->getIndex() - m_numberFacesLimites;
      }
      (*bord)[i]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
      (*cells)[iMailleG]->addBoundary((*bord)[i]);
      (*cells)[iMailleD]->addBoundary((*bord)[i]);

    } //Fin face
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK in " << t1 << " seconds" << endl;
    cout << "... BUILDING GEOMETRY COMPLETE " << endl;
    cout << "------------------------------------------------------" << endl;

    delete[] voisins;
    delete[] voisinsNoeuds;
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MeshUnStruct::initializeGeometrieParallele(Cell ***cells, CellInterface ***bord, string ordreCalcul)
{

  clock_t totalTime(clock());

  //1) Lecture noeuds et elements
  //-----------------------------
  try {
    this->lectureGeometrieGmshParallele(); //Remplissage de m_noeuds et m_elements
    if (rankCpu == 0)
    {
      cout << "------------------------------------------------------" << endl;
      cout << " B) BUILDING GEOMETRY..." << endl;
    }

    //2) Attribution des cells a leur element geometrique
    //------------------------------------------------------
    //Comptage mailles et limites
    if (m_numberElements3D == 0 && m_numberElements2D == 0) //Cas 1D
    {
      m_numberCellsCalcul = m_numberElements1D;
      m_numberFacesLimites = m_numberElements0D;
      m_geometrie = 1;
    }
    else if (m_numberElements3D == 0) //Cas 2D
    {
      m_numberCellsCalcul = m_numberElements2D;
      m_numberFacesLimites = m_numberElements1D;
      m_geometrie = 2;
    }
    else //Cas 3D
    {
      m_numberCellsCalcul = m_numberElements3D;
      m_numberFacesLimites = m_numberElements2D;
      m_geometrie = 3;
    }

    //Dimensionnement du tableau de cells
    m_numberCellsTotal = m_numberCellsCalcul + m_numberElementsFantomes;
    (*cells) = new Cell*[m_numberCellsTotal];

    //Attribution Absorption aux limites manquantes
    unsigned int nbLimites(0);
    for (int i = 0; i < m_numberFacesLimites; i++) {
      unsigned int appPhys(m_elements[i]->getAppartenancePhysique());
      if (appPhys > nbLimites) { nbLimites = appPhys; }
    }
    for (unsigned int i = m_lim.size(); i < nbLimites; i++) {
      m_lim.push_back(new BoundCondAbs);
    }

    //Attribution des cells aux elements et comptage faces internes
    m_numberFacesInternes = 0;
    for (int i = 0; i < m_numberCellsTotal; i++)
    {
      if (ordreCalcul == "FIRSTORDER") { (*cells)[i] = new Cell; }
      else { (*cells)[i] = new CellO2; }
      (*cells)[i]->setElement(m_elements[i + m_numberFacesLimites], i);
      if (i < m_numberCellsCalcul) { m_numberFacesInternes += m_elements[i + m_numberFacesLimites]->getNumberFaces(); }
    }
    m_numberFacesInternes -= m_numberFacesLimites + m_numberFacesParallele; //On enleve les limites et les faces communicantes
    m_numberFacesInternes /= 2; //Les faces internes sont toutes comptees 2 fois => on retabli la verite !
    m_numberFacesTotal = m_numberFacesInternes + m_numberFacesLimites + m_numberFacesParallele;

    //3) Construction de la table de connectivite interne
    //---------------------------------------------------
    //Dimensionnement des tableaux de faces
    (*bord) = new CellInterface*[m_numberFacesTotal];
    m_faces = new FaceNS*[m_numberFacesTotal];

    //On cree un tableau temporaire de faces pour accelerer la recherche d'existance
    int **facesTemp; int *sommeNoeudsTemp;
    facesTemp = new int*[m_numberFacesTotal + 1];
    sommeNoeudsTemp = new int[m_numberFacesTotal + 1];
    //Determination du number de noeuds max pour les faces
    int tailleFace; //Sera initialize a la taille maximale
    if (m_numberElements3D != 0)
    {
      if (m_numberQuadrangles != 0) { tailleFace = 4; }
      else if (m_numberTriangles != 0) { tailleFace = 3; }
      else { Errors::errorMessage("Probleme dans initializeGeometrieMonoCPU pour initialization du tableau facesTemp"); }
    }
    else if (m_numberElements2D != 0) { tailleFace = 2; }
    else { tailleFace = 1; }
    for (int i = 0; i < m_numberFacesTotal + 1; i++) // Le +1 est utilise pour la recherche existance de face
    {
      facesTemp[i] = new int[tailleFace]; // Inconnue sur le number de points d une face (maxi 4 a priori)
    }

    //Faces internes
    //--------------
    int indexMaxFaces(0);
    MPI_Barrier(MPI_COMM_WORLD);
    int frequenceImpressure;
    clock_t tTemp(clock());
    float t1(0.);
    if (rankCpu == 0)
    {
      cout << "  1/Building faces ..." << endl;
      frequenceImpressure = max((m_numberElementsInternes - m_numberFacesLimites) / 10, 1);
    }
    for (int i = m_numberFacesLimites; i < m_numberElementsInternes; i++)
    {
      if (rankCpu == 0 && (i - m_numberFacesLimites) % frequenceImpressure == 0) { cout << "    " << (100 * (i - m_numberFacesLimites) / (m_numberElementsInternes - m_numberFacesLimites)) << "% ... " << endl; }
      //Construction
      m_elements[i]->construitFaces(m_noeuds, m_faces, indexMaxFaces, facesTemp, sommeNoeudsTemp);
    }
    for (int i = 0; i < m_numberFacesTotal + 1; i++) { delete facesTemp[i]; }
    delete[] facesTemp; delete[] sommeNoeudsTemp;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      cout << "    OK in " << t1 << " seconds" << endl;
    }

    //Limites
    //-------
    tTemp = clock();
    if (rankCpu == 0)
    {
      cout << "  2/Boundary elements attribution to boundary faces ..." << endl;
      frequenceImpressure = max(m_numberFacesLimites / 10, 1);
    }
    for (int i = 0; i < m_numberFacesLimites; i++)
    {
      if (rankCpu == 0 && i%frequenceImpressure == 0) { cout << "    " << (100 * i / m_numberFacesLimites) << "% ... " << endl; }
      //Attribution de la limite
      m_elements[i]->attributFaceLimite(m_noeuds, m_faces, indexMaxFaces);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      cout << "    OK in " << t1 << " seconds" << endl;
    }

    //Communications
    //--------------
    tTemp = clock();
    if (rankCpu == 0)
    {
      cout << "  3/Ghost cells attribution to communicating faces ..." << endl;
      frequenceImpressure = max((m_numberElements - m_numberElementsInternes) / 10, 1);
    }
    for (int i = m_numberElementsInternes; i < m_numberElements; i++)
    {
      if (rankCpu == 0 && (i - m_numberElementsInternes) % frequenceImpressure == 0) { cout << "    " << (100 * (i - m_numberElementsInternes) / (m_numberElements - m_numberElementsInternes)) << "% ... " << endl; }
      //Attribution de la limite communicante
      m_elements[i]->attributFaceCommunicante(m_noeuds, m_faces, indexMaxFaces, m_numberNoeudsInternes);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      cout << "    OK in " << t1 << " seconds" << endl;
    }

    //Liaison Geometrie/Bords de compute
    //---------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    tTemp = clock();
    if (rankCpu == 0)
    {
      cout << "  4/Linking Geometries -> Physics ..." << endl;
      frequenceImpressure = max(m_numberFacesTotal / 10, 1);
    }
    int iMailleG, iMailleD;

    for (int i = 0; i < m_numberFacesTotal; i++)
    {
      if (rankCpu == 0 && i%frequenceImpressure == 0) { cout << "    " << 100 * i / m_numberFacesTotal << "% ... " << endl; }
      //Faces limites (limite physique ou communication)
      if (m_faces[i]->getEstLimite())
      {
        //Communication
        if (m_faces[i]->getEstComm())
        {
          if (ordreCalcul == "FIRSTORDER") { (*bord)[i] = new CellInterface; }
          else { (*bord)[i] = new CellInterfaceO2; }
          (*bord)[i]->setFace(m_faces[i]);
          iMailleG = m_faces[i]->getElementGauche()->getNumCellAssociee();
          iMailleD = m_faces[i]->getElementDroite()->getNumCellAssociee();
        }
        //Limite physique
        else
        {
          int appPhys(m_faces[i]->getElementDroite()->getAppartenancePhysique() - 1); //appartenance - 1 pour tableau commencant a zero
          if (appPhys >= static_cast<int>(m_lim.size()) || appPhys < 0) { throw ErrorECOGEN("Number de conditions aux limites non adapte", __FILE__, __LINE__); }
          m_lim[appPhys]->creeLimite(&(*bord)[i]);
          (*bord)[i]->setFace(m_faces[i]);
          iMailleG = m_faces[i]->getElementGauche()->getNumCellAssociee();
          iMailleD = iMailleG;
        }
      }
      //face interne au domain
      else
      {
        if (ordreCalcul == "FIRSTORDER") { (*bord)[i] = new CellInterface; }
        else { (*bord)[i] = new CellInterfaceO2; }
        (*bord)[i]->setFace(m_faces[i]);
        iMailleG = m_faces[i]->getElementGauche()->getNumCellAssociee();
        iMailleD = m_faces[i]->getElementDroite()->getNumCellAssociee();
      }
      (*bord)[i]->initialize((*cells)[iMailleG], (*cells)[iMailleD]);
      (*cells)[iMailleG]->addBoundary((*bord)[i]);
      (*cells)[iMailleD]->addBoundary((*bord)[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      cout << "    OK in " << t1 << " seconds" << endl;
    }

    //4) Construction de la table de connectivite parallele CPUs
    //----------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    tTemp = clock();
    if (rankCpu == 0)
    {
      cout << "  5/Building connectivity tables for CPUs ..." << endl;
      frequenceImpressure = max((m_numberElements - m_numberFacesLimites) / 10, 1);
    }

    vector<int> numberElementsAEnvoyer(Ncpu);
    for (int i = 0; i < Ncpu; i++) { numberElementsAEnvoyer[i] = 0; }
    vector< vector<int> > elementsAEnvoyer(Ncpu);
    vector<int> numberElementsARecevoir(Ncpu);
    for (int i = 0; i < Ncpu; i++) { numberElementsARecevoir[i] = 0; }
    vector< vector<int> > elementsARecevoir(Ncpu);

    for (int i = m_numberFacesLimites; i < m_numberElements; i++)
    {
      int numberAutresCPUs = m_elements[i]->getNumberAutresCPU();
      if (numberAutresCPUs != 0) //Elements qui communique
      {
        if (m_elements[i]->getCPU() == rankCpu) //L element appartient au CPU
        {
          for (int p = 0; p < numberAutresCPUs; p++)
          {
            int numCPU = m_elements[i]->getAutreCPU(p);
            numberElementsAEnvoyer[numCPU]++;
            elementsAEnvoyer[numCPU].push_back(i);
          }
        }
        else // L'element appartient a un autre CPU
        {
          int numCPU = m_elements[i]->getCPU();
          numberElementsARecevoir[numCPU]++;
          elementsARecevoir[numCPU].push_back(i);
        }
      }
    }
    int *buffer;
    for (int v = 0; v < Ncpu; v++)
    {
      string whichCpuAmIForNeighbour("");
      if (numberElementsAEnvoyer[v] != 0) parallel.setNeighbour(v, whichCpuAmIForNeighbour);
      buffer = new int[numberElementsAEnvoyer[v]];
      for (int i = 0; i < numberElementsAEnvoyer[v]; i++) { buffer[i] = elementsAEnvoyer[v][i] - m_numberFacesLimites; }
      parallel.setElementsToSend(v, buffer, numberElementsAEnvoyer[v]);
      delete[] buffer;
      buffer = new int[numberElementsARecevoir[v]];
      for (int i = 0; i < numberElementsARecevoir[v]; i++) { buffer[i] = elementsARecevoir[v][i] - m_numberFacesLimites; }
      parallel.setElementsToReceive(v, buffer, numberElementsARecevoir[v]);
      delete[] buffer;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0)
    {
      tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
      cout << "    OK in " << t1 << " seconds" << endl;
    }

    if (rankCpu == 0)
    {
      cout << "... BUILDING GEOMETRY COMPLETE ";
      totalTime = clock() - totalTime; t1 = static_cast<float>(totalTime) / CLOCKS_PER_SEC;
      cout << " Total time of geometrical building : " << t1 << " seconds" << endl;
      cout << "------------------------------------------------------" << endl;
    }
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MeshUnStruct::pretraitementFichierMeshGmsh()
{
  ElementNS **elementsGlobal;
  Coord *noeudsGlobal;
  int frequenceImpressure(0);
  int numberNoeudsGlobal(0), numberElementsGlobal(0);
  int numberElements0D(0), numberElements1D(0), numberElements2D(0), numberElements3D(0);

  try {
    //Ouverture du file de mesh
    //-------------------------------------
    cout << "------------------------------------------------------" << endl;
    cout << " 0) MESH FILE PRETRAITEMENT " + m_fichierMesh + " IN PROGRESS..." << endl;
    clock_t totalTime(clock());
    m_fichierMesh = "./libMeshes/" + m_fichierMesh;
    ifstream fichierMesh(m_fichierMesh.c_str(), ios::in);
    if (!fichierMesh) { throw ErrorECOGEN("file mesh absent :" + m_fichierMesh, __FILE__, __LINE__); }
    string ligneCourante;
    getline(fichierMesh, ligneCourante);
    getline(fichierMesh, ligneCourante);
    getline(fichierMesh, ligneCourante);
    getline(fichierMesh, ligneCourante);

    //1) Stockage de la grille de vertex dans tableau temporaire de Noeuds
    //-------------------------------------------------------
    cout << "  1/Mesh nodes reading ...";
    fichierMesh >> numberNoeudsGlobal;
    fichierMesh.ignore(1, '\n');
    noeudsGlobal = new Coord[numberNoeudsGlobal];
    int inutile(0); double x, y, z;
    for (int i = 0; i < numberNoeudsGlobal; i++)
    {
      fichierMesh >> inutile >> x >> y >> z;
      noeudsGlobal[i].setXYZ(x, y, z);
    }
    fichierMesh.ignore(1, '\n');
    getline(fichierMesh, ligneCourante);
    getline(fichierMesh, ligneCourante);
    cout << "OK" << endl;

    //2) Recuperation des elements 1D/2D/3D dans le tableau temporaire d'elements et comptage
    //----------------------------------------------------------------------------
    cout << "  2/1D/2D/3D elements reading ...";
    fichierMesh >> numberElementsGlobal;
    fichierMesh.ignore(1, '\n');
    //Allocation tableau d elements
    elementsGlobal = new ElementNS*[numberElementsGlobal];
    //Lecture des elements et attributions proprietes geometriques
    int posDebutElement = static_cast<int>(fichierMesh.tellg()); //reperage debut des elements pour retour rapide
    for (int i = 0; i < numberElementsGlobal; i++)
    {
      this->lectureElementGmshV2(noeudsGlobal, fichierMesh, &elementsGlobal[i]);
      //Comptage des elements
      if (elementsGlobal[i]->getTypeGmsh() == 15) { numberElements0D++; }
      else if (elementsGlobal[i]->getTypeGmsh() == 1) { numberElements1D++; }
      else if (elementsGlobal[i]->getTypeGmsh() <= 3) { numberElements2D++; }
      else if (elementsGlobal[i]->getTypeGmsh() <= 7) { numberElements3D++; }
      else { Errors::errorMessage("Type element du .msh non gere dans ECOGEN"); }
    }
    fichierMesh.close();

    //Comptage mailles et limites
    int numberElementsFrontiere;
    //Comptage mailles et limites
    if (numberElements3D == 0 && numberElements2D == 0) //Cas 1D
    {
      numberElementsFrontiere = numberElements0D;
    }
    else if (numberElements3D == 0) //Cas 2D
    {
      numberElementsFrontiere = numberElements1D;
    }
    else //Cas 3D
    {
      numberElementsFrontiere = numberElements2D;
    }

    cout << "OK" << endl;
    cout << "  -----------------------------------" << endl;
    cout << "    MESH GENERAL INFORMATIONS :" << endl;
    cout << "  -----------------------------------" << endl;
    cout << "    mesh nodes number : " << numberNoeudsGlobal << endl;
    cout << "    elements number : " << numberElementsGlobal << endl;

    //3)Reperage et attribution par CPU
    //---------------------------------
    cout << "  3/Attributing nodes to CPUs ..." << endl;
    clock_t tTemp(clock()); float t1(0);
    frequenceImpressure = max(numberElementsGlobal / 10, 1);
    vector< vector<int> > noeudsCPU(Ncpu);
    vector< vector<int> > elementsCPU(Ncpu);
    int numCPU; int noeudCourant; bool noeudExiste;
    int numCPUMax(0); //pour verification mesh adapte ou non
    //Recherche des noeuds formant le mesh interne (hors fantomes) pour chaque CPU
    for (int i = 0; i < numberElementsGlobal; i++)
    {
      if (i%frequenceImpressure == 0) { cout << "    " << 100 * i / numberElementsGlobal << "% ... " << endl; }
      numCPU = elementsGlobal[i]->getCPU();
      if (numCPU > numCPUMax) numCPUMax = numCPU;
      elementsCPU[numCPU].push_back(i); //Remplissage number des elements propres au CPU
      for (int n = 0; n < elementsGlobal[i]->getNumberNoeuds(); n++)
      {
        noeudCourant = elementsGlobal[i]->getNumNoeud(n);
        noeudExiste = false;
        for (unsigned int j = 0; j < noeudsCPU[numCPU].size(); j++)
        {
          if (noeudsCPU[numCPU][j] == noeudCourant) { noeudExiste = true; break; }
        }
        if (!noeudExiste) { noeudsCPU[numCPU].push_back(noeudCourant); }
      }
    }
    int *numberNoeudsInterne = new int[Ncpu];
    for (int p = 0; p < Ncpu; p++)
    {
      numberNoeudsInterne[p] = noeudsCPU[p].size();
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK in " << t1 << " seconds" << endl;

    //Verification adaptation au mesh
    if (numCPUMax != Ncpu - 1) throw ErrorECOGEN("file mesh .msh non adapte au number de CPU - Generer le mesh et relancer le test", __FILE__, __LINE__);

    //4) Creation Tableau de faces
    //-----------------------------
    //Test estimation number de faces (approximative)
    cout << "  4/Creating faces array ..." << endl;
    tTemp = clock();
    frequenceImpressure = max((numberElementsGlobal - numberElementsFrontiere) / 10, 1);
    int *numberFacesTemp = new int[Ncpu];
    int *iMaxFaces = new int[Ncpu];
    for (int p = 0; p < Ncpu; p++)
    {
      iMaxFaces[p] = 0;
      numberFacesTemp[p] = 0;
      for (unsigned int i = 0; i < elementsCPU[numCPU].size(); i++)
      {
        numberFacesTemp[p] += elementsGlobal[elementsCPU[numCPU][i]]->getNumberFaces();
      }
    }
    int ***facesTemp2 = new int**[Ncpu];
    int **sommeNoeudsTemp2 = new int*[Ncpu];
    for (int p = 0; p < Ncpu; p++)
    {
      sommeNoeudsTemp2[p] = new int[numberFacesTemp[p]];
      facesTemp2[p] = new int*[numberFacesTemp[p]];
      for (int i = 0; i < numberFacesTemp[p]; i++)
      {
        facesTemp2[p][i] = new int[4];
      }
    }
    for (int i = numberElementsFrontiere; i < numberElementsGlobal; i++)
    {
      if ((i - numberElementsFrontiere) % frequenceImpressure == 0) { cout << "    " << (100 * (i - numberElementsFrontiere) / (numberElementsGlobal - numberElementsFrontiere)) << "% ... " << endl; }
      int numCPU(elementsGlobal[i]->getCPU());
      elementsGlobal[i]->construitFacesSimplifie(iMaxFaces[numCPU], facesTemp2[numCPU], sommeNoeudsTemp2[numCPU]);
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK in " << t1 << " seconds" << endl;

    //5)Recherche des elements fantomes qui ont une face commune avec un element interne (communicants)
    //-------------------------------------------------------------------------------------------------
    cout << "  5/Looking for ghosts elements ..." << endl;
    tTemp = clock();
    frequenceImpressure = max(numberElementsGlobal / 10, 1);
    int *numberFacesCommunicantesCPU = new int[Ncpu];
    for (int i = 0; i < Ncpu; i++) { numberFacesCommunicantesCPU[i] = 0; }
    for (int i = 0; i < numberElementsGlobal; i++)
    {
      if (i%frequenceImpressure == 0) { cout << "    " << (100 * i / numberElementsGlobal) << "% ... " << endl; }
      if (elementsGlobal[i]->getNumberAutresCPU() != 0)
      {
        vector<int> CPUAEnlever;
        for (int p = 0; p < elementsGlobal[i]->getNumberAutresCPU(); p++)
        {
          int numCPU(elementsGlobal[i]->getAutreCPU(p));
          //verification si l element est communicant via une face
          //******************************************************
          int numberNoeuds(elementsGlobal[i]->getNumberNoeuds());
          bool *isNoeudInterne = new bool[numberNoeuds];
          for (int n = 0; n < numberNoeuds; n++)
          {
            isNoeudInterne[n] = false;
            noeudCourant = elementsGlobal[i]->getNumNoeud(n);
            for (int j = 0; j < numberNoeudsInterne[numCPU]; j++)
            {
              if (noeudsCPU[numCPU][j] == noeudCourant)
              {
                isNoeudInterne[n] = true;
                break;
              }
            }
          }
          //Determination du number de faces communicantes
          //**********************************************
          //int numberFacesCommunicantes(elementsGlobal[i]->compteFaceCommunicante(facesTemp[numCPU],sommeNoeudsTemp[numCPU]));
          int numberFacesCommunicantes(elementsGlobal[i]->compteFaceCommunicante(iMaxFaces[numCPU], facesTemp2[numCPU], sommeNoeudsTemp2[numCPU]));
          if (numberFacesCommunicantes > 0)
          {
            //L'element est communicant, on l'ajoute ainsi que ces noeuds
            elementsCPU[numCPU].push_back(i);
            for (int n = 0; n < numberNoeuds; n++)
            {
              if (!isNoeudInterne[n]) { noeudsCPU[numCPU].push_back(elementsGlobal[i]->getNumNoeud(n)); }
            }
            numberFacesCommunicantesCPU[numCPU] += numberFacesCommunicantes;
          }
          else
          {
            CPUAEnlever.push_back(p);
          }
          delete[] isNoeudInterne;

        } //Fin CPU
        elementsGlobal[i]->enleveCPUAutres(CPUAEnlever);
      }
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK in " << t1 << " seconds" << endl;

    for (int p = 0; p < Ncpu; p++)
    {
      for (int i = 0; i < numberFacesTemp[p]; i++) { delete facesTemp2[p][i]; }
      delete sommeNoeudsTemp2[p]; delete facesTemp2[p];
    }
    delete[]facesTemp2; delete[]sommeNoeudsTemp2;
    delete[]numberFacesTemp; delete[]iMaxFaces;

    //6) Ecriture des fichiers de mesh pour chaque CPU
    cout << "  6/Printing mesh files for each of " << Ncpu << " CPU ..." << endl;
    tTemp = clock();
    for (int p = 0; p < Ncpu; p++)
    {
      stringstream flux;
      flux << p;
      string fichierMeshCPU("./libMeshes/" + m_nameMesh + "_CPU" + flux.str() + ".msh");
      cout << "    print '" << fichierMeshCPU.c_str() << "' in progress ...' " << endl;
      ofstream fileStream;
      fileStream.open(fichierMeshCPU.c_str());

      fileStream << "$MeshFormat" << endl;
      fileStream << "2.2 0 8" << endl;
      fileStream << "$EndMeshFormat" << endl;
      fileStream << "$Nodes" << endl;
      //Reordonne les noeuds
      //sort(noeudsCPU[p].begin(), noeudsCPU[p].end());
      fileStream << noeudsCPU[p].size() << endl;
      for (unsigned int i = 0; i < noeudsCPU[p].size(); i++)
      {
        fileStream << i + 1 << " " << noeudsGlobal[noeudsCPU[p][i]].getX() << " " << noeudsGlobal[noeudsCPU[p][i]].getY() << " " << noeudsGlobal[noeudsCPU[p][i]].getZ() << endl;
      }
      fileStream << "$EndNodes" << endl;
      fileStream << "$Elements" << endl;
      fileStream << elementsCPU[p].size() << endl;
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
          //Renumbering locale
          noeudCourant = elementsGlobal[e]->getNumNoeud(n);
          for (int j = 0; j < static_cast<int>(noeudsCPU[p].size()); j++)
          {
            if (noeudsCPU[p][j] == noeudCourant)
            {
              fileStream << " " << j + 1; break;
            }
          }
        }
        fileStream << endl;
      }
      fileStream << "$EndElements" << endl;
      //Informations additionelles utiles !!
      fileStream << "Info non lue par Gmsh : number de faces communicante" << endl;
      fileStream << numberFacesCommunicantesCPU[p] << endl;
      fileStream << "Info non lue par Gmsh : number de noeuds internes (hors fantomes)" << endl;
      fileStream << numberNoeudsInterne[p] << endl;
      fileStream.close();
    }
    tTemp = clock() - tTemp; t1 = static_cast<float>(tTemp) / CLOCKS_PER_SEC;
    cout << "    OK in " << t1 << " seconds" << endl;

    cout << "... MESH FILE PRETRAITEMENT COMPLETE ";
    totalTime = clock() - totalTime; t1 = static_cast<float>(totalTime) / CLOCKS_PER_SEC;
    cout << " Total time of pretraitement : " << t1 << " seconds" << endl;

    cout << "------------------------------------------------------" << endl;

    //Desallocations
    for (int i = 0; i < numberElementsGlobal; i++)
    {
      delete elementsGlobal[i];
    }
    delete[] elementsGlobal;
    delete[] noeudsGlobal;
    delete[] numberNoeudsInterne;
    delete[] numberFacesCommunicantesCPU;

  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MeshUnStruct::lectureGeometrieGmsh(vector<ElementNS*>** voisinsNoeuds)
{
  try {
    //Opening mesh file at Gmsh format
    //-----------------------------------
    cout << "------------------------------------------------------" << endl;
    cout << " A) READING MESH FILE " + m_fichierMesh + " IN PROGRESS..." << endl;
    m_fichierMesh = "./libMeshes/" + m_fichierMesh;
    ifstream fichierMesh(m_fichierMesh.c_str(), ios::in);
    if (!fichierMesh){ throw ErrorECOGEN("file mesh absent :" + m_fichierMesh, __FILE__, __LINE__); }
    string currentLine;
    getline(fichierMesh, currentLine);

    //Looking for version
    //-------------------
    stringstream lineToTreat;
    double fileVersion;
    getline(fichierMesh, currentLine);
    lineToTreat << currentLine;
    lineToTreat >> fileVersion;

    cout << "  MeshFile format " << fileVersion << endl;
    if (fileVersion < 4.) { this->readGmshV2(voisinsNoeuds, fichierMesh); }
    else { 
      throw ErrorECOGEN("MeshFile version 4 not yet implemented for Gmsh. Please try with v2", __FILE__, __LINE__);
      this->readGmshV4(voisinsNoeuds, fichierMesh);
    }

    //Closing file and information printing
    //-------------------------------------
    fichierMesh.close();
    cout << "OK" << endl;
    cout << endl << "  --------------------------" << endl;
    cout << "    MESH INFORMATIONS :" << endl;
    cout << "  --------------------------" << endl;
    cout << "    mesh nodes number : " << m_numberNoeuds << endl;
    cout << "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    if (m_numberElements0D != 0)
    {
      cout << endl << "    0D elements number : " << m_numberElements0D << endl;
      cout << "    ~~~~~~~~~~~~~~~~~~~~" << endl;
      if (m_numberPoints != 0) { cout << "      - number vertex : " << m_numberPoints << endl; }
    }
    if (m_numberElements1D != 0)
    {
      cout << endl << "    1D elements number : " << m_numberElements1D << endl;
      cout << "    ~~~~~~~~~~~~~~~~~~~~" << endl;
      if (m_numberSegments != 0) { cout << "      - number segments : " << m_numberSegments << endl; }
    }
    if (m_numberElements2D != 0)
    {
      cout << endl << "    2D elements number : " << m_numberElements2D << endl;
      cout << "    ~~~~~~~~~~~~~~~~~~~~" << endl;
      if (m_numberTriangles != 0) { cout << "      - number triangles   : " << m_numberTriangles << endl; }
      if (m_numberQuadrangles != 0) { cout << "      - number quadrangles : " << m_numberQuadrangles << endl; }
      cout << "      Total surface : " << m_totalSurface << " m2" << endl;
    }
    if (m_numberElements3D != 0)
    {
      cout << endl << "    3D elements number : " << m_numberElements3D << endl;
      cout << "    ~~~~~~~~~~~~~~~~~~~~" << endl;
      if (m_numberTetrahedrons != 0) { cout << "      - number tetraedres   : " << m_numberTetrahedrons << endl; }
      if (m_numberPyramids != 0) { cout << "      - number pyramides    : " << m_numberPyramids << endl; }
      if (m_numberHexahedrons != 0) { cout << "      - number hexaedres    : " << m_numberHexahedrons << endl; }
      cout << "      Total volume : " << m_totalVolume << " m3" << endl;
    }
    cout << endl << "... READING MESH FILE COMPLETE " << endl;
    cout << "------------------------------------------------------" << endl;
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MeshUnStruct::readGmshV2(vector<ElementNS*>** voisinsNoeuds, ifstream &meshFile)
{
  try {
    string currentLine;
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);

    //1) Filling m_noeuds array
    //-------------------------
    cout << "  1/Mesh nodes reading ...";
    meshFile >> m_numberNoeuds;
    meshFile.ignore(1, '\n');
    m_noeuds = new Coord[m_numberNoeuds];
    *voisinsNoeuds = new vector<ElementNS*>[m_numberNoeuds];
    int inutile(0); double x, y, z;
    for (int i = 0; i < m_numberNoeuds; i++)
    {
      meshFile >> inutile >> x >> y >> z;
      m_noeuds[i].setXYZ(x, y, z);
    }
    meshFile.ignore(1, '\n');
    getline(meshFile, currentLine);
    getline(meshFile, currentLine);
    cout << "OK" << endl;

    //2) 1D/2D/3D elements are stored in m_elements array / counting
    //--------------------------------------------------------------
    cout << "  2/0D/1D/2D/3D elements reading ...";
    meshFile >> m_numberElements;
    meshFile.ignore(1, '\n');
    //Allocation tableau d elements
    m_elements = new ElementNS*[m_numberElements];
    //Lecture des elements et attributions proprietes geometriques
    m_numberElements1D = 0, m_numberElements2D = 0, m_numberElements3D = 0;
    int posDebutElement = static_cast<int>(meshFile.tellg()); //reperage debut des elements pour retour rapide
    int noeudG;
    for (int i = 0; i < m_numberElements; i++) {
      this->lectureElementGmshV2(m_noeuds, meshFile, &m_elements[i]);
      if (m_elements[i]->getTypeGmsh() == 15) { m_numberElements0D++; }
      else if (m_elements[i]->getTypeGmsh() == 1) { m_numberElements1D++; }
      else if (m_elements[i]->getTypeGmsh() <= 3) { m_numberElements2D++; m_totalSurface += m_elements[i]->getVolume(); }
      else if (m_elements[i]->getTypeGmsh() <= 7) { m_numberElements3D++; m_totalVolume += m_elements[i]->getVolume(); }
      else { throw ErrorECOGEN("Type element du .msh non gere dans ECOGEN", __FILE__, __LINE__); }
      //Attribution element i voisin pour les noeuds concernes (Ordre 2 muiltislopes)
      for (int n = 0; n < m_elements[i]->getNumberNoeuds(); n++) {
        noeudG = m_elements[i]->getNumNoeud(n);
        (*voisinsNoeuds)[noeudG].push_back(m_elements[i]);
      }
    }
    m_numberElementsInternes = m_numberElements;

    //for (int n = 0; n < m_numberNoeuds; n++) {
    //  for (int e = 0; e < (*voisinsNoeuds)[n].size(); e++) {
    //    cout << (*voisinsNoeuds)[n][e]->getIndex() << " ";
    //  }
    //  cout << endl;
    //}
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MeshUnStruct::readGmshV4(vector<ElementNS*>** voisinsNoeuds, ifstream &meshFile)
{
  try {
    string currentLine;

    //1) Filling m_noeuds array
    //-------------------------
    { 
      cout << "  1/Mesh nodes reading ...";
      while (currentLine != "$Nodes") {
        getline(meshFile, currentLine);
        if (meshFile.eof()) { throw ErrorECOGEN("Nodes block not found in mesh file", __FILE__, __LINE__); }
      }
      int numEntityBlocks;
      { stringstream lineTotreat;
      getline(meshFile, currentLine); lineTotreat << currentLine;
      lineTotreat >> numEntityBlocks >> m_numberNoeuds; }
      m_noeuds = new Coord[m_numberNoeuds];
      //Allocate node array
      *voisinsNoeuds = new vector<ElementNS*>[m_numberNoeuds];
      //Reading nodes
      int tagEntity, dimEntity, parametric, numNodesEntity, tag; double x, y, z;
      for (int e = 0; e < numEntityBlocks; e++) { //Reading entity
        {stringstream lineTotreat;
        getline(meshFile, currentLine); lineTotreat << currentLine;
        lineTotreat >> tagEntity >> dimEntity >> parametric >> numNodesEntity; }
        for (int i = 0; i < numNodesEntity; i++) {
          stringstream lineTotreat;
          getline(meshFile, currentLine); lineTotreat << currentLine;
          lineTotreat >> tag >> x >> y >> z;
          m_noeuds[i].setXYZ(x, y, z);
        }
      }
      getline(meshFile, currentLine);
      if (currentLine != "$EndNodes") { throw ErrorECOGEN("Nodes block not completely read in mesh file", __FILE__, __LINE__); }
      cout << "OK" << endl;
    }

    //2) 1D/2D/3D elements are stored in m_elements array / counting
    //--------------------------------------------------------------
    {
      cout << "  2/0D/1D/2D/3D elements reading ...";
      while (currentLine != "$Elements") {
        getline(meshFile, currentLine);
        if (meshFile.eof()) { throw ErrorECOGEN("Elements block not found in mesh file", __FILE__, __LINE__); }
      }
      int numEntityBlocks;
      { stringstream lineTotreat;
      getline(meshFile, currentLine); lineTotreat << currentLine;
      lineTotreat >> numEntityBlocks >> m_numberElements; }
      //Allocate elements array
      m_elements = new ElementNS*[m_numberElements];
      //Readin elements and geometrical properties attributions
      m_numberElements1D = 0, m_numberElements2D = 0, m_numberElements3D = 0;
      int posDebutElement = static_cast<int>(meshFile.tellg()); //reperage debut des elements pour retour rapide
      int noeudG;
      int tagEntity, dimEntity, typeEle, numElementsEntity;
      int numElement(0);
      for (int e = 0; e < numEntityBlocks; e++) { //Reading entity
        {stringstream lineTotreat;
        getline(meshFile, currentLine); lineTotreat << currentLine;
        lineTotreat >> tagEntity >> dimEntity >> typeEle >> numElementsEntity; }
        for (int i = 0; i < numElementsEntity; i++) {
          this->lectureElementGmshV4(m_noeuds, meshFile, &m_elements[numElement], typeEle,numElement,tagEntity);
          //Les tags entity ne sont pas bon ni l'ordonancement des faces et cellules => a revoir
          if (m_elements[numElement]->getTypeGmsh() == 15) { m_numberElements0D++; }
          else if (m_elements[numElement]->getTypeGmsh() == 1) { m_numberElements1D++; }
          else if (m_elements[numElement]->getTypeGmsh() <= 3) { m_numberElements2D++; m_totalSurface += m_elements[i]->getVolume(); }
          else if (m_elements[numElement]->getTypeGmsh() <= 7) { m_numberElements3D++; m_totalVolume += m_elements[i]->getVolume(); }
          else { throw ErrorECOGEN("Type element du .msh non gere dans ECOGEN", __FILE__, __LINE__); }
          //Attribution element i voisin pour les noeuds concernes (Ordre 2 muiltislopes)
          for (int n = 0; n < m_elements[numElement]->getNumberNoeuds(); n++) {
            noeudG = m_elements[numElement]->getNumNoeud(n);
            (*voisinsNoeuds)[noeudG].push_back(m_elements[numElement]);
          }
          numElement++;
        }
      }
      m_numberElementsInternes = m_numberElements;
      getline(meshFile, currentLine);
      if (currentLine != "$EndElements") { throw ErrorECOGEN("Elements block not completely read in mesh file", __FILE__, __LINE__); }
      cout << "OK" << endl;
    }

  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MeshUnStruct::lectureGeometrieGmshParallele()
{
  int numberNoeudsTotal(0), numberElementsTotal(0);
  int numberCPUVoisins(0), CPUCourant(0);

  try{
    //1) Ouverture du file de mesh
    //-------------------------------------
    if (rankCpu == 0)
    {
      cout << "------------------------------------------------------" << endl;
      cout << " A) READING MESH FILE " + m_nameMesh + "_CPUX.msh" + " IN PROGRESS ..." << endl;
    }
    stringstream flux;
    flux << rankCpu;
    m_fichierMesh = "./libMeshes/" + m_nameMesh + "_CPU" + flux.str() + ".msh";
    ifstream fichierMesh(m_fichierMesh.c_str(), ios::in);
    if (!fichierMesh) { throw ErrorECOGEN("file mesh absent :" + m_fichierMesh, __FILE__, __LINE__); }
    string ligneCourante;
    getline(fichierMesh, ligneCourante);
    getline(fichierMesh, ligneCourante);
    getline(fichierMesh, ligneCourante);
    getline(fichierMesh, ligneCourante);

    //2) Stockage de la grille de vertex dans tableau m_noeuds
    //-------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0) { cout << "  1/Reading mesh nodes ..."; }
    fichierMesh >> m_numberNoeuds;
    fichierMesh.ignore(1, '\n');
    m_noeuds = new Coord[m_numberNoeuds];
    int inutile(0); double x, y, z;
    for (int i = 0; i < m_numberNoeuds; i++)
    {
      fichierMesh >> inutile >> x >> y >> z;
      m_noeuds[i].setXYZ(x, y, z);
    }
    fichierMesh.ignore(1, '\n');
    getline(fichierMesh, ligneCourante);
    getline(fichierMesh, ligneCourante);
    if (rankCpu == 0) { cout << "OK" << endl; }

    //3) Recuperation des elements 1D/2D/3D dans le tableau m_elements et comptage
    //----------------------------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0) { cout << "  2/Reading internal 1D/2D/3D elements ..."; }
    fichierMesh >> m_numberElements;
    fichierMesh.ignore(1, '\n');
    //Allocation tableau d elements
    m_elements = new ElementNS*[m_numberElements];
    //Lecture des elements et attributions proprietes geometriques
    m_numberElements1D = 0, m_numberElements2D = 0, m_numberElements3D = 0;
    int posDebutElement = static_cast<int>(fichierMesh.tellg()); //reperage debut des elements pour retour rapide
    for (int i = 0; i < m_numberElements; i++)
    {
      this->lectureElementGmshV2(m_noeuds, fichierMesh, &m_elements[i]);
      if (m_elements[i]->getCPU() == rankCpu)
      {
        if (m_elements[i]->getTypeGmsh() == 1) { m_numberElements1D++; }
        else if (m_elements[i]->getTypeGmsh() <= 3) { m_numberElements2D++; m_totalSurface += m_elements[i]->getVolume(); }
        else if (m_elements[i]->getTypeGmsh() <= 7) { m_numberElements3D++; m_totalVolume += m_elements[i]->getVolume(); }
        else { throw ErrorECOGEN("Type element du .msh non gere dans ECOGEN", __FILE__, __LINE__); }
      }
      else { m_numberElementsFantomes++; }
    }
    fichierMesh.ignore(1, '\n');
    getline(fichierMesh, ligneCourante);
    //Lecture informations hors Gmsh
    getline(fichierMesh, ligneCourante);
    fichierMesh >> m_numberFacesParallele;
    fichierMesh.ignore(1, '\n');
    getline(fichierMesh, ligneCourante);
    fichierMesh >> m_numberNoeudsInternes;
    fichierMesh.close();
    //Calcul du number d'elements propres au CPU
    m_numberElementsInternes = m_numberElements - m_numberElementsFantomes;
    if (rankCpu == 0) { cout << "OK" << endl; }
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void MeshUnStruct::lectureElementGmshV2(const Coord *TableauNoeuds, ifstream &fichierMesh, ElementNS **element)
{
  int numberElement,numberTags,typeElement,numberEntitePhysique,numberEntiteGeometrique;
  bool creeElement(true);
  fichierMesh >> numberElement >> typeElement;
  //Lecture des tags
  fichierMesh >> numberTags;
  fichierMesh >> numberEntitePhysique;    //Appartenance physique
  fichierMesh >> numberEntiteGeometrique; //Appartenance geometrique

  //1)Affectation du number de vertex selon element
  //----------------------------------------------
  switch (typeElement)
  {
    case 1: //segment (deux points)
      *element = new ElementSegment;// cout << "segment trouve" << endl;
      m_numberSegments++;
      break;
    case 2: //triangle (trois points)
      *element = new ElementTriangle;// cout << "triangle trouve" << endl;
      m_numberTriangles++;
      break;
    case 3: //Quadrangle (quatre points)
      *element = new ElementQuadrangle;// cout << "quadrangle trouve" << endl;
      m_numberQuadrangles++;
      break;
    case 4: //Tetrahedron (quatre points)
      *element = new ElementTetrahedron;// cout << "tetraedre trouve" << endl;
      m_numberTetrahedrons++;
      break;
    //case 7: //Pyramid quadrangulaire (cinq points) // Ce type d'element semble ne pas fonctionner avec GMSH, les volumes des elements du mesh semblent poser probleme...
    //  *element = new ElementPyramid;
    //  m_numberPyramids++;
    //  break;
    case 15: //Point (un vertex)
      *element = new ElementPoint; // cout << "vertex trouve" << endl;
      m_numberPoints++;
      break;
    case 5: //Hexahedron (huit points)
      *element = new ElementHexahedron; // cout << "Hexahedron trouve" << endl;
      m_numberHexahedrons++;
      break;
    case 6: //Prism (six points)
      *element = new ElementPrism; // cout << "Hexahedron trouve" << endl;
      m_numberHexahedrons++;
      break;
    default:
      Errors::errorMessage("Type d element du file .msh inconnu de ECOGEN");
      break;
  } //Fin switch typeElement
 
  //2) Specificite meshs paralleles
  //-----------------------------------
  int numberCPU(0);
  if (numberTags > 2)
  {
    fichierMesh >> numberCPU; //number de partition de mesh auquel appartient l element
    int *numCPU = new int[numberCPU];
    for (int tag = 0; tag < numberCPU; tag++){ fichierMesh >> numCPU[tag]; }
    (*element)->setAppartenanceCPU(numCPU, numberCPU);
    delete[] numCPU;
  }

  //3) Construction de l'element et de ses proprietes
  //-------------------------------------------------
  int noeudCourant;
  int *numNoeud = new int[(*element)->getNumberNoeuds()];
  Coord *noeud = new Coord[(*element)->getNumberNoeuds()];
  for (int i = 0; i < (*element)->getNumberNoeuds(); i++)
  {
    fichierMesh >> noeudCourant;
    numNoeud[i] = noeudCourant-1;         //decalage car tableau commencant a zero
    noeud[i] = TableauNoeuds[noeudCourant - 1];
  }
  int indexElement(numberElement - 1);
  (*element)->construitElement(numNoeud, noeud, numberEntitePhysique, numberEntiteGeometrique, indexElement);

  delete[] noeud;
  delete[] numNoeud;
  
}

//***********************************************************************

void MeshUnStruct::lectureElementGmshV4(const Coord *TableauNoeuds, ifstream &fichierMesh, ElementNS **element, const int &typeElement, int &indiceElement, const int & physicalEntity)
{
  try {
    //1)Number of vertex affectation
    //------------------------------
    switch (typeElement)
    {
    case 1: //segment (deux points)
      *element = new ElementSegment;
      m_numberSegments++;
      break;
    case 2: //triangle (trois points)
      *element = new ElementTriangle;
      m_numberTriangles++;
      break;
    case 3: //Quadrangle (quatre points)
      *element = new ElementQuadrangle;
      m_numberQuadrangles++;
      break;
    case 4: //Tetrahedron (quatre points)
      *element = new ElementTetrahedron;
      m_numberTetrahedrons++;
      break;
    case 7: //Pyramid quadrangulaire (cinq points)
      *element = new ElementPyramid;
      m_numberPyramids++;
      break;
    case 15: //Point (un vertex)
      *element = new ElementPoint;
      m_numberPoints++;
      break;
    case 5: //Hexahedron (huit points)
      *element = new ElementHexahedron;
      m_numberHexahedrons++;
      break;
    case 6: //Prism (six points)
      *element = new ElementPrism;
      m_numberHexahedrons++;
      break;
    default:
      throw ErrorECOGEN("Element type unknown in mesh file", __FILE__, __LINE__);
      break;
    } //Fin switch typeElement

    //2) Element building / properties filling
    //----------------------------------------
    int noeudCourant, tag;
    int *numNoeud = new int[(*element)->getNumberNoeuds()];
    Coord *noeud = new Coord[(*element)->getNumberNoeuds()];
    string currentLine;
    stringstream lineTotreat;
    getline(fichierMesh, currentLine); lineTotreat << currentLine;
    lineTotreat >> tag;
    for (int i = 0; i < (*element)->getNumberNoeuds(); i++)
    {
      lineTotreat >> noeudCourant;
      numNoeud[i] = noeudCourant - 1;         //decalage car tableau commencant a zero
      noeud[i] = TableauNoeuds[noeudCourant - 1];
    }
    (*element)->construitElement(numNoeud, noeud, physicalEntity, 0, indiceElement);

    delete[] noeud;
    delete[] numNoeud;
  }
  catch (ErrorECOGEN &) { throw; }

}

//**************************************************************************
//******************************** ECRITURE ********************************
//**************************************************************************

void MeshUnStruct::ecritHeaderPiece(std::ofstream &fileStream, std::vector<Cell *> *cellsLvl, int lvl) const
{
  fileStream << "    <Piece NumberOfPoints=\"" << m_numberNoeuds << "\" NumberOfCells=\"" << m_numberCellsCalcul - m_numberCellsFantomes << "\">" << endl;
}

//****************************************************************************

void MeshUnStruct::recupereNoeuds(std::vector<double> &jeuDonnees, int lvl) const
{
  for (int noeud = 0; noeud < m_numberNoeuds; noeud++)
  {
    jeuDonnees.push_back(m_noeuds[noeud].getX());
    jeuDonnees.push_back(m_noeuds[noeud].getY());
    jeuDonnees.push_back(m_noeuds[noeud].getZ());
  }
}

//****************************************************************************

void MeshUnStruct::recupereConnectivite(std::vector<double> &jeuDonnees, int lvl) const
{
  for (int i = m_numberFacesLimites; i < m_numberElementsInternes; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      for (int noeud = 0; noeud < m_elements[i]->getNumberNoeuds(); noeud++)
      {
        jeuDonnees.push_back(m_elements[i]->getNumNoeud(noeud));
      }
    }
  }
}

//****************************************************************************

void MeshUnStruct::recupereOffsets(std::vector<double> &jeuDonnees, int lvl) const
{
  int offset(0);
  for (int i = m_numberFacesLimites; i < m_numberElementsInternes; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      offset += m_elements[i]->getNumberNoeuds();
      jeuDonnees.push_back(offset);
    }
  }
}

//****************************************************************************

void MeshUnStruct::recupereTypeCell(std::vector<double> &jeuDonnees, int lvl) const
{
  for (int i = m_numberFacesLimites; i < m_numberElementsInternes; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      jeuDonnees.push_back(m_elements[i]->getTypeVTK());
    }
  }
}

//****************************************************************************

void MeshUnStruct::recupereDonnees(vector<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, const int var, int phase, int lvl) const
{
  jeuDonnees.clear();
  int numCell;
  for (int i = m_numberFacesLimites; i < m_numberElementsInternes; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCellAssociee();
      if (var > 0) { //On veut recuperer les donnees scalars
        if (phase >= 0) { jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnScalar(var)); }      //Donnees de phases
        else if (phase == -1) { jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnScalar(var)); }   //Donnees de mixture
        else if (phase == -2) { jeuDonnees.push_back(cellsLvl[0][numCell]->getTransport(var-1).getValue()); }
        else if (phase == -3) { jeuDonnees.push_back(cellsLvl[0][numCell]->getXi()); }
        else if (phase == -4) { jeuDonnees.push_back(cellsLvl[0][numCell]->getGradient()); }
        else { Errors::errorMessage("MeshUnStruct::recupereDonnees: unknown number of phase: ", phase); }
      }
      else { //On veut recuperer les donnees vectorielles
        if (phase >= 0) { //Phases data
          jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getX());
          jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getY());
          jeuDonnees.push_back(cellsLvl[0][numCell]->getPhase(phase)->returnVector(-var).getZ());
        }
        else if(phase == -1){  //Mixture data
          jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getX());
          jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getY());
          jeuDonnees.push_back(cellsLvl[0][numCell]->getMixture()->returnVector(-var).getZ());
        }
        else { Errors::errorMessage("MeshUnStruct::recupereDonnees: unknown number of phase: ", phase); }
      } //Fin vecteur
    }
  }
}

//****************************************************************************

void MeshUnStruct::setDataSet(std::vector<double> &jeuDonnees, vector<Cell *> *cellsLvl, const int var, int phase, int lvl) const
{
  int iterDataSet(0);
  int numCell;
  Coord vec;
  for (int i = m_numberFacesLimites; i < m_numberElementsInternes; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCellAssociee();
      if (var > 0) { //Scalars data are first set
        if (phase >= 0) { cellsLvl[0][numCell]->getPhase(phase)->setScalar(var, jeuDonnees[iterDataSet++]); } //phases data
        else if (phase == -1) { cellsLvl[0][numCell]->getMixture()->setScalar(var, jeuDonnees[iterDataSet++]); }  //mixture data
        else if (phase == -2) { cellsLvl[0][numCell]->getTransport(var - 1).setValue(jeuDonnees[iterDataSet++]); } //transport data
        else { Errors::errorMessage("MeshUnStruct::setDataSet: unknown phase number: ", phase); }
      }
      else { //On veut recuperer les donnees vectorielles
        if (phase >= 0) { //Phases data
          vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet+1], jeuDonnees[iterDataSet+2]);
          cellsLvl[0][numCell]->getPhase(phase)->setVector(-var, vec);
          iterDataSet += 3;
        }
        else if (phase == -1) {  //Mixture data
          vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet+1], jeuDonnees[iterDataSet+2]);
          cellsLvl[0][numCell]->getMixture()->setVector(-var, vec);
          iterDataSet += 3;
        }
        else { Errors::errorMessage("MeshUnStruct::setDataSet: unknown phase number: ", phase); }
      } //Fin vecteur
    }
  }
}

//****************************************************************************
void MeshUnStruct::extractAbsVeloxityMRF(vector<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, Source *sourceMRF, int lvl) const
{
  jeuDonnees.clear();
  int numCell;
  for (int i = m_numberFacesLimites; i < m_numberElementsInternes; i++)
  {
    if (!m_elements[i]->isFantome())
    {
      numCell = m_elements[i]->getNumCellAssociee();
      Coord absoluteVelocity = sourceMRF->computeAbsVelocity(cellsLvl[0][numCell]->getVelocity(), cellsLvl[0][numCell]->getPosition());
      jeuDonnees.push_back(absoluteVelocity.getX());
      jeuDonnees.push_back(absoluteVelocity.getY());
      jeuDonnees.push_back(absoluteVelocity.getZ());
    }
  }
}

//***********************************************************************