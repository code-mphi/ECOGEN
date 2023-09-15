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

#include "FaceNS.h"

//***********************************************************************

FaceNS::FaceNS(){}

//***********************************************************************

FaceNS::FaceNS(const int& numberNodes) :
Face(),
m_numberNodes(numberNodes),
m_limite(false),
m_comm(false),
m_elementGauche(0),
m_elementDroite(0)
{
  m_numNodes = new int[numberNodes];
}

//***********************************************************************

FaceNS::~FaceNS()
{
  delete[] m_numNodes;
}

//***********************************************************************

ElementNS *FaceNS::getElementGauche() const
{
  return m_elementGauche;
}

//***********************************************************************

ElementNS *FaceNS::getElementDroite() const
{
  return m_elementDroite;
}

//***********************************************************************

void FaceNS::construitFace(const Coord* nodes, const int& numNodeOther, ElementNS *elementNeighbor)
{

  //Calcul position du center de face
  m_position = 0.;
  for (int i = 0; i < m_numberNodes; i++){ m_position += nodes[m_numNodes[i]]; }
  m_position /= static_cast<double>(m_numberNodes);
  //Calculs des proprietes de la face
  this->computeSurface(nodes);
  this->computeRepere(nodes, numNodeOther, elementNeighbor);
}

//***********************************************************************

bool FaceNS::faceExists(FaceNS** faces, const int& indexMaxFaces, int& indexFaceExiste)
{
  int faceTrouvee(1);
  //for (int i = 0; i < indexMaxFaces; i++)
  for (int i = indexMaxFaces-1; i>=0; i--)
  {
    //1) Sans passer par operateur de comparaison
    //-------------------------------------------
    if (faces[i]->m_sumNumNodes != m_sumNumNodes){continue;};
    //Verification node par node
    faceTrouvee = 1;
    for (int n = 0; n < m_numberNodes; n++)
    {
      if (m_numNodes[n] != faces[i]->getNumNode(n)){ faceTrouvee = 0; break; }
    }
    if(faceTrouvee)
    {
      indexFaceExiste = i; return true;
    }

    // //2) Avec operateur de comparaison
    // //--------------------------------
    // if (*faces[i] == *this){ indexFaceExiste = i; return true; }
  }
  indexFaceExiste = 0;
  return false;
}

//***********************************************************************

void FaceNS::addElementNeighbor(ElementNS *elementNeighbor)
{
  if (m_elementGauche == 0)
  {
    m_elementGauche = elementNeighbor;
  }
  else
  {
    m_elementDroite = elementNeighbor;
  }
}

//***********************************************************************

void FaceNS::addElementNeighborLimite(ElementNS *elementNeighbor)
{
  if (m_elementGauche == 0)
  {
    //On inverse pour avoir l element de compute a Gauche (Necessaire pour les limites)
    m_elementGauche = m_elementDroite;
    m_normal.changeSign();
    m_binormal = Coord::crossProduct(m_normal, m_tangent);
  }
  m_elementDroite = elementNeighbor;
  m_limite = true;
}

//***********************************************************************

void FaceNS::setEstLimite(const bool &estLimite)
{
  m_limite = estLimite;
}

//***********************************************************************

void FaceNS::setEstComm(const bool &estComm)
{
  m_comm = estComm;
}

//***********************************************************************

void FaceNS::getInfoNodes(int* numNodes, int& sumNumNodes) const
{
  for(int i=0; i<m_numberNodes; i++)
  {
    numNodes[i] = m_numNodes[i];
  }
  sumNumNodes = m_sumNumNodes;
}

//***********************************************************************

void FaceNS::printInfo() const
{
  std::cout << "----------------" << std::endl;
  //cout << "Face" << endl << " Nodes : ";
  //for (int i = 0; i < m_numberNodes; i++)
  //{
  //  cout << m_numNodes[i] << " ";
  //}
  //cout << endl;
  std::cout << "center : " << m_position.getX() << " " << m_position.getY() << " " << m_position.getZ() << std::endl;
  //cout << " limites : " << m_limite << ", comm : " << m_comm << endl;
  //if (m_elementGauche != 0)
  //{
  //  cout << " Element gauche " << m_elementGauche->getindex();
  //}
  //if (m_elementDroite != 0)
  //{
  //  cout << " Element Droite " << m_elementDroite->getindex();
  //}
  //cout << endl;

  //cout << " normal : "; m_normal.printInfo();
  //cout << " tangent : "; m_tangent.printInfo();
  //cout << " binormal : "; m_binormal.printInfo();
}

//***********************************************************************

void FaceNS::printNodes() const
{
  std::cout << "Face" << " Nodes : ";
  for (int i = 0; i < m_numberNodes; i++)
  {
    std::cout << m_numNodes[i] << " ";
  }
  std::cout << std::endl;
}

//*********************************************************************

int FaceNS::searchFace(int* face, int& sumNodes, int** arrayFaces, int* arraySumNodes, int numberNodes, int& indexMaxFaces)
{
    for (int i = indexMaxFaces-1; i>=0; i--)
    {
      int existe = 1;
      //On utilise la sum des number des nodes pour faire un premier tri : permet d accelerer la recherche
      if (arraySumNodes[i] != sumNodes){ continue; };
      //Verification node par node
      for (int n=0; n<numberNodes; n++)
      {
        if (face[n] != arrayFaces[i][n]){ existe = 0; break; }
      }
      if (existe) return i; //Face trouvee, on renvoi son number
    }
    return -1;  //Face non trouvee
}

//*********************************************************************

int FaceNS::searchFace(int* face, int& sumNodes, std::vector<int*> arrayFaces, std::vector<int> arraySumNodes, int numberNodes, int& indexMaxFaces)
{
    for (int i = indexMaxFaces-1; i>=0; i--)
    {
      int existe = 1;
      //On utilise la sum des number des nodes pour faire un premier tri : permet d accelerer la recherche
      if (arraySumNodes[i] != sumNodes){ continue; };
      //Verification node par node
      for (int n=0; n<numberNodes; n++)
      {
        if (face[n] != arrayFaces[i][n]){ existe = 0; break; }
      }
      if (existe) return i; //Face trouvee, on renvoi son number
    }
        // cout << "AHAHAHAHAHA nouveau" << endl;
    return -1;  //Face non trouvee
}

//*********************************************************************
//Surcharge operateur externe a la classe car prends deux arguments

bool operator==(const FaceNS &a, const FaceNS &b)
{
  //Sortie directe si m_sumNumNodes differents
  if (a.getSumNumNodes() != b.getSumNumNodes()){ return false; }
  //Sortie directe si number nodes differents
  if (a.getNumberNodes() != b.getNumberNodes()){ return false; }
  //Verification node par node
  for (int i = 0; i < a.getNumberNodes(); i++)
  {
    if (a.getNumNode(i) != b.getNumNode(i)){ return false; }
  }
  //Par defaut : il sont identiques
  return true;
}

//*********************************************************************

bool operator!=(const FaceNS &a, const FaceNS &b)
{
  return !(a == b);
}

//*********************************************************************