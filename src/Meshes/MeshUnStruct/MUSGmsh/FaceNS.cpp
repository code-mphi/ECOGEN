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

//! \file      FaceNS.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "FaceNS.h"

//***********************************************************************

FaceNS::FaceNS(){}

//***********************************************************************

FaceNS::FaceNS(const int &numberNoeuds) :
Face(),
m_numberNoeuds(numberNoeuds),
m_limite(false),
m_comm(false),
m_elementGauche(0),
m_elementDroite(0)
{
  m_numNoeuds = new int[numberNoeuds];
}

//***********************************************************************

FaceNS::~FaceNS()
{
  delete[] m_numNoeuds;
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

void FaceNS::construitFace(const Coord *noeuds, const int &numNoeudAutre, ElementNS *elementVoisin)
{

  //Calcul position du center de face
  m_position = 0.;
  for (int i = 0; i < m_numberNoeuds; i++){ m_position += noeuds[m_numNoeuds[i]]; }
  m_position /= static_cast<double>(m_numberNoeuds);
  //Calculs des proprietes de la face
  this->computeSurface(noeuds);
  this->computeRepere(noeuds, numNoeudAutre, elementVoisin);
}

//***********************************************************************

bool FaceNS::faceExiste(FaceNS **faces, const int &indexMaxFaces, int &indexFaceExiste)
{
  int faceTrouvee(1);
  //for (int i = 0; i < indexMaxFaces; i++)
  for (int i = indexMaxFaces-1; i>=0; i--)
  {
    //1) Sans passer par operateur de comparaison
    //-------------------------------------------
    if (faces[i]->m_sommeNumNoeuds != m_sommeNumNoeuds){continue;};
    //Verification noeud par noeud
    faceTrouvee = 1;
    for (int n = 0; n < m_numberNoeuds; n++)
    {
      if (m_numNoeuds[n] != faces[i]->getNumNoeud(n)){ faceTrouvee = 0; break; }
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

void FaceNS::ajouteElementVoisin(ElementNS *elementVoisin)
{
  if (m_elementGauche == 0)
  {
    m_elementGauche = elementVoisin;
  }
  else
  {
    m_elementDroite = elementVoisin;
  }
}

//***********************************************************************

void FaceNS::ajouteElementVoisinLimite(ElementNS *elementVoisin)
{
  if (m_elementGauche == 0)
  {
    //On inverse pour avoir l element de compute a Gauche (Necessaire pour les limites)
    m_elementGauche = m_elementDroite;
    m_normal.changeSign();
    m_binormal = Coord::crossProduct(m_normal, m_tangent);
  }
  m_elementDroite = elementVoisin;
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

void FaceNS::getInfoNoeuds(int *numNoeuds, int &sommeNumNoeuds) const
{
  for(int i=0; i<m_numberNoeuds; i++)
  {
    numNoeuds[i] = m_numNoeuds[i];
  }
  sommeNumNoeuds = m_sommeNumNoeuds;
}

//***********************************************************************

void FaceNS::printInfo() const
{
  std::cout << "----------------" << std::endl;
  //cout << "Face" << endl << " Noeuds : ";
  //for (int i = 0; i < m_numberNoeuds; i++)
  //{
  //  cout << m_numNoeuds[i] << " ";
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

void FaceNS::afficheNoeuds() const
{
  std::cout << "Face" << " Noeuds : ";
  for (int i = 0; i < m_numberNoeuds; i++)
  {
    std::cout << m_numNoeuds[i] << " ";
  }
  std::cout << std::endl;
}

//*********************************************************************

int FaceNS::rechercheFace(int *face, int &sommeNoeuds, int **tableauFaces, int *tableauSommeNoeuds, int numberNoeuds, int &indexMaxFaces)
{
    for (int i = indexMaxFaces-1; i>=0; i--)
    {
      int existe = 1;
      //On utilise la somme des number des noeuds pour faire un premier tri : permet d accelerer la recherche
      if (tableauSommeNoeuds[i] != sommeNoeuds){ continue; };
      //Verification noeud par noeud
      for (int n=0; n<numberNoeuds; n++)
      {
        if (face[n] != tableauFaces[i][n]){ existe = 0; break; }
      }
      if (existe) return i; //Face trouvee, on renvoi son number
    }
    return -1;  //Face non trouvee
}

//*********************************************************************

int FaceNS::rechercheFace(int *face, int &sommeNoeuds, std::vector<int*> tableauFaces, std::vector<int> tableauSommeNoeuds, int numberNoeuds, int &indexMaxFaces)
{
    for (int i = indexMaxFaces-1; i>=0; i--)
    {
      int existe = 1;
      //On utilise la somme des number des noeuds pour faire un premier tri : permet d accelerer la recherche
      if (tableauSommeNoeuds[i] != sommeNoeuds){ continue; };
      //Verification noeud par noeud
      for (int n=0; n<numberNoeuds; n++)
      {
        if (face[n] != tableauFaces[i][n]){ existe = 0; break; }
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
  //Sortie directe si m_sommeNumNoeuds differents
  if (a.getSommeNumNoeuds() != b.getSommeNumNoeuds()){ return false; }
  //Sortie directe si number noeuds differents
  if (a.getNumberNoeuds() != b.getNumberNoeuds()){ return false; }
  //Verification noeud par noeud
  for (int i = 0; i < a.getNumberNoeuds(); i++)
  {
    if (a.getNumNoeud(i) != b.getNumNoeud(i)){ return false; }
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