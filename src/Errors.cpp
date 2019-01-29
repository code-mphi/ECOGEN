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

//! \file      Errors.cpp
//! \author    F. Petitpas, S. Le Martelot
//! \version   1.0
//! \date      December 20 2017

#include "Errors.h"
#include "Run.h"

using namespace std;

std::vector<Errors> errors;

//***********************************************************************

Errors::Errors() : m_state(0), m_ligne(0), m_value(0.), m_message("non renseigne")
{}

//***********************************************************************

Errors::Errors(const string &message, const char* fichierSource, int numberLigne) :
  m_message(message), m_fichier(fichierSource), m_ligne(numberLigne), m_value(0.)
{
  m_state = 1;
}

//***********************************************************************

Errors::~Errors(){}

//***********************************************************************

void Errors::errorMessage(const string &message)
{
  if (rankCpu == 0)
  { 
    stringstream numCPU;
    numCPU << rankCpu;

    cerr << endl << "-------------------------------------------------------" << endl;
    cerr << "ERREUR sur CPU " + numCPU.str() + " : " << message.c_str() << endl; 
    cerr << "Fix it and try again ;-)" << endl;
    cerr << "-------------------------------------------------------" << endl;
  }
  //run.finalize();
  exit(0);
}

//***********************************************************************

void Errors::errorMessage(const std::string &message, double value)
{
  if (rankCpu == 0)
  {
    stringstream numCPU;
    numCPU << rankCpu;

    cout << endl << "-------------------------------------------------------" << endl;
    cout << "ERREUR sur CPU " + numCPU.str() + " : " << message.c_str() << " " << value << endl;
    cout << "Corriger error et relancer le code" << endl;
    cout << "-------------------------------------------------------" << endl;
  }
  //run.finalize();
  exit(0);
}

//***********************************************************************

void Errors::setError(const std::string &message, const char* fichierSource, int numberLigne)
{
  m_message = message;
  m_state = 1;
  m_fichier = fichierSource;
  m_ligne = numberLigne;
}

//***********************************************************************

void Errors::setError(const std::string &message, const double value)
{
  m_message = message;
  m_state = 1;
  m_value = value;
}

//***********************************************************************

int Errors::getEtat()
{
  return m_state;
}

//***********************************************************************

void Errors::afficheError()
{
  cout << endl << "-------------------------------------------------------" << endl;
  cout << "        ERREUR NECESSITANT ARRET DU PROGRAMME" << endl;
  cout << " - number du CPU : " << rankCpu << endl;
  cout << " - file :  " << m_fichier.c_str() << " ligne : " << m_ligne << endl;
  cout << " - message : " << m_message.c_str() << " " << m_value << endl;
  cout << "=> Corriger error et relancer le run" << endl;
  cout << "-------------------------------------------------------" << endl;
}

//***********************************************************************

void Errors::ecritErrorFichier()
{
  ofstream fileStream;
  stringstream num;

  num << "error_CPU" << rankCpu;
  fileStream.open((num.str()).c_str());

  fileStream << "-------------------------------------------------------" << endl;
  fileStream << "        ERREUR NECESSITANT ARRET DU PROGRAMME" << endl;
  fileStream << " - number du CPU : " << rankCpu << endl;
  fileStream << " - file :  " << m_fichier.c_str() << " ligne : " << m_ligne << endl;
  fileStream << " - message : " << m_message.c_str() << " " << m_value << endl;
  fileStream << "=> Corriger error et relancer le run" << endl;
  fileStream << "-------------------------------------------------------" << endl;

  fileStream.close();
}

//***********************************************************************

void Errors::arretCodeApresError(vector<Errors> &errors)
{
  //Affichage ecran des errors
  for (unsigned int e = 0; e < errors.size(); e++) {
    errors[e].afficheError();
  }
  //Arret propre
  //run.finalize();
  exit(0);
}

//***********************************************************************
