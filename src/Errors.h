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

#ifndef ERRORS_H
#define ERRORS_H

//! \file      Errors.h
//! \author    F. Petitpas, S. Le Martelot, K. Schmidmayer, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

//Definitions de classes d'Errors tout type
//Error de base
//Exceptions sur lecture XML

#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>

// Macro permettant de faire un assert avec un message d'error.
# ifndef NDEBUG
#   define assertM(condition, message) \
  do {\
    if (!(condition)) {  \
    std::cerr << "---------------------------------------------------------" \
      << std::endl << "Error assertion non verify" << std::endl << \
      "  file : " << __FILE__ << std::endl << \
      "  ligne : " << __LINE__ << std::endl << \
      "  assertion : `" #condition "` a echouee" << std::endl \
      << "  " << message << std::endl \
      << "---------------------------------------------------------" << std::endl; \
      std::exit(EXIT_FAILURE); \
    } \
  } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

class Errors
{
public:
  Errors();
  Errors(const std::string &message, const char* fichierSource = "non renseigne", int numberLigne = -1);
  virtual ~Errors();

  static void errorMessage(const std::string &message);
  static void errorMessage(const std::string &message, double value);

  void setError(const std::string &message, const char* fichierSource = "", int numberLigne = -1);
  void setError(const std::string &message, const double value);
  void afficheError();
  void ecritErrorFichier();
  static void arretCodeApresError(std::vector<Errors> &errors);

  //Accesseur
  int getEtat();

  static constexpr int defaultInt = 0;
  static constexpr int defaultIntNeg = -1;
  static constexpr double defaultDouble = 0.;
  static const std::string defaultString;

private:
  std::string m_message;
  int m_state;
  std::string m_fichier;
  int m_ligne;
  double m_value; //!< permet de faire remonter une information en plus
};

extern std::vector<Errors> errors;

//Gestion des exceptions sur error code ECOGEN
//---------------------------------------------
class ErrorECOGEN : public std::exception
{
public:
  //***************
  ErrorECOGEN(std::string infoError = "", const char* fichierSource = "", int numberLigne = -1):
    std::exception(), m_infoError(infoError), m_numberLigne(numberLigne), m_fichierSource(fichierSource) {}
  virtual ~ErrorECOGEN() throw() {}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception ECOGEN : corriger et relancer le run";
  }
  //***************
  std::string infoError(void) const throw()
  {
    std::stringstream message;
    message << "--------------------------------------------------" << std::endl;
    message << this->what() << std::endl;
    message << "****************************************" << std::endl;
    if (m_fichierSource != "")
    {
      message << " infos sur exception code source :" << std::endl;
      message << "  file : '" << m_fichierSource << "'" << std::endl;
      if (m_numberLigne != -1) message << "  ligne : " << m_numberLigne << std::endl;
    }
    message << this->infosAdditionelles() << std::endl;
    message << "--------------------------------------------------" << std::endl;
    return message.str();
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    return m_infoError;
  }
  //***************

private:
  int m_numberLigne;
  std::string m_fichierSource;
  std::string m_infoError;
};


//Gestion des exceptions sur fichiers entrees XML
//---------------------------------------------------------------

class ErrorXML : public ErrorECOGEN
{
public:
  //***************
  ErrorXML(std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorECOGEN(), m_fichierXML(fichierXML), m_fichierSource(fichierSource), m_numberLigne(numberLigne){}
  virtual ~ErrorXML() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur lecture file XML : file introuvable ou structure incorrecte";
  }
  //***************
  std::string infoError(void) const throw()
  {
    std::stringstream message;
    message << "--------------------------------------------------" << std::endl;
    message << this->what() << std::endl;
    message << "****************************************" << std::endl;
    if (m_fichierXML != "") { message << " file XML concerne : '" << m_fichierXML << "'" << std::endl; }
    if (m_fichierSource != "")
    {
      message << " infos sur exception code source :" << std::endl;
      message << "  file : '" << m_fichierSource << "'" << std::endl;
      if (m_numberLigne != -1) message << "  ligne : " << m_numberLigne << std::endl;
    }
    message << this->infosAdditionelles();
    message << "--------------------------------------------------" << std::endl;
    return message.str();
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    return "";
  }
  //***************
  int numberLigne() const throw () { return m_numberLigne; }
  //***************
  std::string fileName() const throw () { return m_fichierSource; }
  //***************
private:
  int m_numberLigne;
  std::string m_fichierSource;
  std::string m_fichierXML;
};

//---------------------------------------------------------------

class ErrorXMLRacine : public ErrorXML
{
public:
  //***************
  ErrorXMLRacine(std::string racine = "", std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne), m_racine(racine){}
  virtual ~ErrorXMLRacine() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : racine introuvable";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " name de la racine recherchee : '" << m_racine << "'" << std::endl;
    return message.str();
  }
  //***************
  std::string racine() const throw () { return m_racine; }
  //***************
private:
  std::string m_racine;
};

//---------------------------------------------------------------

class ErrorXMLElement : public ErrorXML
{
public:
  //***************
  ErrorXMLElement(std::string element = "", std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne), m_element(element){}
  virtual ~ErrorXMLElement() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : element introuvable ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " name de l element recherche : '" << m_element << "'" << std::endl;
    return message.str();
  }
  //***************
  std::string element() const throw () { return m_element; }
  //***************
private:
  std::string m_element;
};

//---------------------------------------------------------------

class ErrorXMLAttribut : public ErrorXML
{
public:
  //***************
  ErrorXMLAttribut(std::string attribut = "", std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne), m_attribut(attribut){}
  virtual ~ErrorXMLAttribut() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : attribut introuvable ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " name de l attribut recherche : '" << m_attribut << "'" << std::endl;
    return message.str();
  }
  //***************
  std::string attribut() const throw () { return m_attribut; }
  //***************
private:
  std::string m_attribut;
};

//---------------------------------------------------------------

class ErrorXMLDev : public ErrorXML
{
public:
  //***************
  ErrorXMLDev(std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne){}
  virtual ~ErrorXMLDev() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : morceaux de code en cours de developpement ";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErrorXMLLimite : public ErrorXML
{
public:
  //***************
  ErrorXMLLimite(std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne){}
  virtual ~ErrorXMLLimite() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : error sur conditions aux limites ";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErrorXMLTermeSource : public ErrorXML
{
public:
  //***************
  ErrorXMLTermeSource(std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne){}
  virtual ~ErrorXMLTermeSource() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : error sur choix des termes sources ";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErrorXMLEOS : public ErrorXML
{
public:
  //***************
  ErrorXMLEOS(std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne){}
  virtual ~ErrorXMLEOS() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : error sur equation state ";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErrorXMLEOSInconnue : public ErrorXML
{
public:
  //***************
  ErrorXMLEOSInconnue(std::string typeEOS = "", std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne), m_typeEOS(typeEOS){}
  virtual ~ErrorXMLEOSInconnue() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : error sur equation state ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " type EOS : '" << m_typeEOS << "' inconnu" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_typeEOS;
};

//---------------------------------------------------------------

class ErrorXMLDomaineInconnu : public ErrorXML
{
public:
  //***************
  ErrorXMLDomaineInconnu(std::string typeDomaine = "", std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne), m_typeDomaine(typeDomaine){}
  virtual ~ErrorXMLDomaineInconnu() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : error sur domain CI ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " type Domaine : '" << m_typeDomaine << "' inconnu" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_typeDomaine;
};

//---------------------------------------------------------------

class ErrorXMLBoundCondInconnue : public ErrorXML
{
public:
  //***************
  ErrorXMLBoundCondInconnue(std::string typeBoundCond = "", std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne), m_typeBoundCond(typeBoundCond) {}
  virtual ~ErrorXMLBoundCondInconnue() throw() {}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : error sur Condition Limite CL ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " type Limite : '" << m_typeBoundCond << "' inconnu" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_typeBoundCond;
};

//---------------------------------------------------------------

class ErrorXMLEtat : public ErrorXML
{
public:
  //***************
  ErrorXMLEtat(std::string nameEtat = "", std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne), m_nameEtat(nameEtat){}
  virtual ~ErrorXMLEtat() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : error sur state ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " Etat : '" << m_nameEtat << "' non trouve ou incomplet" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_nameEtat;
};

//---------------------------------------------------------------

class ErrorXMLMateriauInconnu : public ErrorXML
{
public:
  //***************
  ErrorXMLMateriauInconnu(std::string nameMateriau = "", std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne), m_nameMateriau(nameMateriau){}
  virtual ~ErrorXMLMateriauInconnu() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : error sur state CI ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " type Materiau : '" << m_nameMateriau << "' inconnu" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_nameMateriau;
};

//---------------------------------------------------------------

class ErrorXMLStretching : public ErrorXML
{
public:
  //***************
  ErrorXMLStretching(std::string fichierXML = "", const char* fichierSource = "", int numberLigne = -1) :
    ErrorXML(fichierXML, fichierSource, numberLigne) {}
  virtual ~ErrorXMLStretching() throw() {}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur file XML : error on stretching definition ";
  }
  //***************
private:
  std::string m_nameMateriau;
};

//---------------------------------------------------------------

#endif // ERRORS_H 