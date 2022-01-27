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

#ifndef ERRORS_H
#define ERRORS_H

//Definition of Error classes of all kind
//Exceptions on input XML files

#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <list>
#include <string>
#include <algorithm>

//! \brief     Enumeration for the type of error (warning, error)
enum TypeError { WARNING = 0, ERROR = 1 };

// Macro allowing to use assert with error message
# ifndef NDEBUG
#   define assertM(condition, message) \
  do {\
    if (!(condition)) {  \
    std::cerr << "---------------------------------------------------------" \
      << std::endl << "Error assertion not verified" << std::endl << \
      "  file: " << __FILE__ << std::endl << \
      "  line: " << __LINE__ << std::endl << \
      "  assertion: `" #condition "` failed" << std::endl \
      << "  " << message << std::endl \
      << "---------------------------------------------------------" << std::endl; \
      std::exit(EXIT_FAILURE); \
    } \
  } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

//***************************************************************************************
//--------------------------Computing errors/warnings definition-------------------------
//***************************************************************************************

class Errors
{
public:
  Errors();
  Errors(const std::string& message, const char* sourceFile = "unknown", int lineNumber = -1);
  virtual ~Errors();

  static void errorMessage(const std::string& message);
  static void errorMessage(const std::string& message, double value);
  static void prepareErrorFiles(const std::string& folder);

  void setError(const std::string& message, const char* sourceFile = "", int lineNumber = -1);
  void setError(const std::string& message, const double value);
  void displayError(const int& num);
  void writeErrorInFile(const int& num, const std::string& folder, const int& ErrorType);

  //Accessor
  int getState();
  
  static constexpr int defaultInt = 0;
  static constexpr int defaultIntNeg = -1;
  static constexpr double defaultDouble = 0.;
  static const std::string defaultString;

private:
  std::string m_message;
  int m_state;
  std::string m_file;
  int m_line;
  double m_value; //!< Allows you to send an additionnal piece of information
}; 

extern std::vector<Errors> errors;
extern std::vector<Errors> warnings; //FP//TODO// To improve by a specific class

//***************************************************************************************
//--------------------------------------EXCEPTIONS---------------------------------------
//***************************************************************************************

//Exception handling on error code ECOGEN
//---------------------------------------
class ErrorECOGEN : public std::exception
{
public:
  //***************
  ErrorECOGEN(std::string infoError = "", const char* sourceFile = "", int lineNumber = -1, int errorCode = 1):
    std::exception(), m_errorCode(errorCode), m_lineNumber(lineNumber), m_sourceFile(sourceFile), m_infoError(infoError) {}
  virtual ~ErrorECOGEN() throw() {}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception ECOGEN: fix and restart run";
  }
  //***************
  std::string infoError(void) const throw()
  {
    std::stringstream message;
    message << "--------------------------------------------------" << std::endl;
    message << this->what() << std::endl;
    message << "****************************************" << std::endl;
    if (m_sourceFile != "")
    {
      message << " infos on exception in code source :" << std::endl;
      message << "  file: '" << m_sourceFile << "'" << std::endl;
      if (m_lineNumber != -1) message << "  line: " << m_lineNumber << std::endl;
    }
    message << this->additionalInfo() << std::endl;
    message << "--------------------------------------------------" << std::endl;
    return message.str();
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    return m_infoError;
  }
  //***************
  int getErrorCode() { return m_errorCode; }
  //***************
private:
  int m_errorCode;
  int m_lineNumber;
  std::string m_sourceFile;
  std::string m_infoError;
};


//Exception handling on input XML files
//-------------------------------------

class ErrorXML : public ErrorECOGEN
{
public:
  //***************
  ErrorXML(std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1, int errorCode = 2) :
    ErrorECOGEN(), m_errorCode(errorCode), m_lineNumber(lineNumber), m_sourceFile(sourceFile), m_fileXML(fileXML){}
  virtual ~ErrorXML() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception during reading XML file: file not found or incorrect structure";
  }
  //***************
  std::string infoError(void) const throw()
  {
    std::stringstream message;
    message << "--------------------------------------------------" << std::endl;
    message << this->what() << std::endl;
    message << "****************************************" << std::endl;
    if (m_fileXML != "") { message << " XML file concerned: '" << m_fileXML << "'" << std::endl; }
    if (m_sourceFile != "")
    {
      message << " infos on exception in code source :" << std::endl;
      message << "  file: '" << m_sourceFile << "'" << std::endl;
      if (m_lineNumber != -1) message << "  line: " << m_lineNumber << std::endl;
    }
    message << this->additionalInfo();
    message << "--------------------------------------------------" << std::endl;
    return message.str();
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    return "";
  }
  //***************
  int lineNumber() const throw () { return m_lineNumber; }
  //***************
  std::string fileName() const throw () { return m_sourceFile; }
  //***************
  int getErrorCode() { return m_errorCode; }
  //***************
private:
  int m_errorCode;
  int m_lineNumber;
  std::string m_sourceFile;
  std::string m_fileXML;
};

//---------------------------------------------------------------

class ErrorXMLMessage : public ErrorXML
{
public:
  //***************
  ErrorXMLMessage(std::string message = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber), m_message(message){}
  virtual ~ErrorXMLMessage() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: input files are not compatible";
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    std::stringstream message;
    message << " details: " << m_message << std::endl;
    return message.str();
  }
  //***************
  std::string message() const throw () { return m_message; }
  //***************
private:
  std::string m_message;
};

//---------------------------------------------------------------

class ErrorXMLRacine : public ErrorXML
{
public:
  //***************
  ErrorXMLRacine(std::string root = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber), m_root(root){}
  virtual ~ErrorXMLRacine() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: root not found";
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    std::stringstream message;
    message << " root name searched: '" << m_root << "'" << std::endl;
    return message.str();
  }
  //***************
  std::string root() const throw () { return m_root; }
  //***************
private:
  std::string m_root;
};

//---------------------------------------------------------------

class ErrorXMLElement : public ErrorXML
{
public:
  //***************
  ErrorXMLElement(std::string element = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber), m_element(element){}
  virtual ~ErrorXMLElement() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: element not found";
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    std::stringstream message;
    message << " element name searched: '" << m_element << "'" << std::endl;
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
  ErrorXMLAttribut(std::string attribute = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber), m_attribute(attribute){}
  virtual ~ErrorXMLAttribut() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: attribute not found";
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    std::stringstream message;
    message << " attribute name searched: '" << m_attribute << "'" << std::endl;
    return message.str();
  }
  //***************
  std::string attribute() const throw () { return m_attribute; }
  //***************
private:
  std::string m_attribute;
};

//---------------------------------------------------------------

class ErrorXMLDev : public ErrorXML
{
public:
  //***************
  ErrorXMLDev(std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber){}
  virtual ~ErrorXMLDev() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: this portion of code is under development";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErrorXMLLimite : public ErrorXML
{
public:
  //***************
  ErrorXMLLimite(std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber){}
  virtual ~ErrorXMLLimite() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: error in boundary conditions";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErrorXMLTermeSource : public ErrorXML
{
public:
  //***************
  ErrorXMLTermeSource(std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber){}
  virtual ~ErrorXMLTermeSource() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: error in source term selection";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErrorXMLEOS : public ErrorXML
{
public:
  //***************
  ErrorXMLEOS(std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber){}
  virtual ~ErrorXMLEOS() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: error on equation of state";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErrorXMLEOSInconnue : public ErrorXML
{
public:
  //***************
  ErrorXMLEOSInconnue(std::string typeEOS = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber), m_typeEOS(typeEOS){}
  virtual ~ErrorXMLEOSInconnue() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: error on equation of state";
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    std::stringstream message;
    message << " EOS type: '" << m_typeEOS << "' unknown" << std::endl;
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
  ErrorXMLDomaineInconnu(std::string typeDomain = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber), m_typeDomain(typeDomain){}
  virtual ~ErrorXMLDomaineInconnu() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: error on domain CI";
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    std::stringstream message;
    message << " domain type: '" << m_typeDomain << "' unknown" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_typeDomain;
};

//---------------------------------------------------------------

class ErrorXMLBoundCondInconnue : public ErrorXML
{
public:
  //***************
  ErrorXMLBoundCondInconnue(std::string typeBoundCond = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber), m_typeBoundCond(typeBoundCond) {}
  virtual ~ErrorXMLBoundCondInconnue() throw() {}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: error in boundary condition CL";
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    std::stringstream message;
    message << " boundary type: '" << m_typeBoundCond << "' unknown" << std::endl;
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
  ErrorXMLEtat(std::string nameEtat = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber), m_nameEtat(nameEtat){}
  virtual ~ErrorXMLEtat() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: error in state";
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    std::stringstream message;
    message << " state: '" << m_nameEtat << "' not found or incomplete" << std::endl;
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
  ErrorXMLMateriauInconnu(std::string nameMateriau = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber), m_nameMateriau(nameMateriau){}
  virtual ~ErrorXMLMateriauInconnu() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: error in state CI";
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    std::stringstream message;
    message << " material type: '" << m_nameMateriau << "' unknown" << std::endl;
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
  ErrorXMLStretching(std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber) {}
  virtual ~ErrorXMLStretching() throw() {}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: error on stretching definition";
  }
  //***************
private:
  std::string m_nameMateriau;
};

//---------------------------------------------------------------

class ErrorXMLRelaxation : public ErrorXML
{
public:
  //***************
  ErrorXMLRelaxation(std::string nameRelaxation = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber = -1) :
    ErrorXML(fileXML, sourceFile, lineNumber), m_nameRelaxation(nameRelaxation){}
  virtual ~ErrorXMLRelaxation() throw(){}
  //***************
  virtual const char* what(void) const throw()
  {
    return "Exception on XML file: error on relaxation";
  }
  //***************
  virtual std::string additionalInfo(void) const throw()
  {
    std::stringstream message;
    message << " relaxation type: '" << m_nameRelaxation << "' unknown or already added" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_nameRelaxation;
};

//---------------------------------------------------------------

#endif // ERRORS_H