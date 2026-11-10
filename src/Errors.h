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
enum TypeError
{
  WARNING = 0,
  ERROR   = 1
};

// Macro allowing to use assert with error message
#ifndef NDEBUG
  #define assertM(condition, message)                                                          \
    do {                                                                                       \
      if (!(condition)) {                                                                      \
        std::cerr << "---------------------------------------------------------" << std::endl  \
                  << "Error assertion not verified" << std::endl                               \
                  << "  file: " << __FILE__ << std::endl                                       \
                  << "  line: " << __LINE__ << std::endl                                       \
                  << "  assertion: `" #condition "` failed" << std::endl                       \
                  << "  " << message << std::endl                                              \
                  << "---------------------------------------------------------" << std::endl; \
        std::exit(EXIT_FAILURE);                                                               \
      }                                                                                        \
    } while (false)
#else
  #define ASSERT(condition, message) \
    do {                             \
    } while (false)
#endif

//***************************************************************************************
//--------------------------Computing errors/warnings definition-------------------------
//***************************************************************************************

class Errors
{
  public:
    Errors();
    Errors(const std::string& message_, const char* sourceFile = "unknown", int lineNumber_ = -1);
    virtual ~Errors();

    static void errorMessage(const std::string& message_);
    static void errorMessage(const std::string& message_, double value);
    static void errorMessage(const std::string& message_, const char* sourceFile, int lineNumber_);

    void setError(const std::string& message_, const char* sourceFile = "", int lineNumber_ = -1);
    void setError(const std::string& message_, const double value);
    std::string errorToString(const int& num, const int& ErrorType, const int iteration = -1) const;
    void displayError(const int& num) const;

    //Accessor
    int getState();

    static constexpr int defaultInt       = 0;
    static constexpr int defaultIntNeg    = -1;
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
    ErrorECOGEN(std::string infoError_ = "", const char* sourceFile = "", int lineNumber_ = -1, int errorCode = 1) :
      std::exception(), m_errorCode(errorCode), m_lineNumber(lineNumber_), m_sourceFile(sourceFile), m_infoError(infoError_)
    {}
    virtual ~ErrorECOGEN() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception ECOGEN: fix and restart run"; }
    //***************
    std::string infoError(void) const throw()
    {
      std::stringstream message_;
      message_ << "--------------------------------------------------" << std::endl;
      message_ << this->what() << std::endl;
      message_ << "****************************************" << std::endl;
      if (m_sourceFile != "") {
        message_ << " infos on exception in code source :" << std::endl;
        message_ << "  file: '" << m_sourceFile << "'" << std::endl;
        if (m_lineNumber != -1) message_ << "  line: " << m_lineNumber << std::endl;
      }
      message_ << this->additionalInfo() << std::endl;
      message_ << "--------------------------------------------------" << std::endl;
      return message_.str();
    }
    //***************
    virtual std::string additionalInfo(void) const throw() { return m_infoError; }
    //***************
    int getErrorCode() { return m_errorCode; }
    //***************
  private:
    int m_errorCode;
    int m_lineNumber;
    std::string m_sourceFile;
    std::string m_infoError;
};

//---------------------------------------------------------------

//Exception handling on input files (XML, mesh)
//---------------------------------------------

class ErrorInput : public ErrorECOGEN
{
  public:
    //***************
    ErrorInput(std::string fileInput = "", const char* sourceFile = "", int lineNumber_ = -1, int errorCode = 2) :
      ErrorECOGEN(), m_errorCode(errorCode), m_lineNumber(lineNumber_), m_sourceFile(sourceFile), m_fileInput(fileInput)
    {}
    virtual ~ErrorInput() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception during reading input file"; }
    //***************
    std::string infoError(void) const throw()
    {
      std::stringstream message_;
      message_ << "--------------------------------------------------" << std::endl;
      message_ << this->what() << std::endl;
      message_ << "****************************************" << std::endl;
      if (m_fileInput != "") {
        message_ << " Input file concerned: '" << m_fileInput << "'" << std::endl;
      }
      if (m_sourceFile != "") {
        message_ << " Infos on exception in source code:" << std::endl;
        message_ << "  file: '" << m_sourceFile << "'" << std::endl;
        if (m_lineNumber != -1) message_ << "  line: " << m_lineNumber << std::endl;
      }
      message_ << this->additionalInfo();
      message_ << "--------------------------------------------------" << std::endl;
      return message_.str();
    }
    //***************
    virtual std::string additionalInfo(void) const throw() { return ""; }
    //***************
    int lineNumber() const throw() { return m_lineNumber; }
    //***************
    std::string fileName() const throw() { return m_sourceFile; }
    //***************
    int getErrorCode() { return m_errorCode; }
    //***************
  protected:
    int m_errorCode;
    int m_lineNumber;
    std::string m_sourceFile;
    std::string m_fileInput;
};

//---------------------------------------------------------------

//Exception handling on unstructured mesh file
//--------------------------------------------
class ErrorMeshNS : public ErrorInput
{
  public:
    //***************
    ErrorMeshNS(std::string fileMesh = "", const char* sourceFile = "", int lineNumber_ = -1, int errorCode = 2) :
      ErrorInput(fileMesh, sourceFile, lineNumber_, errorCode)
    {}
    virtual ~ErrorMeshNS() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception during reading unstructured mesh file: file not found or incorrect structure"; }
    //***************
    std::string infoError(void) const throw()
    {
      std::stringstream message_;
      message_ << "--------------------------------------------------" << std::endl;
      message_ << this->what() << std::endl;
      message_ << "****************************************" << std::endl;
      if (m_fileInput != "") {
        message_ << " mesh file concerned: '" << m_fileInput << "'" << std::endl;
      }
      if (m_sourceFile != "") {
        message_ << " infos on exception in code source :" << std::endl;
        message_ << "  file: '" << m_sourceFile << "'" << std::endl;
        if (m_lineNumber != -1) message_ << "  line: " << m_lineNumber << std::endl;
      }
      message_ << this->additionalInfo();
      message_ << "--------------------------------------------------" << std::endl;
      return message_.str();
    }
};

//---------------------------------------------------------------

//Exception handling on input XML files
//-------------------------------------

class ErrorXML : public ErrorInput
{
  public:
    //***************
    ErrorXML(std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1, int errorCode = 2) :
      ErrorInput(fileXML, sourceFile, lineNumber_, errorCode)
    {}
    virtual ~ErrorXML() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception during reading XML file: file not found or incorrect structure"; }
    //***************
    std::string infoError(void) const throw()
    {
      std::stringstream message_;
      message_ << "--------------------------------------------------" << std::endl;
      message_ << this->what() << std::endl;
      message_ << "****************************************" << std::endl;
      if (m_fileInput != "") {
        message_ << " XML file concerned: '" << m_fileInput << "'" << std::endl;
      }
      if (m_sourceFile != "") {
        message_ << " infos on exception in code source :" << std::endl;
        message_ << "  file: '" << m_sourceFile << "'" << std::endl;
        if (m_lineNumber != -1) message_ << "  line: " << m_lineNumber << std::endl;
      }
      message_ << this->additionalInfo();
      message_ << "--------------------------------------------------" << std::endl;
      return message_.str();
    }
};

//---------------------------------------------------------------

class ErrorXMLMessage : public ErrorXML
{
  public:
    //***************
    ErrorXMLMessage(std::string message_ = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXML(fileXML, sourceFile, lineNumber_), m_message(message_)
    {}
    virtual ~ErrorXMLMessage() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: input files are not compatible"; }
    //***************
    virtual std::string additionalInfo(void) const throw()
    {
      std::stringstream message_;
      message_ << " details: " << m_message << std::endl;
      return message_.str();
    }
    //***************
    std::string message() const throw() { return m_message; }
    //***************
  private:
    std::string m_message;
};

//---------------------------------------------------------------

class ErrorXMLRacine : public ErrorXML
{
  public:
    //***************
    ErrorXMLRacine(std::string root_ = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXML(fileXML, sourceFile, lineNumber_), m_root(root_)
    {}
    virtual ~ErrorXMLRacine() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: root not found"; }
    //***************
    virtual std::string additionalInfo(void) const throw()
    {
      std::stringstream message_;
      message_ << " root name searched: '" << m_root << "'" << std::endl;
      return message_.str();
    }
    //***************
    std::string root() const throw() { return m_root; }
    //***************
  private:
    std::string m_root;
};

//---------------------------------------------------------------

class ErrorXMLElement : public ErrorXML
{
  public:
    //***************
    ErrorXMLElement(std::string element_ = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXML(fileXML, sourceFile, lineNumber_), m_element(element_)
    {}
    virtual ~ErrorXMLElement() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: element not found"; }
    //***************
    virtual std::string additionalInfo(void) const throw()
    {
      std::stringstream message_;
      message_ << " element name searched: '" << m_element << "'" << std::endl;
      return message_.str();
    }
    //***************
    std::string element() const throw() { return m_element; }
    //***************
  private:
    std::string m_element;
};

//---------------------------------------------------------------

class ErrorXMLAttribut : public ErrorXML
{
  public:
    //***************
    ErrorXMLAttribut(std::string attribute_ = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXML(fileXML, sourceFile, lineNumber_), m_attribute(attribute_)
    {}
    virtual ~ErrorXMLAttribut() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: attribute not found"; }
    //***************
    virtual std::string additionalInfo(void) const throw()
    {
      std::stringstream message_;
      message_ << " attribute name searched: '" << m_attribute << "'" << std::endl;
      return message_.str();
    }
    //***************
    std::string attribute() const throw() { return m_attribute; }
    //***************
  private:
    std::string m_attribute;
};

//---------------------------------------------------------------

class ErrorXMLDev : public ErrorXML
{
  public:
    //***************
    ErrorXMLDev(std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) : ErrorXML(fileXML, sourceFile, lineNumber_) {}
    virtual ~ErrorXMLDev() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: this portion of code is under development"; }
    //***************
  private:
};

//---------------------------------------------------------------

class ErrorXMLLimite : public ErrorXML
{
  public:
    //***************
    ErrorXMLLimite(std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) : ErrorXML(fileXML, sourceFile, lineNumber_) {}
    virtual ~ErrorXMLLimite() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: error in boundary conditions"; }
    //***************
  private:
};

//---------------------------------------------------------------

class ErrorXMLTermeSource : public ErrorXML
{
  public:
    //***************
    ErrorXMLTermeSource(std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) : ErrorXML(fileXML, sourceFile, lineNumber_) {}
    virtual ~ErrorXMLTermeSource() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: error in source term selection"; }
    //***************
  private:
};

//---------------------------------------------------------------

class ErrorXMLEOS : public ErrorXML
{
  public:
    //***************
    ErrorXMLEOS(std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) : ErrorXML(fileXML, sourceFile, lineNumber_) {}
    virtual ~ErrorXMLEOS() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: error on equation of state"; }
    //***************
  private:
};

//---------------------------------------------------------------

class ErrorXMLEOSUnknown : public ErrorXML
{
  public:
    //***************
    ErrorXMLEOSUnknown(std::string typeEOS = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXML(fileXML, sourceFile, lineNumber_), m_typeEOS(typeEOS)
    {}
    virtual ~ErrorXMLEOSUnknown() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: error on equation of state"; }
    //***************
    virtual std::string additionalInfo(void) const throw()
    {
      std::stringstream message_;
      message_ << " EOS type: '" << m_typeEOS << "' unknown" << std::endl;
      return message_.str();
    }
    //***************
  private:
    std::string m_typeEOS;
};

//---------------------------------------------------------------

class ErrorXMLDomaineUnknown : public ErrorXML
{
  public:
    //***************
    ErrorXMLDomaineUnknown(std::string typeDomain = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXML(fileXML, sourceFile, lineNumber_), m_typeDomain(typeDomain)
    {}
    virtual ~ErrorXMLDomaineUnknown() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: error on domain CI"; }
    //***************
    virtual std::string additionalInfo(void) const throw()
    {
      std::stringstream message_;
      message_ << " domain type: '" << m_typeDomain << "' unknown" << std::endl;
      return message_.str();
    }
    //***************
  private:
    std::string m_typeDomain;
};

//---------------------------------------------------------------

class ErrorXMLBoundCondUnknown : public ErrorXML
{
  public:
    //***************
    ErrorXMLBoundCondUnknown(std::string typeBoundCond = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXML(fileXML, sourceFile, lineNumber_), m_typeBoundCond(typeBoundCond)
    {}
    virtual ~ErrorXMLBoundCondUnknown() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: error in boundary condition"; }
    //***************
    virtual std::string additionalInfo(void) const throw()
    {
      std::stringstream message_;
      message_ << " boundary type: '" << m_typeBoundCond << "' unknown" << std::endl;
      return message_.str();
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
    ErrorXMLEtat(std::string nameEtat = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXML(fileXML, sourceFile, lineNumber_), m_nameEtat(nameEtat)
    {}
    virtual ~ErrorXMLEtat() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: error in state"; }
    //***************
    virtual std::string additionalInfo(void) const throw()
    {
      std::stringstream message_;
      message_ << " state: '" << m_nameEtat << "' not found or incomplete" << std::endl;
      return message_.str();
    }
    //***************
  private:
    std::string m_nameEtat;
};

//---------------------------------------------------------------

class ErrorXMLMaterialUnknown : public ErrorXML
{
  public:
    //***************
    ErrorXMLMaterialUnknown(std::string nameMateriau = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXML(fileXML, sourceFile, lineNumber_), m_nameMateriau(nameMateriau)
    {}
    virtual ~ErrorXMLMaterialUnknown() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: error in state CI"; }
    //***************
    virtual std::string additionalInfo(void) const throw()
    {
      std::stringstream message_;
      message_ << " material type: '" << m_nameMateriau << "' unknown" << std::endl;
      return message_.str();
    }
    //***************
  private:
    std::string m_nameMateriau;
};

//---------------------------------------------------------------

class ErrorXMLStretching : public ErrorXMLMessage
{
  public:
    //***************
    ErrorXMLStretching(std::string message_, std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXMLMessage(message_, fileXML, sourceFile, lineNumber_)
    {}
    virtual ~ErrorXMLStretching() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: error on stretching definition"; }
};

//---------------------------------------------------------------

class ErrorXMLRelaxation : public ErrorXML
{
  public:
    //***************
    ErrorXMLRelaxation(std::string nameRelaxation = "", std::string fileXML = "", const char* sourceFile = "", int lineNumber_ = -1) :
      ErrorXML(fileXML, sourceFile, lineNumber_), m_nameRelaxation(nameRelaxation)
    {}
    virtual ~ErrorXMLRelaxation() throw() {}
    //***************
    virtual const char* what(void) const throw() { return "Exception on XML file: error on relaxation"; }
    //***************
    virtual std::string additionalInfo(void) const throw()
    {
      std::stringstream message_;
      message_ << " relaxation type: '" << m_nameRelaxation << "' unknown or already added" << std::endl;
      return message_.str();
    }
    //***************
  private:
    std::string m_nameRelaxation;
};

//---------------------------------------------------------------

#endif // ERRORS_H
