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

#include "Errors.h"
#include "Run.h"

std::vector<Errors> errors;
std::vector<Errors> warnings;

//***********************************************************************

Errors::Errors() : m_message("not specified"), m_state(0), m_line(0), m_value(0.)
{}

//***********************************************************************

Errors::Errors(const std::string& message, const char* sourceFile, int lineNumber) :
  m_message(message), m_file(sourceFile), m_line(lineNumber), m_value(0.)
{
  m_state = 1;
}

//***********************************************************************

Errors::~Errors(){}

//***********************************************************************

void Errors::errorMessage(const std::string& message)
{
  std::stringstream numCPU;
  numCPU << rankCpu;

  std::cout << std::endl << "-------------------------------------------------------" << std::endl;
  std::cout << "ERROR on CPU " + numCPU.str() + ": " << message.c_str() << std::endl; 
  std::cout << "=> Fix error and restart run" << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
}

//***********************************************************************

void Errors::errorMessage(const std::string& message, double value)
{
  std::stringstream numCPU;
  numCPU << rankCpu;

  std::cout << std::endl << "-------------------------------------------------------" << std::endl;
  std::cout << "ERROR on CPU " + numCPU.str() + " : " << message.c_str() << " " << value << std::endl;
  std::cout << "=> Fix error and restart run" << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
}

//***********************************************************************

void Errors::prepareErrorFiles(const std::string& folder)
{
  std::ofstream fileStream;
  std::stringstream myStream;
  std::string fileName;

  fileName = "error";
  myStream << folder << "errorsAndWarnings/" << fileName << "_CPU" << rankCpu << ".out";
  fileStream.open((myStream.str()).c_str());
  fileStream.close();

  myStream.str("");
  myStream.clear();
  fileName = "warning";
  myStream << folder << "errorsAndWarnings/" << fileName << "_CPU" << rankCpu << ".out";
  fileStream.open((myStream.str()).c_str());
  fileStream.close();
}

//***********************************************************************

void Errors::setError(const std::string& message, const char* sourceFile, int lineNumber)
{
  m_message = message;
  m_state = 1;
  m_file = sourceFile;
  m_line = lineNumber;
}

//***********************************************************************

void Errors::setError(const std::string& message, const double value)
{
  m_message = message;
  m_state = 1;
  m_value = value;
}

//***********************************************************************

int Errors::getState()
{
  return m_state;
}

//***********************************************************************

void Errors::displayError(const int& num)
{
  std::cout << std::endl << "-------------------------------------------------------" << std::endl;
  std::cout << "        ERROR " << num << " REQUIRING PROGRAM SHUTDOWN" << std::endl;
  std::cout << " - CPU number: " << rankCpu << std::endl;
  std::cout << " - file: " << m_file.c_str() << " line : " << m_line << std::endl;
  std::cout << " - message: " << m_message.c_str() << std::endl;
  std::cout << "=> Fix error and restart run" << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
}

//***********************************************************************

void Errors::writeErrorInFile(const int& num, const std::string& folder, const int& ErrorType)
{
  std::ofstream fileStream;
  std::stringstream myStream;

  std::string fileName;
  if      (ErrorType == ERROR)   { fileName = "error"; }
  else if (ErrorType == WARNING) { fileName = "warning"; }
  myStream << folder <<"errorsAndWarnings/" <<  fileName << "_CPU" << rankCpu << ".out";
  fileStream.open((myStream.str()).c_str(), std::ofstream::app);
  fileStream << "-------------------------------------------------------" << std::endl;
  if      (ErrorType == ERROR)   { fileStream << "        ERROR " << num << " REQUIRING PROGRAM SHUTDOWN" << std::endl; }
  else if (ErrorType == WARNING) { fileStream << "        WARNING " << num << " REQUIRING ATTENTION" << std::endl; }
  fileStream << " - CPU number: " << rankCpu << std::endl;
  fileStream << " - file: " << m_file.c_str() << " line: " << m_line << std::endl;
  fileStream << " - message: " << m_message.c_str() << std::endl;
  if      (ErrorType == ERROR)   { fileStream << "=> Fix error and restart run" << std::endl; }
  fileStream << "-------------------------------------------------------" << std::endl;
  fileStream.close();
}

//***********************************************************************

constexpr int Errors::defaultInt;
constexpr int Errors::defaultIntNeg;
constexpr double Errors::defaultDouble;
const std::string Errors::defaultString = "NA";

//***********************************************************************
