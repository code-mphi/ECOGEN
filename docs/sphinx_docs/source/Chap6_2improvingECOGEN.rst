********************
Contribute to ECOGEN
********************

Before contributing, please read the section :ref:`Sec:Licence` informations.

API documentation
=================
An API_ documentation is available. It is generated thanks to Doxygen_ and special comments in the source code. Interested developer will find informations on the different classes and functions used in ECOGEN.

When contributing to the project and to automatise code documentation using Doxygen, comments should be insert as follow in header files.

At the head of files
--------------------

.. code-block:: c++

	//! \file      file name
	//! \author    authors names
	//! \version   1.0
	//! \date      12 Novembre 2009
	//! \brief     brief description

Before class definition
-----------------------

.. code-block:: c++

	//! \class     class name
	//! \brief     brief descritption
	//! \details    Detailed description
	//!                  continue description

Before function and method prototypes
-------------------------------------

.. code-block:: c++

	//! \brief      brief description
	//! \details    Detailed description
	//!                  continue description
	//! \param      parameter name         description
	//! \param      parameter name         description
	//! \return     return description

Class member descrption
-----------------------

.. code-block:: c++

	type m_variable; //!< member description


Coding constraints
==================

Developers are thanks to respect some constraints when contributing to the project.

Variable and Classe names
-------------------------

  1. Variable name should began with lowercase letter and be self-understandable. Each new word began with a uppercase letter.
  
  .. code-block:: c++

  	int myInteger; vector<double *> vectorOfDoublePointer; etc.

  2. Class attribute should began by "m\_". 
    
  .. code-block:: c++

		int m_myInteger; double * m_doublePointer; etc.

  3. Class name should began with an uppercase Letter.
    
  .. code-block:: c++

		class MyClass;

Developer personnal comments - flags
-------------------------------------

Developer personnal comments should be included using the following template:	

//Developer//KeyWord// comments
  
.. code-block:: c++

	//FP//DEV// comment, description

Here is the list of keyword to use :

  
.. code-block:: c++

	//DEV//    in developement
	//Q//      question to dig
	//TODO//   should be done in the future
	//ERR//    error : to correct ASAP
	//ID//     idea
	//ICI//    Stop developement position
	//VERIF//  to verify : is it needed ?
	//TEST//   test : To delete ASAP

Git-hub submit
==============

For each modification, a comment should be prepared to be included to the commit message for Git.

.. _API: https://code-mphi.github.io/ECOGEN/docs/doxygen_docs/index.html
.. _Doxygen: http://www.doxygen.nl/