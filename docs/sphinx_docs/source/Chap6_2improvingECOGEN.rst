********************
Contribute to ECOGEN
********************

Before contributing, please read the section :ref:`Sec:Licence` information.

API documentation
=================
An API_ documentation is available. It is generated thanks to Doxygen_ and special comments in the source code. Interested developers will find information on the different classes, methods and functions used in ECOGEN.

When contributing to the project and to automatize the code documentation using Doxygen, comments should be inserted as follow in header files.

At the head of files
--------------------

.. code-block:: c++

	//! \file      file name
	//! \author    authors names
	//! \version   1.1
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

Class member description
------------------------

.. code-block:: c++

	type m_variable; //!< member description


Coding constraints
==================

Developers are thanks to respect some constraints when contributing to the project.

Variable and class names
------------------------

  1. Variable names should began with lowercase letter and be self-understandable. Each new word began with a uppercase letter.
  
  .. code-block:: c++

  	int myInteger; vector<double *> vectorOfDoublePointer; etc.

  2. Class attributes should began by "m\_". 
    
  .. code-block:: c++

		int m_myInteger; double * m_doublePointer; etc.

  3. Class names should began with an uppercase Letter.
    
  .. code-block:: c++

		class MyClass;

Developer personnal comments - flags
------------------------------------

Developer personnal comments should be included using the following template:	

//DeveloperInitials//KeyWord// comments
  
.. code-block:: c++

	//FP//DEV// comment, description

Here is the list of keyword to use :

  
.. code-block:: c++

	//DEV//    in development
	//Q//      question to dig
	//TODO//   should be done in the future
	//ERR//    error: To correct ASAP
	//ID//     idea
	//ICI//    stop development position
	//VERIF//  to verify: Is it needed?
	//TEST//   test: To delete ASAP

GitHub submit
=============

For each modification, a comment should be prepared to be included to the commit message for Git.
Issues and contributions to ECOGEN project are possible directly on GitHub_ using pull requests.

.. _API: https://code-mphi.github.io/ECOGEN/docs/doxygen_docs/index.html
.. _Doxygen: http://www.doxygen.nl/
.. _`GitHub`: https://github.com/code-mphi/ECOGEN