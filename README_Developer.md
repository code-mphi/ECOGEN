*********************************
Developpers information - ECOGEN
*********************************

Variable and Classe names:
--------------------------
A) Variable name should began with lowercase letter and be self-understandable. Each new word began with a uppercase letter. 
Exemples : int myInteger; vector<double *> vectorOfDoublePointer; etc.

B) Class attribute should began by "m_". 
Exemples : int m_myInteger; double * m_doublePointer; etc.

C) Class name should began with an uppercase Letter.
Exemples : class MyClass;


Standart Doxygen comments
-------------------------
To automatise code documentation, comments should be insert as follow in header files.

A)Before class definition
//! \class     class name
//! \brief     brief descritption
//! \details    Detailed description
//!                  continue description

B)Before function and method prototypes
//! \brief      brief description
//! \details    Detailed description
//!                  continue description
//! \param      parameter name         description
//! \param      parameter name         description
//! \return     return description

C)class member descrption
type m_variable; //!< member description



Developper personnal comments - flags
-------------------------------------
	
	//Developper//KeyWord// comments, ex: "//FP//DEV// comment, description"

Key word list :      //DEV//    in developpement
		     //Q//      question to dig
		     //TODO//   should be done in the future
		     //VS2015// problem on Visual Studio 2015
                     //ERR//    error : to correct ASAP
                     //ID//     idea
                     //ICI//    Stop developement position
		     //VERIF//  to verify : is it needed ?
		     //TEST//   test : To delete ASAP

-----For each modification, a comment should be prepared to be included to the commit message for Git.
