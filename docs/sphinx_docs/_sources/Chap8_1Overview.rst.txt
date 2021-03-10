.. role:: xml(code)
  :language: xml

.. _Sec:GUI:overview:

********
Overview
********

Editing a test case requires deep understanding of input-file structure. To simplify this editing, a beta version of graphical user interface ECOGEN_GUI is proposed.

Possible actions are:

	* Open, edit and save existing test cases.
	* Create new test cases.
	* Generate the main input file *ECOGEN.xml* to prepare runs.
	* Edit essential settings for one or several test cases: Output format, flow model, fluid models, mesh dimensions (or mesh file choice), initial conditions and boundary conditions.
	* Visualize equation-of-state parameter library (editing these parameters is possible directely through the interface).

.. admonition:: Limitations

	Since this is a beta version of ECOGEN_GUI, only essential computation settings appear. Specific settings are thus not available but can be edited directly via the input files. Please refer to section :ref:`Chap:input` for details.

************
Installation
************

From self-extracting files
==========================

ECOGEN_GUI can be installed automatically under Windows 10 (64 bits) from the precompiled binary installer *ECOGEN_GUI_Setup.exe*.

To install the application, simply follow the installer's instructions (figure :numref:`Fig:GUI:install`).

.. _Fig:GUI:install:

.. figure:: ./_static/GUI/install.png
  :scale: 70%
  :align: center

  ECOGEN_GUI installation window under Windows 10.

It is preferable to install ECOGEN_GUI in a user folder (by default) so that the application can save its configuration files.

Source compilation
==================

The source files of the ECOGEN_GUI application can also be compiled. The compilation nevertheless requires an environment able to compile C++ code using the `Qt`_ library.

.. _`Qt`: https://www.qt.io/

********************
General presentation
********************

First-time running
==================

The first time you run ECOGEN_GUI, the configuration window appears (figure :numref:`Fig:GUI:config`.)). Afterwards, this window will be accessible via the menu bar of the main application window.

.. _Fig:GUI:config:

.. figure:: ./_static/GUI/config.png
  :scale: 70%
  :align: center

  Configuration window of the ECOGEN_GUI application.

Before being able to use ECOGEN_GUI, it is necessary to specify the directory in which the ECOGEN code is located, in order to make the association between the application and the code input files. 

.. Caution:: 

	The ECOGEN simulation tool must be installed independently. Please refer to the corresponding section :ref:`Chap:Start` for detailed instructions on ECOGEN installation.

The other fields (*Mesh software*, *Text editor software*, etc.) can be filled in optionally. They will then allow the launching of these applications directly from the ECOGEN_GUI interface.

Main window
===========

The main window of ECOGEN_GUI is presented, as visible at the launch of the application, on the figure :numref:`Fig:GUI:mainWindow`. 

.. _Fig:GUI:mainWindow:

.. figure:: ./_static/GUI/mainWindow.png
  :scale: 30%
  :align: center

  ECOGEN_GUI application main window.

It consists of several areas:

	1. Menu bar,
	2. **Run controls** dock,
	3. **Models** dock,
	4. **Mesh properties** dock,
	5. **Boundary conditions** dock,
	6. **Physical domains** dock,
	7. **Equations of state** dock.

Apart from the **Boundary conditions** dock, the other docks are greyed out when the application is launched and will only be activated when a test case is opened. To start editing a test case, choose *File* :math:`\rightarrow` *Open test case* (Ctrl+O) or *File* :math:`\rightarrow` *New test case* (Ctrl+N) from the menu bar (1).

A first example
---------------
As an example, open the reference test case located in *./libTests/referenceTestCases/euler/1D/transport/positiveVelocity*. The application then loads the parameters from the ECOGEN input files corresponding to the requested test case (figure :numref:`Fig:GUI:casTestRef`).

.. _Fig:GUI:casTestRef:

.. figure:: ./_static/GUI/casTestRef.png
  :scale: 30%
  :align: center

  Opening the reference test case located in *./libTests/referenceTestCases/euler/1D/transport/positiveVelocity*.

Loading the test case causes a new tab to appear in the central part of the window (8):

	* The tab title corresponds to the name of the test case being edited. It is also the name of the folder that will contain the results of the corresponding ECOGEN simulation (see section :ref:`Sec:input:main:runName` for more details).
	* The white zone of the tab will record all the operations performed on this test case (log).
	* If the test case parameters are modified, the tab title will be followed by a star which means that the modifications must be recorded to be applied to the ECOGEN code input files. To do this, select *File* :math:`\rightarrow` *Save test case* (Ctrl+S).

Multiple test case downloads and preparation of ECOGEN
------------------------------------------------------
It is possible to open multiple test cases simultaneously. In this case, as many tabs will be present in the central area (figure :numref:`Fig:GUI:zoomOnglets`).

.. _Fig:GUI:zoomOnglets:

.. figure:: ./_static/GUI/zoomOnglets.png
  :scale: 70%
  :align: center

  Here, 3 test cases are opened simultaneously. The test case being edited *pressureVelocityEq2DrichtmyerMeshkov* has been modified and not saved.

Once the test cases have been loaded into ECOGEN_GUI and all saved, it is possible to generate the main input file of the ECOGEN code from the menu bar (*Run* :math:`\rightarrow` *Prepare* or Ctrl+E). The simulations are ready to run (refer to section :ref:`Sec:installation:compileAndExecute`).

Reloading a test case
---------------------
At any time, the input files can be reloaded again within the ECOGEN_GUI interface via the menu command (*View* :math:`\rightarrow` *Refresh* or Ctrl+R). This command is particularly useful to cancel unsaved modifications of the parameters of a test case or during a "manual" modification of the input files.

*******************
Desription of docks
*******************

The main window of ECOGEN_GUI consists of docks that can be arranged according to the user's preferences. These different docks gather the settings of a test case by "family". A help on each parameter can be obtained via the appearance of a tooltip on mouse over.

Any change to a setting of a test case is accompanied by a line in the *log* field (central white area (8)) of the corresponding tab. This allows you to keep a history of the work performed on each test case.

.. Caution::

	Some parameters are not yet implemented within ECOGEN_GUI. The specific requirements for these parameters can nevertheless be edited manually within the input files. In this case, refer to section :ref:`Chap:UserGuide`.

**Run Controls** dock
=====================
The dock **Run controls** (2) gathers all the computation settings initially contained in the input file *mainV5.xml* of the considered test case.
The details of the parameters are presented in section :ref:`Sec:input:main`.

**Models** dock
===============
The dock **Models** (3) gathers all the parameters of the mathematical flow model initially contained in the input file *modelV4.xml* of the considered test case.
The details of the parameters are presented in section :ref:`Sec:input:model`. This dock allows, among other things, to:

	* Change the model on the fly.
	* Change the number of phases.
	* Select for each phase an equation of state present in the ECOGEN directory.

An interesting feature of ECOGEN_GUI is that any modification linked to the model will automatically be reported on the other docks (**Boundary conditions** and **Physical domains**) and consequently on the input files of the test case during saving.

**Mesh properties** dock
========================
The dock **Mesh properties** (4) gathers all the settings related to the mesh. These parameters are initially contained in the input file *meshV5.xml* of the considered test case.
The details of the parameters are presented in section :ref:`Sec:input:mesh`.

The appearance of this dock depends on the type of mesh considered (Cartesian or unstructured). The two possible aspects of the dock **Mesh properties** are presented on figure :numref:`Fig:GUI:meshDock`.

.. _Fig:GUI:meshDock:

.. figure:: ./_static/GUI/meshDock.png
  :scale: 70%
  :align: center

  Two different aspects of the dock **Mesh properties**. Left: The standalone Cartesian mesh version is fully adaptable. Right: The unstructured version which requires to specify a mesh file from a third-party application.

It is not possible to change the nature of the mesh on the fly. To change the mesh type, a new test case with the correct initial mesh structure must be created.

Cartesian mesh
--------------
When the mesh is Cartesian (figure :numref:`Fig:GUI:meshDock`, left), it is fully adaptable via ECOGEN_GUI (dimensions, number of computation cells, use of AMR for adaptative mesh refinement of discontinuities).

Unstructured mesh
-----------------
When the mesh is unstructured (figure :numref:`Fig:GUI:meshDock`, right), the mesh must first be prepared using a third-party mesh application. The resulting mesh file must be specified here. Details on the mesh files accepted by ECOGEN are available in section :ref:`Sec:tuto:generatingMeshes`.

.. Caution::

	In case of modification of the mesh file, it may be necessary to edit the boundary conditions of the input file *initialConditionsV4.xml* manually, these boundary conditions being dependent on the mesh file. This will be recalled in the *log* area of the test case if needed.

**Boundary conditions** dock
============================
The dock **Boundary conditions** (5) gathers all the settings related to the boundary conditions. These parameters are initially contained in the input file *initialConditionsV4.xml* of the considered test case.
The details of the parameters are presented in section :ref:`Sec:input:boundaryConditions`.

The dock will automatically adjust to the parameters of the other docks (**Mesh properties** and **Models**).
The different boundary conditions available are selectable within the dock list. Once the condition is selected within the list, its parameters can be modified.

**Physical domains** dock
=========================
The dock **Physical domains** (6) gathers all the settings related to the physical domains to initialize (initial conditions of the computation). These parameters are initially contained in the input file *initialConditionsV4.xml* of the considered test case.
The details of the parameters are presented in section :ref:`Sec:input:physicalDomains`.

In ECOGEN_GUI, the dock **Physical domains** will automatically adapt to the parameters present in the dock **Models**.

The **Physical domains** dock is used to define different initialization regions for the thermo-mechanical variables of fluids within the computational domain. These regions operate by accumulation and can therefore overlap each other. It is then possible to create a base region initializing the whole domain and then add as many domains as desired using the :math:`+` button.
Once added, the geometric domains can be removed using the :math:`-` button or moved up or down (the domain at the top of the list being the first initialized, the following ones will be overlapped in order).

**Equations of state** dock
===========================
This dock is independent. It indexes all the files of thermodynamic parameters of the fluids available in the directory of ECOGEN. These parameters can be directly modified within the ECOGEN interface by checking the box within the dock (used as a lock). Information on the equations of state can be found in section :ref:`Sec:IO:materials`.

*********************************
Modify ECOGEN_GUI's configuration
*********************************

At any time, ECOGEN_GUI can be reconfigured via the menu bar *Edit* :math:`\rightarrow` *Configure*). In this case, the user has access to the configuration window again (:numref:`Fig:GUI:config`).

ECOGEN's working directory
==========================
In this window, the path to ECOGEN's working directory must be correctly specified to ensure the operation of the ECOGEN_GUI interface. At the time of validation, if the chosen directory does not contain the main input file of the code *ECOGEN.xml*, an error message will appear and the user will again be prompted to modify the directory (see section :ref:`Sec:tuto:mainXML` for details on the main input file *ECOGEN.xml*).

Links to external tools
=======================
It is an option you can use in ECOGEN_GUI. It is possible to specify the link to external application executables which can be frequently used during a simulation session via ECOGEN. Once these links are effective, it will be possible to call the corresponding applications via the menu (1) of the main window (*Tools*).