<h1>Welcome to the Organism-Tissue Simulator software page</h1>

The Organism-Tissue Simulator is a C++ software for simulating biological systems. It
mainly focuses on systems with multiple cells/compartments, and includes
biochemical as well as mechanical and cell proliferation rules for numerically simulating the system dynamics. 

The software is currently divided into two separate codes:

<h2>Organism <a href="https://gitlab.com/slcu/teamHJ/organism"> [main git repository] </a>: </h2> used for simulating cells with geometries that can be described by a fixed number of variables and parameters. 
Examples are cells described by Spheres, Capped cylinders for rod-shaped cells, and Budding (yeast) shapes. It also can be used for simulations on geometries extracted from confocal data, where it uses information on cell volumes cell-cell connections (walls) and their areas. Example papers using these different approaches are:
<br></br>
Jönsson et al (2006) <i>PNAS</i>
<br></br>
Cho, Jönsson et al (2006) <i>PLoS Biology</i>
<br></br>
Jönsson and Levchenko (2005) <i>Multiscale Modeling and Simulation</i>
<br></br>
Gruel et al (2016) <i>Science Advances</i>
<br></br>
<br></br>
<h2>Tissue [this git repository]: </h2> used for simulating vertex-based cell geometries including finite element mechanical models. A cell wall (cell in 2.5D) is defined by a list of vertices defining a polygon. Example papers using this software are:
<br></br>
Hamant, Heisler, Jönsson et al (2008) <i>Science</i>
<br></br>
Bozorg et al (2014) <i>PLoS Comp Biol</i>
<br></br>
Bhatia et al (2016) <i>Current Biology</i>
<br></br>
<br></br>
<h2>Features</h2>

Some of the features of the Organism-Tissue software are:

* It includes several numerical solvers for the ordinary differential equations, also with noise, and for stochastic simulations (Organism).

* It includes a library of common biochemical, growth, division and mechanical rules.

* Models are defined within a specific model file - no need for recompilation when models are changed.

* The code is compartmentalised such that it should be easy for a programmer to define
additional user-specific reaction- division- etc. rules.

* It is and has been used to simulate growing plant tissue, and yeast and
bacterial cell colonies.

* It includes visualisation and analysis tools.

* it includes an optimisation environment where data from multiple mutants can be used to define a cost function (only Organism). 

<h2>Development</h2>

The software is mainly developed at the Computational Biology & Biological
Physics group at Lund University and the Computational Morphodynamics group at the Sainsbury Laboratory, University of Cambridge. Contact: henrik.jonsson@slcu.cam.ac.uk.

Developers (incomplete list): Henrik Jönsson, Patrik Sahlin, Pontus Melke, Pau Formosa-Jordan, Laura Brown, Pawel Krupinski, Behruz Bozorg 

<h2>Binaries</h2>

There are three main binaries. ''simulator'' which is simulating a single organism-tissue
model. For the Organism code there is an ''optimizer'' binary which is used to optimise model parameters towards a data template, and ''newman'' which is used to visualise model output.
For the Tissue code, the preferred viualisation tool is Paraview (http://www.paraview.org).

<i>Simulator</i>

The <i>simulator</i> uses three input files:

* model file: this file defines the model. Information of its format can be
found in the documentation of Organism::readModel() [for tissue it is Tissue::readModel()].

* init file: this file holds the initial variable values. Information of
its structure can be found at the documentation for Organism::readInit().

* solver parameter file: With this file parameters for the solver is
provided. Information of its structure can be found in
BaseSolver::getSolver().

These three files are required to be included among the command line arguments
and in correct order. In addition different flags to the simulator can be set,
either on the command line or in a file ''$HOME/.organism''. Information on these
additional flags can be found in the documentation for the namespace myConfig.

A common command line execution of the simulator binary is then:

simulator example.model example.init example.rk5 > example.data

which generates the system output in the file example.data. Note that the simulator binary is in the bin directory of the organism (or tissue) main directory. This needs to be in your PATH or you will need to specify the full path to the binary.

<i>Optimizer</i>

The <i>optimizer</i> uses three input files:

* model file: this file defines the model. Information of its format can be
found in the documentation of Organism::readModel().

* solver parameter file: In this file parameters for the solver is
provided. Information of its structure can be found in
BaseSolver::getSolver().

* optimizer parameter file: Information about the optimizer to be used is
given, including parameters to be optimized, template file(s) to be compared
to, and cost function to be used. There is more documentation on the file
content in the documentation for the class BaseOptimizer and its sub-classes.

<h2>Compilation</h2>

The code can be compiled by using make. A Makefile can be found in the 'src' directory (to compile: 'cd src; make' will generate the binaries in the bin directory) . The
code is ANSI-compatible and has been tested to compile on Linux (default for the provided Makefile), Mac OS X and
Windows (using Cygwin) platforms. It is dependent on boost and for plotting OpenGL. 

<h2>Documentation</h2>

Documentation is generated via Doxygen (http://www.doxygen.org) and a Makefile is found in the 'doc' directory. 'cd doc; make' will generate a html version in the doc/html directory and the introductory page is doc/html/index.html. PDF documentation can be generated by 'cd doc/latex; make' which generates doc/latex/refman.pdf. Other formats are also available via Doxygen.

<h2>Contact</h2>

The latest version of the code is avialiable via the 
<a href="https://gitlab.com/slcu/teamHJ/organism"> Organism </a> and 
<a href="https://gitlab.com/slcu/teamHJ/tissue"> Tissue </a> gitlab repositories. 
Main contact: henrik.jonsson@slcu.cam.ac.uk.