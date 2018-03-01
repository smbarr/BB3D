BB3D
====
Blade Builder 3D - 3-Dimensional Blade Surface CAD Generator  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1186451.svg)](http://dx.doi.org/10.5281/zenodo.1186451)

BB3D is a utility to generate a surface model of a blade based on a table of cross-sectional airfoil data. The code is built on the OpenCascade CAD package. The surface model is exported as an IGES file.

Quick start
-----------
1. Download and install OpenCascade (>=7.0.0) https://www.opencascade.com/content/latest-release
2. Clone the repository (git clone https://github.com/szb5100/BB3D)
3. Edit the "data.dat" file or airfoil coordinate files
4. Execute "./bin/make_blade" to generate the IGES file (Takes about a minute depending on system)

data.dat
--------
The data.dat file contains the airfoil cross-section and spar location information needed to build the blade. Sample data.dat and corresponding airfoil cross-section files defining the NREL 5-MW blade [1] are included in the BB3d/data/ directory.

The first three lines specify the number of blade sections along the span, the number of lofts to perform and the number of spars, respectively.  
Next, for each span location, the span number and airfoil coordinate file name.  
After the airfoil coordinate filenames is a header line with column names.  
Then the cross-section data is listed for each span section.  
Finally, after another header line, the span locations are listed on separate lines as a percentage of chord.

The airfoil coordinate files contain the x and y coordinate data scaled between x=0.0 and x=1.0. The order of the file should be from leading edge->trailing edge for the top of the blade, and then trailing edge->leading edge for the bottom of the blade.

1. J. Jonkman, S. Butterfield, W. Musial, and G. Scott. Definition of a 5-mw reference wind turbine for offshore system development. Technical Report NREL/TP-500-38060, National Renewable Energy Laboratory, 2009. 
