
## How to build FEMAddOns


#### Prerequisites
* Version 11.3 or greater of Mathematica or Wolfram Desktop.
* Wolfram Workbench

#### Building FEMAddOns 
* Import FEMAddOns in Workbench. This needs to be done only once.

Build the documentation:
* In the FEMAddOns folder right click on docbuild.xml
  * Choose Run As...
  * Choose 2 Ant Build...
  * Deselect all and select *clean*
  * Run
* Repeat the above and choose *DistMesh* in stead of clean 
* Repeat the above and choose *DomainDecomposition* in stead of clean 

This will create a folder build, which will contain a folder FEMAddOns that contains the build documentation.

Copy the code:
* Run the script copyContent.sh once to copy the code into the package.

This will copy the source code into the build/FEMAddOns package.

Create a Paclet:
* Possibly increment the VersionNumber in PacletInfo.m
* Load the paclet manager: Needs["PacletManager`"] 
* Set the directory into the build folder: SetDirectory[...] 
* Check: FileNames[] should return {"FEMAddOns"}
* Create the paclet: PackPaclet["FEMAddOns"]

This will leave you with an path to a FEMAddOns-1.0.paclet in the build folder.
