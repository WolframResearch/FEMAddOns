
## How to build FEMAddOns


#### Prerequisites
* Version 11.3 or greater of Mathematica or Wolfram Desktop.
* Wolfram Workbench

#### Building FEMAddOns 
First, import FEMAddOns in Workbench:
* Select Import...
* Git - Projects from Git (Next)
* Existing local repository (Next)
* Add... (browse to FEMAddOns, select, Next)
* Import as general project (Finish)

The importing of the FEMAddOns needs to be done only once.


Next, build the documentation:
* In the FEMAddOns folder right click on docbuild.xml
  * Choose Run As...
  * Choose 2 Ant Build...
  * Deselect all 
  * Select *clean*
  * Run
* Repeat the above and choose *DistMesh* instead of clean 
* Repeat the above and choose *DomainDecomposition* instead of clean 
* Repeat the above and choose *FEMUtils* instead of clean 

This will create a folder named build, which will contain a folder FEMAddOns that contains the build documentation of package.

Copy the code into the package:
* Run the script copyContent.sh once to copy the code into the package.

This will copy the source code into the build/FEMAddOns package.

Last, create the Paclet:
* Possibly increment the VersionNumber in PacletInfo.m
* Load the paclet manager: Needs["PacletManager`"] 
* Set the directory into the build folder: SetDirectory[...] 
* Check: FileNames[] should return {"FEMAddOns"}
* Create the paclet: PackPaclet["FEMAddOns"]

This will leave you with a path to a FEMAddOns-1.0.paclet in the build folder.

