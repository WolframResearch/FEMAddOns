
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

Copy source code and create and install the paclet:
* Open CopyContent.nb
* Follow instructions to copy the source code
* Follow instructions to create and install the paclet

This will leave you with a path to a FEMAddOns-X.Y.paclet in the build folder.
