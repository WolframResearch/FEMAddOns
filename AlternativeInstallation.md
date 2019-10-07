
**Alternative Install** 

The FEMAddOns release comes in the form of a `.paclet` file, which contains the entire package and its documentation. Download the latest release from the [Github repo's releases page](https://github.com/WolframResearch/FEMAddOns/releases). To install, run the following command in the Wolfram Language:

    PacletInstall["/full/path/to/FEMAddOns.paclet"]

**Subsequent installs** 

If you have a previous version (>=1.3.0) of the FEMAddOns installed you can use

	Needs["FEMAddOns`"]
	UpdateFEMAddOns[]

to update to the latest version.


