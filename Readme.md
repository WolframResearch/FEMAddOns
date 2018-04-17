
# FEMAddOns for Wolfram Language

FEMAddOns is a package that provides additional functionality to the [Wolfram Language](https://www.wolfram.com/language/) build in Finite Element Method. FEMAddOns supports 11.3 and later versions of Wolfram Language deployments for the desktop, including [Wolfram Desktop](https://www.wolfram.com/desktop/) and [Mathematica](https://www.wolfram.com/mathematica/).

### Installing the FEMAddOns release

The FEMAddOns release comes in the form of a `.paclet` file, which contains the entire package and its documentation. Download the latest release from the [Github repo's releases page](https://github.com/WolframResearch/FEMAddOns/releases). To install, run the following command in the Wolfram Language:

    PacletInstall["/full/path/to/FEMAddOns.paclet"]

This will permanently install the FEMAddOns paclet. The Wolfram Language will always use the latest installed version of FEMAddOns. Installed versions can be enumerated using the command:

    PacletFind["FEMAddOns"]

And all versions can be uninstalled using the command:

    PacletUninstall["FEMAddOns"]

To make use of the documentation it may be necessary to restart.

### Using FEMAddOns

To access the documentation, open the notebook interface help viewer, and search for FEMAddOns. The first hit will be a summary page enumerating the most commonly used functions in FEMAddOns. 


### More...

See the following files for more information:

* [License.md](License.md) - FEMAddOns license
* [Contributing.md](Contributing.md) - Guidelines for contributing to FEMAddOns
* [HowToBuild.md](HowToBuild.md) - Instructions for building and debugging FEMAddOns
