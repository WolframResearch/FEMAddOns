
BeginPackage["FEMAddOns`", {
		"NDSolve`FEM`",
		"DistMesh`",
		"DomainDecomposition`",
		"FEMUtils`",
		"ImportMesh`"
}]


FEMAddOns`$FEMAddOnsInstallationDirectory::usage="";

Begin["`Private`"];

FEMAddOns`$FEMAddOnsInstallationDirectory = DirectoryName[ $InputFileName]

End[]

EndPackage[]

