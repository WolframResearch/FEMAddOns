(* Paclet Info File *)

(* created 2017/11/6*)

Paclet[
Name -> "FEMAddOns",
Version -> "1.0",
MathematicaVersion -> "11+",
Extensions -> {
	{"Kernel", Context -> "FEMAddOns`"},
	{"Kernel", Root -> "DistMesh/Kernel", Context -> {"DistMesh`"}},
	{"Kernel", Root -> "DomainDecomposition/Kernel", Context -> {"DomainDecomposition`"}},

	{"Documentation", Root -> "DistMesh/Documentation"},
	{"Documentation", Root -> "DomainDecomposition/Documentation"},
	{"Documentation", Root -> "FEMUtils/Documentation", MainPage -> "Guides/FEMAddOns"}
}
]


