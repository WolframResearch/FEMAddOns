(* Paclet Info File *)

Paclet[
    Name -> "FEMAddOns",
    Version -> "1.4.9",
    MathematicaVersion -> "11+",
    Description -> "Package provides additional finite element method functionality.",
    Extensions -> 
        {
            {"Kernel", Context -> "FEMAddOns`"}, 
            {"Kernel", Root -> "Applications/Kernel", Context -> 
                {"Applications`"}
            }, 
            {"Kernel", Root -> "DistMesh/Kernel", Context -> 
                {"DistMesh`"}
            }, 
            {"Kernel", Root -> "DomainDecomposition/Kernel", Context -> 
                {"DomainDecomposition`"}
            }, 
            {"Kernel", Root -> "FEMUtils/Kernel", Context -> 
                {"FEMUtils`"}
            }, 
            {"Kernel", Root -> "ImportMesh/Kernel", Context -> 
                {"ImportMesh`"}
            }, 
            {"Documentation", Root -> "Applications/Documentation"}, 
            {"Documentation", Root -> "DistMesh/Documentation"}, 
            {"Documentation", Root -> "DomainDecomposition/Documentation"}, 
            {"Documentation", Root -> "FEMUtils/Documentation", MainPage -> "Guides/FEMAddOns"},
            {"Documentation", Root -> "ImportMesh/Documentation"}

        }
]

