(* ::Section:: *)
(* Package Setup *)
BeginPackage["DomainDecomposition`", {"NDSolve`FEM`"}]

ElementMeshDecomposition;
DecompositionNDSolveValue;
SetMaxCellMeasure;
Subdomain;
SubdomainGroup;

Begin["`Private`"]

(* ::Section:: *)
(* Data Structures *)

(*
The domain decomposition method code uses three custom data structures, MethodData, Subdomain, and SubdomainGroup. MethodData takes the form MethodData[{pdespec}, u, vars]. 
It encapsulates the problem information to make this readily available throughout the code. There are three different functions that can be used to extract
data from this data structure, e.g. MethodDataDependentVariables[MethodData[...]] will extract the vars components.

The Subdomain data structure takes the form Subdomain[mesh, artificialBoundaryIncidents, sharedIncidents, localToGlobal, solutions, bcs]. The mesh is the mesh corresponding to
the subdomain. The artificial boundary incidents are the new boundary nodes created by splitting the global mesh, these are the nodes on which the Schwarz method
imposes Dirichlet conditions. Values corresponding to sharedIncidents are needed by other subdomains to update their Dirichlet conditions. localToGlobal is a mapping
that takes incidents in the subdomain mesh, and converts that incident to the corresponding incident in the global mesh. solutions is a list of lists, where each list
is the solution to the problem over the subdomain mesh. There is one solution for each independent variable in the problem. bcs is the boundary conditions imposed on the
artificial boundary incidents.

The sharedIncidents are not necessary, and in fact that is not used for the parallel implementation anymore since the convergence criterion requires getting all values
from each subkernel in each iteration. It was originally used to not have to assemble the full solution in each iteration. (Now it is still used in some places, including
the tutorial, so it cannot be removed without going through the code.)

There are several functions available to extract information from a Subdomain object. There are also functions that update various components. For example, updates
boundary conditions, or updates the solutions.

SubdomainGroup takes the form SubdomainGroup[mesh, subdomains], where mesh is the mesh of the global problem. ElementMeshDecomposition returns a SubdomainGroup expression, this can
then be passed directly into DecompositionNDSolveValue.
*)

CreateSubdomain[mesh_, boundaryIncidents_, sharedIncidents_, localToGlobal_] := Subdomain[mesh, boundaryIncidents, sharedIncidents, localToGlobal, {}, {}]
CreateMethodData[{pdespec__}, u_, vars_] := MethodData[{pdespec}, Flatten[{u}], Flatten[{vars}]]

SubdomainMesh[Subdomain[mesh_, _, _, _, _, _]] := mesh
SubdomainCoordinates[Subdomain[mesh_, _, _, _, _, _]] := mesh["Coordinates"]
SubdomainGlobalIncidents[subdomain_Subdomain] := DDMCache[SubdomainLocalToGlobal[subdomain] /@ SubdomainLocalIncidents[subdomain], SubdomainGlobalIncidents, Hash[SubdomainMesh[subdomain]]]
SubdomainLocalIncidents[subdomain:Subdomain[mesh_, _, _, _, _, _]] := DDMCache[GetElementIncidents[mesh], SubdomainLocalIncidents, Hash[SubdomainMesh[subdomain]]]
SubdomainGlobalInteriorIncidents[subdomain_Subdomain] := DDMCache[Complement[SubdomainGlobalIncidents[subdomain], SubdomainGlobalArtificialBoundaryIncidents[subdomain]], SubdomainGlobalInteriorIncidents, Hash[SubdomainMesh[subdomain]]]
SubdomainLocalInteriorIncidents[subdomain_Subdomain] := DDMCache[Complement[SubdomainLocalIncidents[subdomain], SubdomainLocalArtificialBoundaryIncidents[subdomain]], SubdomainLocalInteriorIncidents, Hash[SubdomainMesh[subdomain]]]
SubdomainGlobalBoundaryIncidents[subdomain_Subdomain] := DDMCache[SubdomainLocalToGlobal[subdomain] /@ SubdomainLocalBoundaryIncidents[subdomain], SubdomainGlobalBoundaryIncidents, Hash[SubdomainMesh[subdomain]]]
SubdomainLocalBoundaryIncidents[subdomain: Subdomain[mesh_, _, _, _, _, _]] := DDMCache[GetBoundaryElementIncidents[mesh], SubdomainLocalBoundaryIncidents, Hash[SubdomainMesh[subdomain]]]
SubdomainLocalArtificialBoundaryIncidents[Subdomain[_, artificialBoundaryIncidents_, _, _, _, _]] := artificialBoundaryIncidents;
SubdomainGlobalArtificialBoundaryIncidents[subdomain_Subdomain] := DDMCache[SubdomainLocalToGlobal[subdomain] /@ SubdomainLocalArtificialBoundaryIncidents[subdomain], SubdomainGlobalArtificialBoundaryIncidents, Hash[SubdomainMesh[subdomain]]];
SubdomainLocalSharedIncidents[Subdomain[_, _, sharedIncidents_, _, _, _]] := sharedIncidents
SubdomainGlobalSharedIncidents[subdomain_Subdomain] := DDMCache[SubdomainLocalToGlobal[subdomain] /@ SubdomainLocalSharedIncidents[subdomain], SubdomainGlobalSharedIncidents, Hash[SubdomainMesh[subdomain]]]
SubdomainGlobalToLocal[subdomain:Subdomain[_, _, _, localToGlobal_, _, _]] := DDMCache[GeneralUtilities`AssociationInvert[localToGlobal], subdomain, Hash[SubdomainMesh[subdomain]]]
SubdomainLocalToGlobal[Subdomain[_, _, _, localToGlobal_, _, _]] := localToGlobal
SubdomainLocalSolutions[Subdomain[_, _, _, _, solutions_, _]] := solutions
SubdomainBoundaryConditions[Subdomain[_, _, _, _, _, bcs_]] := bcs

SetLocalSolutions[Subdomain[mesh_, artificialBoundaryIncidents_, sharedIncidents_, localToGlobal_, _, bcs_], solutions_] := Subdomain[mesh, artificialBoundaryIncidents, sharedIncidents, localToGlobal, solutions, bcs]
SetBoundaryConditions[Subdomain[mesh_, artificialBoundaryIncidents_, sharedIncidents_, localToGlobal_, solutions_, _], bcs_] := Subdomain[mesh, artificialBoundaryIncidents, sharedIncidents, localToGlobal, solutions, bcs]

SubdomainGroupMesh[SubdomainGroup[mesh_, _]] := mesh
SubdomainGroupSubdomains[SubdomainGroup[_, subdomains_]] := subdomains

MethodDataPDESpecification[MethodData[{pdespec__}, _, _]] := pdespec
MethodDataDependentVariables[MethodData[_, u_, _]] := u
MethodDataIndependentVariables[MethodData[_, _, vars_]] := vars

FormatSubdomainGroup[SubdomainGroup[mesh_, subdomains_]] := SubdomainGroup@BoxForm`SurroundWithAngleBrackets[Tooltip[Length[subdomains], "Number of subdomains"]];
SubdomainGroup /: Format[subdomainGroup_SubdomainGroup] := FormatSubdomainGroup[subdomainGroup]
SubdomainGroup /: MakeBoxes[subdomainGroup_SubdomainGroup, fmt_] := BoxForm`SuppressedBoxForm[FormatSubdomainGroup[subdomainGroup], subdomainGroup, fmt]

FormatSubdomain[Subdomain[mesh_, __]] := Subdomain[mesh["Bounds"]];
Subdomain /: Format[subdomain_Subdomain] := FormatSubdomain[subdomain]
Subdomain /: MakeBoxes[subdomain_Subdomain, fmt_] := BoxForm`SuppressedBoxForm[FormatSubdomain[subdomain], subdomain, fmt]

(* The following are properties, for users who may want to implement their own domain decomposition code *)
Subdomain[mesh_, _, _, _, _, _]["Mesh"] := mesh
Subdomain[mesh_, _, _, _, _, _]["Coordinates"] := mesh["Coordinates"]
Subdomain[mesh_, _, _, _, _, _]["Incidents"] := GetElementIncidents[mesh]
Subdomain[args__]["InteriorIncidents"] := Complement[SubdomainLocalIncidents[Subdomain[args]], SubdomainLocalArtificialBoundaryIncidents[Subdomain[args]]]
Subdomain[mesh_, _, _, _, _, _]["BoundaryIncidents"] := GetBoundaryElementIncidents[mesh]
Subdomain[_, artificialBoundaryIncidents_, _, _, _, _]["ArtificialBoundaryIncidents"] := artificialBoundaryIncidents
Subdomain[_, _, sharedIncidents_, _, _, _]["SharedIncidents"] := sharedIncidents
Subdomain[_, _, _, localToGlobal_, _, _]["LocalToGlobal"] := localToGlobal
Subdomain[_, _, _, _, solutions_, _]["Solutions"] := solutions
Subdomain[_, _, _, _, _, bcs_]["BoundaryConditions"] := bcs

Subdomain[args__]["GlobalIncidents"] := With[{s = Subdomain[args]}, s["LocalToGlobal"] /@ s["Incidents"]]
Subdomain[args__]["GlobalInteriorIncidents"] := With[{s = Subdomain[args]}, s["LocalToGlobal"] /@ s["InteriorIncidents"]]
Subdomain[args__]["GlobalBoundaryIncidents"] := With[{s = Subdomain[args]}, s["LocalToGlobal"] /@ s["BoundaryIncidents"]]
Subdomain[args__]["GlobalArtificialBoundaryIncidents"] := With[{s = Subdomain[args]}, s["LocalToGlobal"] /@ s["ArtificialBoundaryIncidents"]]
Subdomain[args__]["GlobalSharedIncidents"] := With[{s = Subdomain[args]}, s["LocalToGlobal"] /@ s["SharedIncidents"]]

Subdomain[__]["Properties"] := Sort[{"Mesh", "Coordinates", "Incidents", "InteriorIncidents", "BoundaryIncidents", "ArtificialBoundaryIncidents", "SharedIncidents", "LocalToGlobal", "Solutions", "BoundaryConditions"}] 

SubdomainGroup[mesh_, _]["Mesh"] := mesh
SubdomainGroup[_, subdomains_]["Subdomains"] := subdomains
SubdomainGroup[__]["Properties"] = Sort[{"Mesh", "Subdomains", "Wireframe"}]

SubdomainGroup[mesh_ElementMesh, subdomains : {__Subdomain}]["Wireframe"] := SubdomainGroup[mesh, subdomains]["Wireframe"[]]
SubdomainGroup[mesh_ElementMesh, subdomains : {__Subdomain}]["Wireframe"[pos:{__Integer}]] := SubdomainGroup[mesh, Part[subdomains, pos]]["Wireframe"[]]

SubdomainGroup[_ElementMesh, subdomains : {__Subdomain}]["Wireframe"[opts___]] := Module[
	{meshes, colors, meshElementStyle, defaults, optionNames, optionValues, wireframeOptions, graphicsOptions},
	
	meshes = SubdomainMesh /@ subdomains;
	colors = PadRight[{}, Length[subdomains], ColorData[97, "ColorList"]];
	meshElementStyle = Directive[FaceForm[{Opacity[0.4], #}], EdgeForm[{Opacity[0.4], Darker[#]}]] & /@ colors;
	
	defaults = {
		"MeshElementStyle" -> meshElementStyle,
		"MeshElementIDStyle" -> {None},
		"MeshElementMarkerStyle" -> {None},
		"MeshElement" -> "MeshElements"
	};
	
	optionNames = {"MeshElementStyle", "MeshElementIDStyle", "MeshElementMarkerStyle", "MeshElement"};
	optionValues = optionNames /. Join[{opts}, defaults];
	optionValues = PadRight[{}, Length[subdomains], #] & /@ optionValues;
	wireframeOptions = MapThread[Thread[# -> #2] &, {optionNames, optionValues}];
	graphicsOptions = FilterRules[{opts}, Options[Graphics]];
	Show[MapThread[#["Wireframe"[##2]] &, {meshes, Sequence @@ wireframeOptions}], graphicsOptions]
]

(* DDMCache serves to cache certain computations that are done repeatedly *)
(*SetAttributes[DDMCache, HoldFirst];
DDMCache[action_, symbol_, hash_] := DDMCache[_, symbol, hash] = action*)

(* Disabled caching because of problems *)
DDMCache[action_, _, _] := action

(* ::Section:: *)
(* ElementMeshDecomposition *)
Options[ElementMeshDecomposition] = Join[{
	"Subdomains" -> 4,
	"Overlap" -> Automatic
}];

ElementMeshDecomposition::overlap = "Overlap must be a non-negative integer or a real value between 0 and 1.";
ElementMeshDecomposition::noverlap = "Number of overlaps specified `1` does not match the number of subdomains `2`.";
ElementMeshDecomposition::subdomains = "The subdomains specification `1` is not a valid subdomains specification for the given mesh.";

ElementMeshDecomposition[mesh_ElementMesh, opts : OptionsPattern[]] := Block[{
		overlap, elements, nrOfElements, undirectedEdges, sa, graph, subdomains, coordinates, meshElements, globalToLocal, localToGlobal,
		subdomainMeshes, incidents, boundaryIncidents, artificialBoundaryIncidents, sharedIncidents, layers
	},
	
	elements = Join @@ ElementIncidents[mesh["MeshElements"]];
	nrOfElements = Length[elements];

	(* Creates the graph without unpacking *)
	undirectedEdges = GetUndirectedEdge[mesh["ContinuousElementConnectivity"]];
	sa = SparseArray[undirectedEdges -> Pattern, {nrOfElements, nrOfElements}];
    graph = Graph[Automatic, {Null, sa}];

    (* The "Subdomains" option can take any form that can be used by FindGraphPartition. If
    FindGraphPartition fails, this code returns $Failed. *)
	subdomains = Check[
		Quiet@FindGraphPartition[graph, OptionValue["Subdomains"]],
		Message[ElementMeshDecomposition::subdomains, OptionValue["Subdomains"]];
		Return[$Failed]
	];
	If[
		Head[subdomains] === FindGraphPartition,
		Message[ElementMeshDecomposition::subdomains, OptionValue["Subdomains"]];
		Return[$Failed]
	];
	subdomains = DeleteCases[subdomains, {}];

	overlap = OptionValue["Overlap"];
	If[overlap === Automatic, overlap = ComputeOverlap[#, undirectedEdges, 0.1] & /@ subdomains];
	
	(* If overlap is given as a list of fractions *)
	If[MatchQ[overlap, _Real?(0 < # < 1 &)], overlap = ConstantArray[overlap, Length[subdomains]]];
	If[
		MatchQ[overlap, {__Real?(0 < # < 1 &)}] && Length[overlap] == Length[subdomains],
		overlap = MapThread[ComputeOverlap[#, undirectedEdges, #2] &, {subdomains, overlap}];
	];
	
	(* Overlap can be given as an integer or a list of integer. *)
	If[MatchQ[overlap, Except[_List]],overlap = ConstantArray[overlap, Length[subdomains]]];
	If[
		! MatchQ[overlap, {__Integer?(# >= 0 &)}],
		Message[ElementMeshDecomposition::overlap];
		Return[$Failed]
	];
	If[
		Length[overlap] != Length[subdomains],
		Message[ElementMeshDecomposition::noverlap, Length[overlap], Length[subdomains]];
		Return[$Failed]
	];
	(* For comment about layers, see the function ExtendSubdomain. *)
	{layers, subdomains} = Transpose@MapThread[ExtendSubdomain[mesh["ContinuousElementConnectivity"], #, #2] &, {subdomains, overlap}];
	
	coordinates = Part[mesh["Coordinates"], Union @@ Part[elements, #]] & /@ subdomains;
	meshElements = SelectMeshElements[mesh["MeshElements"], #] & /@ subdomains;

	(* For comment about reindexing, see the function ReindexMeshElements *)
	{globalToLocal, meshElements} = Transpose[ReindexMeshElements /@ meshElements];
	localToGlobal = GeneralUtilities`AssociationInvert /@ globalToLocal;
  	
  	subdomainMeshes = MapThread[
		ToElementMesh["Coordinates" -> #, "MeshElements" -> #2] &,
		{coordinates, meshElements}
	]; (*TODO what about edge element markers? point element markers? *)

	incidents = LocalToGlobal[GetElementIncidents /@ subdomainMeshes, localToGlobal];
	boundaryIncidents = LocalToGlobal[GetBoundaryElementIncidents /@ subdomainMeshes, localToGlobal];

	artificialBoundaryIncidents = ArtificialBoundaryIncidents[layers, mesh, {incidents, boundaryIncidents}, {localToGlobal, globalToLocal}];

	(* This code finds the incidents in each subdomain which are used to set Dirichlet conditions in other subdomains. *)
	sharedIncidents = ConstantArray[{}, Length[incidents]];
	Do[
		sharedIncidents[[i]] = Intersection[incidents[[i]], Union @@ Delete[artificialBoundaryIncidents, i]];
		, {i, Length[incidents]}
     ];
	sharedIncidents = GlobalToLocal[sharedIncidents, globalToLocal];	
	
	artificialBoundaryIncidents = GlobalToLocal[artificialBoundaryIncidents, globalToLocal];

	subdomains = MapThread[CreateSubdomain, {subdomainMeshes, artificialBoundaryIncidents, sharedIncidents, localToGlobal}];
	SubdomainGroup[mesh, subdomains]
 ]
 
 
 
(* ::Section:: *)
(* DecompositionNDSolveValue *)

Options[DecompositionNDSolveValue] = {
	"Tolerance" -> 10^(-6),
	"MaxIterations" -> 100,
	"InitialGuess" -> None,
	"Method" -> Automatic,
	"EvaluationMonitor" -> None,
	"Kernels" -> None,
	"DomainDecomposition" -> {
		"MeshOptions" -> {},
		"Subdomains" -> 4,
		"Overlap" -> 0.1
	},
	"LinearSolveMethod" -> "Pardiso"
};

DecompositionNDSolveValue::method = "`` is not one of the available methods (MultiplicativeSchwarz, AdditiveSchwarz)";
DecompositionNDSolveValue::maxiter = "Failed to converge to the requested accuracy or precision within `` iterations."
DecompositionNDSolveValue::nkernels = "The number of kernels `` does not match the number of subdomains ``.";
DecompositionNDSolveValue::tmkernels = "There are `` kernels and `` subdomains. All kernels will not be used.";

DecompositionNDSolveValue[{pdespec__}, u_, Element[vars_, region: (_SubdomainGroup | _ElementMesh | _?RegionQ)], opts : OptionsPattern[]] := Module[
	{subdomains, nrOfSubdomains, dof, guess, globalSolutions, md, evaluationMonitor, i, multiplicative, additive, 
	previousGlobalSolutions, sequential, parallel, kernels, partialSolutions, solutions, mesh, iterations, restrictedPDEspec, maxIterations, subdomainGroup,
	domainDecomposition, meshOptions, subdomainOptions, linearSolveMethod},
	
	If[
		!MatchQ[OptionValue["Method"], Automatic | "MultiplicativeSchwarz" | "AdditiveSchwarz"],
		Message[DecompositionNDSolveValue::method, OptionValue["Method"]]	
	];
	
	domainDecomposition = DeleteDuplicatesBy[Join[OptionValue["DomainDecomposition"], "DomainDecomposition" /. Options[DecompositionNDSolveValue]], First];
	meshOptions = FilterRules["MeshOptions" /. domainDecomposition, Options[ToElementMesh]];
	subdomainOptions = FilterRules[domainDecomposition, Options[ElementMeshDecomposition]];

	If[MatchQ[OptionValue["Kernels"], _List], PrependTo[subdomainOptions, "Subdomains" -> Length[OptionValue["Kernels"]]]];

	If[
		MatchQ[region, _SubdomainGroup],
		subdomainGroup = region,
		If[
			MatchQ[region, _ElementMesh],
			subdomainGroup = ElementMeshDecomposition[region, Sequence @@ subdomainOptions],
			subdomainGroup = ElementMeshDecomposition[ToElementMesh[region, Sequence @@ meshOptions], Sequence @@ subdomainOptions]
		]
	];
	If[subdomainGroup === $Failed, Return[$Failed]];

	mesh = SubdomainGroupMesh[subdomainGroup];
	dof = Length[mesh["Coordinates"]];
	
	subdomains = SubdomainGroupSubdomains[subdomainGroup];
	nrOfSubdomains = Length[subdomains];
	
	guess = OptionValue["InitialGuess"];
	If[guess == None, guess = ConstantArray[0 &, Length[Flatten@{u}]]];
	subdomains = SetApproximateSolutions[subdomains, Flatten@{guess}];
	
	globalSolutions = ConstantArray[0, {Length[Flatten@{u}], dof}];
	globalSolutions = UpdateGlobalSolutions[subdomains, globalSolutions];
	
	restrictedPDEspec = RestrictBoundaryConditions[{pdespec}, vars, ToBoundaryMesh[mesh]];
	
	md = CreateMethodData[restrictedPDEspec, u, vars];
	subdomains = NewBoundaryConditions[subdomains, md];
	subdomains = UpdateBoundaryConditions[subdomains, globalSolutions];
	subdomains = SetLocalSolutions[#, globalSolutions[[All, SubdomainGlobalIncidents@#]]] & /@ subdomains;
	
	evaluationMonitor = OptionValue["EvaluationMonitor"];
	maxIterations = OptionValue["MaxIterations"];
	linearSolveMethod = OptionValue["LinearSolveMethod"];
	
	multiplicative := (
		iterations = 0;
		While[
			iterations < maxIterations,
			Do[
				subdomains[[i]] = UpdateBoundaryConditions[subdomains[[i]], globalSolutions];
				subdomains[[i]] = UpdateLocalSolutions[subdomains[[i]], md, linearSolveMethod];
				globalSolutions = UpdateGlobalSolutions[subdomains[[i]], globalSolutions];
				,
				{i, nrOfSubdomains}
			];
			
			previousGlobalSolutions = globalSolutions;
			globalSolutions = GetGlobalSolutions[subdomains, globalSolutions];
			If[
				evaluationMonitor =!= None,
				evaluationMonitor[globalSolutions, previousGlobalSolutions];
			];
			If[Max[Norm /@ (globalSolutions - previousGlobalSolutions)] < OptionValue["Tolerance"], Break[]];
			
			iterations++;
		];
		If[iterations >= maxIterations, Message[DecompositionNDSolveValue::maxiter, maxIterations]];
	);
	
	additive := (
		iterations = 0;
		While[
			iterations < maxIterations,
			previousGlobalSolutions = globalSolutions;
			Do[
				subdomains[[i]] = UpdateBoundaryConditions[subdomains[[i]], previousGlobalSolutions];
				subdomains[[i]] = UpdateLocalSolutions[subdomains[[i]], md, linearSolveMethod];
				globalSolutions = UpdateGlobalSolutions[subdomains[[i]], globalSolutions];
				,
				{i, nrOfSubdomains}
			];
			
			globalSolutions = GetGlobalSolutions[subdomains, globalSolutions];
			If[
				evaluationMonitor =!= None,
				evaluationMonitor[globalSolutions, previousGlobalSolutions];
			];
			If[Max[Norm /@ (globalSolutions - previousGlobalSolutions)] < OptionValue["Tolerance"], Break[]];
			
			iterations++
		];
		If[iterations >= maxIterations, Message[DecompositionNDSolveValue::maxiter, maxIterations]];
	);
		
	sequential := (
		Switch[
			OptionValue["Method"],
			"AdditiveSchwarz", additive,
			"MultiplicativeSchwarz" | Automatic | _, multiplicative
		];
		GetGlobalSolutions[subdomains, globalSolutions]
	);
	
	parallel := (
		iterations = 0;
		MapThread[InitializeKernel[#, md, globalSolutions, SubdomainLocalSharedIncidents[#], #2, linearSolveMethod] &, {subdomains, kernels}];
		While[
			iterations < maxIterations,
			previousGlobalSolutions = globalSolutions;	
			partialSolutions = ParallelEvaluate[Iterate, kernels];
			MapThread[(globalSolutions[[All, #]] = #2) &, {SubdomainGlobalInteriorIncidents /@ subdomains, partialSolutions}];
			Map[UpdateGlobalSolutionsOnKernel[globalSolutions, #] &, kernels];
			
			If[
				evaluationMonitor =!= None,
				evaluationMonitor[globalSolutions, previousGlobalSolutions];
			];
			If[Max[Norm /@ (globalSolutions - previousGlobalSolutions)] < OptionValue["Tolerance"] && iterations > 2, Break[]];
			
			iterations++;
		];
		If[iterations >= maxIterations, Message[DecompositionNDSolveValue::maxiter, maxIterations]];
		GetGlobalSolutionsFromKernels[globalSolutions, kernels]
	);
	
	(* If kernels were specified, use the parallel implementation. Otherwise use a sequential implementation. *)
	kernels = OptionValue["Kernels"];
	solutions = If[
		MatchQ[kernels, {__Parallel`Kernels`kernel}],
		If[
			Length[kernels] == Length[subdomains],
			parallel, 
			If[
				Length[kernels] > Length[subdomains],
				Message[DecompositionNDSolveValue::tmkernels, Length[kernels], Length[subdomains]];
				kernels = Take[kernels, Length[subdomains]];
				parallel,
				Message[DecompositionNDSolveValue::nkernels, Length[kernels], Length[subdomains]]; Return[$Failed]
			]
		],
		sequential
	];
	
	(* Following the practice of NDSolveValue *)
	If[
		Length[solutions] == 1,
		ElementMeshInterpolation[{mesh}, First[solutions]],
		ElementMeshInterpolation[{mesh}, #] & /@ solutions	
	]
]

(* ::Section:: *)
(* Auxiliary functions primarily for ElementMeshDecomposition *)

GetUndirectedEdge = Compile[{{cec, _Integer, 2}}, 
	Block[{rows, cols, this, index, vals},
		{rows, cols} = Dimensions[cec];
		index = Internal`Bag[Most[{0}]];
		vals = Internal`Bag[Most[{0}]];
		Do[
			Do[
				this = Compile`GetElement[cec, r, c];
				If[this == 0, Continue[]];
				Internal`StuffBag[index, r];
				Internal`StuffBag[vals, Abs[this]];
				,
				{c, cols}
			];
			,
			{r, rows}
		];
		Sort /@ Transpose[{Internal`BagPart[index, All], Internal`BagPart[vals, All]}]
	]
];

(* One "layer" consists of all elements connected to a subdomain, which are not currently part of it. A subdomain can be extended
by increasing its size recursively one layer at a time. ExtendSubdomain returns a list of n layers and the original subdomain.*)
(* The argument "connectivity" should be a ContinuousElementConnectivity list. *)
ExtendSubdomain[_, subdomain_, 0] := {{}, subdomain};
ExtendSubdomain[connectivity_, subdomain_, n_] := Module[{layers},
	layers = NestWhile[NextLayer[connectivity, #] &, {{}, subdomain}, UnsameQ[#, {}] &, 1, n];
	{layers, Union @@ layers}
]

ElementMeshDecomposition::size = "The subdomain including overlap is the same as the original mesh.";
NextLayer[connectivity_, previousLayers_] := Module[{neighbors},
	neighbors = Union @@ Abs@connectivity[[Last@previousLayers]]; (* All elements connected to the previous layer. *)
	If[neighbors == {}, Message[ElementMeshDecomposition::size]; Return[previousLayers]];
	neighbors = If[First[neighbors] === 0, Rest[neighbors], neighbors];
	Append[previousLayers, Complement[neighbors, Sequence @@ previousLayers[[{-1, -2}]]]]
]

SelectMeshElements[elems : {__?MeshElementQ}, part : {__?IntegerQ}] := Block[{toBeSelected},
	toBeSelected = GetMeshElementParts[ElementIncidents[elems], part];
	Select[MapThread[MeshElementByPart, {elems, toBeSelected}], MeshElementQ]	
]

(* Gets the parts as required by MeshElementByPart, defined in the NDSolve`FEM package. *)
GetMeshElementParts[incidents_, subdomain1_] := Block[
	{subdomain, subdomainLength, bags, continuousIndex, currentIncident},
	subdomain = Sort[subdomain1];
	subdomainLength = Length[subdomain];
	bags = Table[Internal`Bag[Most[{0}]], {Length[incidents]}];
	continuousIndex = 1;
	currentIncident = 1;
	Do[
		Do[
 			If[
				continuousIndex == subdomain[[currentIncident]],
				Internal`StuffBag[bags[[i]], j];
				If[(currentIncident += 1) > subdomainLength, Break[]];
			];
			continuousIndex++;
			, 
			{j, Length@incidents[[i]]}
		];
		If[currentIncident > subdomainLength, Break[]];
		, 
		{i, Length@incidents}
	];
	Internal`BagPart[#, All] & /@ bags
]

(* If the original mesh has incidents {1, 2, 3, 4, 5}, then a subdomain may e.g. have the incidents {1, 2, 5}. To be compatible
with ToElementMesh, {1, 2, 5} has to be mapped to {1, 2, 3} (a continuous sequence of integers starting at 1) before creating
the subdomain mesh. This is what ReindexMeshElements does. The map {1 -> 1, 2 -> 2, 3 -> 5} is stored so that the solution from the 
subdomain can be used to reconstruct the solution of the original problem later. *)
ReindexMeshElements[elems : {__?MeshElementQ}] := Block[{incidents, incidentsUnion, dict, reindexed, markers},
	incidents = ElementIncidents /@ elems;
	incidentsUnion = Union@Flatten[incidents];
	dict = AssociationThread[incidentsUnion -> Range[Length[incidentsUnion]]];
	reindexed = Map[dict, #, {2}] & /@ incidents;
	{dict, MapThread[ReindexMeshElement, {elems, reindexed}]}
]

ReindexMeshElement[el_?MeshElementQ, reindexedIncidents_] := Block[{markers},
	markers = Sequence[];
	If[ElementMarkersQ[el], markers = ElementMarkers[el]];
	MeshElementCreate[Head[el], reindexedIncidents, markers]
]

GetElementIncidents[mesh_] := Union@Flatten@ElementIncidents[mesh["MeshElements"]]
GetBoundaryElementIncidents[mesh_] := Union@Flatten@ElementIncidents[mesh["BoundaryElements"]]

LocalToGlobal[incidents : {__List}, localToGlobal : {__Association}] := MapThread[LocalToGlobal, {incidents, localToGlobal}]
LocalToGlobal[incidents_List, localToGlobal_Association] := Map[localToGlobal, incidents]
GlobalToLocal[incidents : {__List}, globalToLocal : {__Association}] := MapThread[GlobalToLocal, {incidents, globalToLocal}]
GlobalToLocal[incidents_List, globalToLocal_Association] := Map[globalToLocal, incidents]

(* OverlappingBoundaryIncidents selects the boundary incidents of one mesh which also belong to another mesh. In other words,
it selects the boundary incidents which lie in the overlap between two or more elements. *)
OverlappingBoundaryIncidents[incidents_, boundaryIncidents_] := Module[{permutations, complements},
	permutations = RotateLeft@Partition[incidents, Length[incidents] - 1, 1, {1, 1}];
	permutations = MapThread[Prepend, {permutations, boundaryIncidents}];
	complements = Complement @@@ permutations;
	MapThread[Complement, {boundaryIncidents, complements}]
]

(* The artificial boundary is made up of the incidents in the subdomain mesh on which Schwarz methods impose Dirichlet conditions. 
This boundary is called artificial since it mostly lies in the interior of the original, global mesh. *)
ArtificialBoundaryIncidents[layers_, mesh_, {incidents_, boundaryIncidents_}, {localToGlobal_, globalToLocal_}] := Module[
	{elements, elementIncidents, selectedIncidents, overlappingBoundaryIncidents},
	If[
		Length[First[layers]] < 4, (* no overlap *)
		Return[OverlappingBoundaryIncidents[incidents, boundaryIncidents]]
	];
	elements = Apply[Union, layers[[All, {-1}]], {1}];
	elementIncidents = Join @@ ElementIncidents[mesh["MeshElements"]];
	selectedIncidents = Map[Union @@ Part[elementIncidents, #] &, elements];
	overlappingBoundaryIncidents = OverlappingBoundaryIncidents[incidents, boundaryIncidents];
	selectedIncidents = MapThread[Intersection, {overlappingBoundaryIncidents, selectedIncidents}]
]

 SetMaxCellMeasure[{mesh_, symbolicRegion_}, maxCellMeasure_] := Module[
 	{bmeshcoarse, nrcoarse},
 	bmeshcoarse = ToBoundaryMesh[mesh];
 	nrcoarse = ToNumericalRegion[symbolicRegion];
 	SetNumericalRegionElementMesh[nrcoarse, bmeshcoarse];
 	ToElementMesh[nrcoarse, MaxCellMeasure -> maxCellMeasure, "SteinerPoints" -> False]
 ]

(* computeOverlap takes a subdomain and suggests what an appropriate amount of overlap might be,
by using the graph diameter of the element connectivity graph. *)
ComputeOverlap[subdomain_, undirectedEdges_, beta_] := Module[{graph, edges},
	edges = Pick[undirectedEdges, And @@@ (lookup[subdomain] /@ undirectedEdges), True];
	graph = Graph[UndirectedEdge @@@ edges];
    	Ceiling[beta 2 GraphDiameter[graph, Method -> "PseudoDiameter"]]
];
(* Multiplies diameter by two because PseudoDiamer returns half the value of Dijkstra and Johnson.*)
(* PseudoDiameter is necessary because the other methods consume orders of magnitude more memory. *)

lookup[subdomain_] := With[{asc = AssociationThread[subdomain -> True]},
   Function[edges, Block[{Missing},
     Missing["KeyAbsent", _] = False;
     asc /@ edges
     ]]
   ];

(*
There is a potential one-liner:
ComputeOverlap[graph_, subdomain_, beta_] := Ceiling[beta 2 GraphDiameter[Subgraph[graph, subdomain], Method -> "PseudoDiameter"]]
*)

(* ::Section:: *)
(* Parallel Additive Schwarz implementation *)
(*UpdateGlobalSolutionsOnKernel[solutions_, kernel_] := ParallelEvaluate[globalSolutions[[All, boundaryIncidents]] = solutions, kernel]*)
UpdateGlobalSolutionsOnKernel[solutions_, kernel_] := ParallelEvaluate[globalSolutions = solutions, kernel]
  
InitializeKernel[subdomain0_, methodData_, globalSolutions0_, partialSolutionsIncidents_, kernel_, linearSolveMethod_] := (
	DistributeDefinitions[SubdomainGlobalIncidents, SubdomainLocalArtificialBoundaryIncidents, UpdateBoundaryConditions, UpdateLocalSolutions, SubdomainLocalSolutions, SubdomainBoundaryConditions,SubdomainLocalInteriorIncidents,SubdomainGlobalInteriorIncidents];
	ParallelEvaluate[
		subdomain = subdomain0;
		globalSolutions = globalSolutions0;
		KernelUpdateBoundaryConditions := subdomain = UpdateBoundaryConditions[subdomain, globalSolutions];
		KernelUpdateLocalSolutions := subdomain = UpdateLocalSolutions[subdomain, methodData, linearSolveMethod];
		KernelLocalSolutions := Part[SubdomainLocalSolutions[subdomain], All, SubdomainLocalInteriorIncidents[subdomain]];
		KernelGlobalIncidents := SubdomainGlobalIncidents[subdomain];
		KernelGlobalInteriorIncidents := SubdomainGlobalInteriorIncidents[subdomain];
		KernelBoundaryConditions := SubdomainBoundaryConditions[subdomain];
		
		Iterate := (
			KernelUpdateBoundaryConditions;
			KernelUpdateLocalSolutions;
			KernelLocalSolutions
		);
		,
		kernel,
		DistributedContexts -> None
	]
);

GetGlobalSolutionsFromKernels[globalSolutions_, kernels_] := Module[{incidents, values, newGlobalSolutions},
	newGlobalSolutions = globalSolutions;
	incidents = ParallelEvaluate[KernelGlobalInteriorIncidents, kernels];
	values = ParallelEvaluate[KernelLocalSolutions, kernels];
	MapThread[(newGlobalSolutions[[All, #]] = #2) &, {incidents, values}];
	newGlobalSolutions
]

(* ::Section:: *)
(* Auxiliary Functions for DecompositionNDSolveValue *)

SetApproximateSolutions[subdomain_Subdomain, coarseGlobalSolutions_] := Module[{coordinates, solutions, localIncidents},
	localIncidents = SubdomainLocalIncidents[subdomain];
	coordinates = SubdomainCoordinates[subdomain];
	solutions = Apply[#, coordinates, {1}] & /@ coarseGlobalSolutions;
	SetLocalSolutions[subdomain, solutions[[All, localIncidents]]]
]
SetApproximateSolutions[subdomains : {__Subdomain}, coarseGlobalSolution_] := SetApproximateSolutions[#, coarseGlobalSolution] & /@ subdomains

(* This updates precisely those components in the global solutions which are required by other subdomains to update their Dirichlet conditions. *)
UpdateGlobalSolutions[subdomain_Subdomain, globalSolutions_] := Module[{newGlobalSolutions, globalIncidents, localSolutions, localIncidents},
	globalIncidents = SubdomainGlobalSharedIncidents[subdomain];
	localIncidents = SubdomainLocalSharedIncidents[subdomain];
	localSolutions = SubdomainLocalSolutions[subdomain][[All, localIncidents]];
	newGlobalSolutions = globalSolutions;
	newGlobalSolutions[[All, globalIncidents]] = localSolutions;
	newGlobalSolutions
]
UpdateGlobalSolutions[subdomains : {__Subdomain}, globalSolution_] := Fold[UpdateGlobalSolutions[#2, #] &, globalSolution, subdomains]

(* This merges all solutions in all subdomains, which results in the global solution *)
GetGlobalSolutions[subdomain_Subdomain, globalSolutions_] := Module[{newGlobalSolutions, globalIncidents, localSolutions, localIncidents},
	globalIncidents = SubdomainGlobalInteriorIncidents[subdomain];
	localIncidents = SubdomainLocalInteriorIncidents[subdomain];
	localSolutions = SubdomainLocalSolutions[subdomain][[All, localIncidents]];
	newGlobalSolutions = globalSolutions;
	newGlobalSolutions[[All, globalIncidents]] = localSolutions;
	newGlobalSolutions
]
GetGlobalSolutions[subdomains : {__Subdomain}, globalSolution_] := Fold[GetGlobalSolutions[#2, #] &, globalSolution, subdomains]

(* Creates a template for every node on the artificial boundary, e.g. DirichletCondition[u[x,y] = 0, "ElementIncidents" -> {n}] *)
NewBoundaryConditions[subdomain_Subdomain, md_MethodData] := Module[{bcs, localIncidents, u, vars},
	localIncidents = SubdomainLocalArtificialBoundaryIncidents[subdomain];
	u = MethodDataDependentVariables[md];
	vars = MethodDataIndependentVariables[md];
	bcs = Outer[DirichletCondition[#[Sequence @@ vars] == 0, "ElementIncidents" -> {{#2}}] &, u, localIncidents, 1];
	SetBoundaryConditions[subdomain, bcs]
]
NewBoundaryConditions[subdomains : {__Subdomain}, md_MethodData] := NewBoundaryConditions[#, md] & /@ subdomains

(* Updates the templates given by NewBoundaryCondition by setting a new value for the Dirichlet condition *)
UpdateBoundaryConditions[subdomains : {__Subdomain}, globalSolutions_] := UpdateBoundaryConditions[#, globalSolutions] & /@ subdomains
UpdateBoundaryConditions[subdomain_Subdomain, globalSolutions_] := Module[{bcs, newBCs, globalBoundaryIncidents},
	globalBoundaryIncidents = SubdomainGlobalArtificialBoundaryIncidents[subdomain];
	bcs = SubdomainBoundaryConditions[subdomain];
	newBCs = MapThread[UpdateBoundaryConditions, {bcs, globalSolutions[[All, globalBoundaryIncidents]]}];
	SetBoundaryConditions[subdomain, newBCs]
]
UpdateBoundaryConditions[{}, _] := {}
UpdateBoundaryConditions[bcs : {__DirichletCondition}, newValues_] := MapThread[UpdateBoundaryConditions, {bcs, newValues}]
UpdateBoundaryConditions[DirichletCondition[lhs_ == rhs_, pred_], newValue_] := DirichletCondition[lhs == newValue, pred]

UpdateLocalSolutions[subdomain_Subdomain, md_MethodData, linearSolveMethod_:Automatic] := Module[{solutions, pdespec, bcs, u, vars, mesh},
	pdespec = MethodDataPDESpecification[md];
	bcs = Sequence @@ Flatten[SubdomainBoundaryConditions[subdomain]];
	u = MethodDataDependentVariables[md];
	vars = MethodDataIndependentVariables[md];
	mesh = SubdomainMesh[subdomain];
	If[
		linearSolveMethod === Automatic,
		solutions = Quiet[NDSolveValue[{pdespec, bcs}, u, Element[vars, mesh]], {NDSolveValue::bcnop}],
		solutions = Quiet[NDSolveValue[{pdespec, bcs}, u, Element[vars, mesh], Method -> {"PDEDiscretization" -> {"FiniteElement", "LinearSolveMethod" -> {linearSolveMethod}}}], {NDSolveValue::bcnop}];
	];
	(*solutions = {ndsolve[{pdespec, bcs}, u, Element[vars, mesh]]};*)
	(*solutions = {ndsolve2[{"Poisson", {bcs}}, u, Element[vars, mesh]]};*)
	solutions = #["ValuesOnGrid"] & /@ solutions;
	SetLocalSolutions[subdomain, solutions]
]

(* If a user specifies a boundary condition such as DirichletCondition[u[x,y] = 0, x > 5] this will also affect the artificial boundary if
nothing is done about it. RestrictBoundaryConditions modifies the condition so that it only applies to the physical boundary. *)
RestrictBoundaryConditions[pdespec_, vars_, subregion_] := Module[{meshRegion, regionMemberQ},
	meshRegion = MeshRegion[subregion];
	regionMemberQ = RegionMember[meshRegion];
	ReplaceAll[
		pdespec,
		DirichletCondition[bc_, pred_] :> DirichletCondition[bc, And[regionMemberQ[vars], pred]]
	]
]

(* Used to compare timings. Should be removed when no longer needed. *)
ndsolve[{pdespec__}, u_, Element[vars_, mesh_], 
  opts : OptionsPattern[]] := Module[
  {state, femdata, initCoeffs, methodData, sd, discretePDE, initBCs,discreteBCs, load, stiffness, damping, mass, sol},
  {state} = NDSolve`ProcessEquations[{pdespec}, u, Element[vars, mesh], Method -> {"PDEDiscretization" -> "FiniteElement"}];
  femdata = state["FiniteElementData"];
  
  initCoeffs = femdata["PDECoefficientData"];
  methodData = femdata["FEMMethodData"];
  sd = methodData["SolutionData"];
  discretePDE = DiscretizePDE[initCoeffs, methodData, sd];
  
  initBCs = femdata["BoundaryConditionData"];
  discreteBCs = DiscretizeBoundaryConditions[initBCs, methodData, sd];
  
  {load, stiffness, damping, mass} = discretePDE["SystemMatrices"];
  DeployBoundaryConditions[{load, stiffness}, discreteBCs];
  
  sol = LinearSolve[stiffness, load, Method -> "Pardiso"];
  ElementMeshInterpolation[{mesh}, sol]
]

(* Used to compare timings. Should be removed when no longer needed. *)
ndsolve2[{"Poisson", artificialBCs_}, u_, Element[vars_, mesh_], opts : OptionsPattern[]] := Module[
	{nr, vd, sd, mdata, cdata, bcdata, dpde, dbcs, l, s, d, m, sol, dep},
	
	dep = Flatten[{u}];
	
	nr = ToNumericalRegion[mesh];
	vd = NDSolve`VariableData[{"DependentVariables", "Space"} -> {dep, vars}];
	sd = NDSolve`SolutionData["Space" -> nr];
	
	mdata = InitializePDEMethodData[vd, sd];
	cdata = InitializePDECoefficients[vd, sd, "DiffusionCoefficients" -> {{-IdentityMatrix[2]}}, "LoadCoefficients" -> {{1}}];
	bcdata = InitializeBoundaryConditions[vd, sd, {Join[{DirichletCondition[Through[dep @@ vars][[1]] == 0., Global`x == 0 || Global`x == 5 || Global`y == 0 || Global`y == 1]}, artificialBCs]}];
	
	dpde = DiscretizePDE[cdata, mdata, sd];
	dbcs = DiscretizeBoundaryConditions[bcdata, mdata, sd];
	
	{l, s, d, m} = dpde["All"];
	DeployBoundaryConditions[{l, s}, dbcs];
	
	sol = LinearSolve[s, l];
	ElementMeshInterpolation[{mesh}, sol]
]

(* ::Section:: *)
(* Package Postscript *)
End[]
EndPackage[]
