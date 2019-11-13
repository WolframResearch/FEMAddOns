(* ::Package:: *)

BeginPackage["FEMUtils`", {"NDSolve`FEM`"}];


ClearAll[ElementMeshSmoothing];
ElementMeshSmoothing::usage = "ElementMeshSmoothing[mesh] smoothes an ElementMesh.";
Options[ElementMeshSmoothing] = {Method -> Automatic};

ClearAll[StructuredMesh];
StructuredMesh::usage="StructuredMesh[raster,{nx,ny}] creates structured mesh of quadrilaterals.
StructuredMesh[raster,{nx,ny,nz}] creates structured mesh of hexahedra.";

ClearAll[ToQuadMesh]
ToQuadMesh::usage="ToQuadMesh[mesh] converts triangular mesh to quadrilateral mesh.";

ClearAll[UpdateFEMAddOns]
UpdateFEMAddOns::usage="UpdateFEMAddOns[] downloads and installes the latest version of the FEMAddOns paclet.";

ClearAll[
	BoundaryElementMeshDifference,
	BoundaryElementMeshIntersection,
	BoundaryElementMeshJoin,
	BoundaryElementMeshRotation,
	BoundaryElementMeshSymmetricDifference,
	BoundaryElementMeshTranslate,
	BoundaryElementMeshUnion
]

BoundaryElementMeshDifference::usage="";
BoundaryElementMeshIntersection::usage="";
BoundaryElementMeshJoin::usage="";
BoundaryElementMeshRotation::usage="";
BoundaryElementMeshSymmetricDifference::usage="";
BoundaryElementMeshTranslate::usage="";
BoundaryElementMeshUnion::usage="";



Begin["`Private`"];

pack = Developer`ToPackedArray;

$defaultMaxIter = 150;
ClearAll[IterativeElementMeshSmoothing];
Options[IterativeElementMeshSmoothing] = {MaxIterations -> $defaultMaxIter};

mkFindConnectedNodes[mesh_ElementMesh] := 
Block[
	{vec, vecTvec, nzv, nzp, neigbours, temp2, sa}, 

	vec = mesh["VertexElementConnectivity"];
	vecTvec = vec.Transpose[vec];
	nzv = vecTvec["NonzeroValues"];
	nzp = vecTvec["NonzeroPositions"];
	neigbours = Join @@ Position[nzv, 2];
	temp2 = Transpose[nzp[[neigbours]]];
	sa = SparseArray[ Transpose[ Developer`ToPackedArray[{temp2[[1]],
		Range[Length[temp2[[1]]]]}]] -> 1,
			{Length[mesh["Coordinates"]], Length[neigbours]}];

	With[{lookup = sa, ids = temp2[[2]]},
		Function[theseNodes, 
			ids[[Flatten[lookup[[#]]["NonzeroPositions"]]]] & /@ theseNodes]
	]
]

findInteriorNodes[mesh_ElementMesh] := 
With[
	{allNodes = Range[Length[mesh["Coordinates"]]],
	boundaryNodes = 
		Join @@ (Flatten /@ ElementIncidents[mesh["BoundaryElements"]])
	},
	Complement[allNodes, boundaryNodes]
]


IterativeElementMeshSmoothing[m_ElementMesh,
	opts : OptionsPattern[IterativeElementMeshSmoothing]] := 
Block[
	{inodes, findConnectedNodes2, temp, movedCoords, oldCoords, iter},

	iter = OptionValue["MaxIterations"];
	If[ !(IntegerQ[iter] && iter > 0), iter = $defaultMaxIter];
  
	inodes = findInteriorNodes[m];
	findConnectedNodes2 = mkFindConnectedNodes[m];
	temp = findConnectedNodes2[inodes];
	oldCoords = m["Coordinates"];
	Do[
		movedCoords = Developer`ToPackedArray[Mean[oldCoords[[#]]] & /@ temp];
		oldCoords[[inodes]] = movedCoords;
	, {iter}];

	ToElementMesh["Coordinates" -> oldCoords,
		"MeshElements" -> m["MeshElements"],
		"BoundaryElements" -> m["BoundaryElements"],
		"PointElements" -> m["PointElements"],
		"RegionHoles" -> m["RegionHoles"]
	]
]


DirectElementMeshSmoothing[mesh_, 
	opts : OptionsPattern[DirectElementMeshSmoothing]] := 
Block[
	{n, vec, mat, adjacencymatrix2, mass2, laplacian2, typoOpt, 
	bndvertices2, interiorvertices, stiffness, load, newCoords}, 

	n = Length[mesh["Coordinates"]];
	vec = mesh["VertexElementConnectivity"];
	mat = Unitize[vec.Transpose[vec]];
	vec = Null;
	adjacencymatrix2 = mat - DiagonalMatrix[Diagonal[mat]];
	mass2 = DiagonalMatrix[SparseArray[Total[adjacencymatrix2, {2}]]];
	stiffness = N[mass2 - adjacencymatrix2];
	adjacencymatrix2 = Null;
	mass2 = Null;

	bndvertices2 =  Flatten[Join @@ ElementIncidents[mesh["PointElements"]]];
	interiorvertices = Complement[Range[1, n], bndvertices2];

	stiffness[[bndvertices2]] = IdentityMatrix[n, SparseArray][[bndvertices2]];

	load = ConstantArray[0., {n, mesh["EmbeddingDimension"]}];
	load[[bndvertices2]] = mesh["Coordinates"][[bndvertices2]];

	newCoords = LinearSolve[stiffness, load];

	typoOpt = If[$VersionNumber <= 11.3,
			"CheckIncidentsCompletness" -> False
			,
			"CheckIncidentsCompleteness" -> False
		];

	ToElementMesh["Coordinates" -> newCoords, 
		"MeshElements" -> mesh["MeshElements"], 
		"BoundaryElements" -> mesh["BoundaryElements"], 
		"PointElements" -> mesh["PointElements"], 
		typoOpt,
		"CheckIntersections" -> False, 
		"DeleteDuplicateCoordinates" -> False,
		"RegionHoles" -> mesh["RegionHoles"]
	]
]


(*
https://mathematica.stackexchange.com/questions/156611/improve-the-mesh-smoothing-procedure
*)

ElementMeshSmoothing[mesh_ElementMesh, 
	opts : OptionsPattern[ElementMeshSmoothing]] := 
Block[
	{method, fun, subOpts},

	method = Flatten[{OptionValue["Method"]}];

	Switch[method[[1]],
		"Iterative",
			fun = IterativeElementMeshSmoothing;
		,_,
			(* default, includes "Direct" *)
			fun = DirectElementMeshSmoothing;
	 ];

	subOpts = {};
	If[ Length[method] > 1, subOpts = method[[2 ;; -1]]	];

	fun[mesh, FilterRules[subOpts, Options[fun]]]
]




StructuredMesh::array = "Raster of input points must be a full array of numbers with depth of `1`.";
Options[StructuredMesh] = {InterpolationOrder->1, "MeshOrder" -> 1};
SyntaxInformation[StructuredMesh]={"ArgumentsPattern"->{_,_,OptionsPattern[]}};


getElementConnectivity[nx_,ny_]:=Flatten[
	Table[{
		j+(i-1)(ny+1),
		j+i(ny+1),
		j+i(ny+1)+1,
		j+(i-1)(ny+1)+1
		},
		{i,1,nx},
		{j,1,ny}
	],
	1
];

getElementConnectivity[nx_,ny_,nz_]:=Flatten[
	Table[{
		k+(j-1)(nz+1)+(i-1)(nz+1)(ny+1),
		k+(j-1)(nz+1)+i(nz+1)(ny+1),
		k+j(nz+1)+i(nz+1)(ny+1),
		k+j(nz+1)+(i-1)(nz+1)(ny+1),
		k+(j-1)(nz+1)+(i-1)(nz+1)(ny+1)+1,
		k+(j-1)(nz+1)+i(nz+1)(ny+1)+1,
		k+j(nz+1)+i(nz+1)(ny+1)+1,
		k+j(nz+1)+(i-1)(nz+1)(ny+1)+1
		},
		{i,1,nx},
		{j,1,ny},
		{k,1,nz}
	],
	2
];


unitStructuredMesh[nx_,ny_]:=With[{
	(* Vectorized (the fastest) way to get coordinates of unit square. *)
	nodes=Flatten[
		Outer[List,Subdivide[0.,1.,nx],Subdivide[0.,1.,ny]],
		1
	],
	connectivity=getElementConnectivity[nx,ny]
	},
	(* We can disable all error checks for this mesh, because we know it is correct. *)
	ToElementMesh[
		"Coordinates"->nodes,
		"MeshElements"->{QuadElement[connectivity]},
		"CheckIntersections"->False,
		"CheckQuality"->False,
		"DeleteDuplicateCoordinates"->False
	]
];

unitStructuredMesh[nx_,ny_,nz_]:=With[{
	nodes=Flatten[
		Outer[List,Subdivide[0.,1.,nx],Subdivide[0.,1.,ny],Subdivide[0.,1.,nz]],
		2
	],
	connectivity=getElementConnectivity[nx,ny,nz]
	},
	ToElementMesh[
		"Coordinates"->nodes,
		"MeshElements"->{HexahedronElement[connectivity]},
		"CheckIntersections"->False,
		"CheckQuality"->False,
		"DeleteDuplicateCoordinates"->False
	]
];


(* We create structured mesh on unit square (or cube) and then re-interpolate nodal 
coordinates over given raster points. Benefit is that it is easier to change "MeshOrder" 
or do some other mesh manipulation (refinement of edges). *)

StructuredMesh[raster_, {nx_, ny_}, opts:OptionsPattern[]
	] /; Length[raster] > 1 :=
Module[
    {order, dim, restructured, xInt, yInt, zInt, unitMesh, unitCrds, crds},
    If[ Not @ ArrayQ[raster, 3, NumericQ],
		Message[StructuredMesh::array, 3+1];
		Return[$Failed, Module]
	];

    order = OptionValue[InterpolationOrder] /. Automatic->1;
    dim = Last @ Dimensions[raster];
    
	(* It makes sense to work only with machine precison numbers from here on. *)
    restructured = Transpose[N@raster, {3, 2, 1}];
    xInt = ListInterpolation[restructured[[1]], {{0.,1.},{0.,1.}}, InterpolationOrder->order];
    yInt = ListInterpolation[restructured[[2]], {{0.,1.},{0.,1.}}, InterpolationOrder->order];

	unitMesh = MeshOrderAlteration[
		unitStructuredMesh[nx,ny],
		OptionValue["MeshOrder"]/.(Except[1|2]->1)
	];
	
	unitCrds = unitMesh["Coordinates"];
	crds=If[
		dim == 3,
		zInt = ListInterpolation[restructured[[3]], {{0.,1.},{0.,1.}}, InterpolationOrder->order];
		Transpose@{xInt@@@unitCrds, yInt@@@unitCrds, zInt@@@unitCrds},
		Transpose@{xInt@@@unitCrds, yInt@@@unitCrds}
	];
	
	(* We still check for quality and duplicate nodes, because of possible mesh distortions. *)
	If[
		dim == 3,
		ToBoundaryMesh[
			"Coordinates" -> crds,
			"BoundaryElements" -> unitMesh["MeshElements"]
		],
		ToElementMesh[
			"Coordinates" -> crds,
			"MeshElements" -> unitMesh["MeshElements"]
		]
	]
];

StructuredMesh[raster_, {nx_, ny_, nz_}, opts:OptionsPattern[]
	] /; Length[raster] > 1 :=
Module[
    {order, restructured, xInt, yInt, zInt, unitMesh, unitCrds, crds},
    If[ Not @ ArrayQ[raster, 4, NumericQ],
		Message[StructuredMesh::array, 4+1];
		Return[$Failed, Module]
	];

    order = OptionValue[InterpolationOrder] /. Automatic->1;
       
    restructured = Transpose[N@raster, {4, 3, 2, 1}];
    xInt = ListInterpolation[restructured[[1]], {{0.,1.},{0.,1.},{0.,1.}}, InterpolationOrder->order];
    yInt = ListInterpolation[restructured[[2]], {{0.,1.},{0.,1.},{0.,1.}}, InterpolationOrder->order];
    zInt = ListInterpolation[restructured[[3]], {{0.,1.},{0.,1.},{0.,1.}}, InterpolationOrder->order];
    
	unitMesh=MeshOrderAlteration[
		unitStructuredMesh[nx,ny,nz],
		OptionValue["MeshOrder"]/.(Except[1|2]->1)
	];
	
	unitCrds=unitMesh["Coordinates"];
	crds=Transpose@{xInt@@@unitCrds, yInt@@@unitCrds, zInt@@@unitCrds};
	
	ToElementMesh[
		"Coordinates" -> crds,
		"MeshElements" -> unitMesh["MeshElements"]
	]
];


(*
Source of algorithm for mesh conversion is:
	Houman Borouchaki, Pascal J. Frey;
	Adaptive triangular\[Dash]quadrilateral mesh generation;
	International Journal for Numerical Methods in Engineering;
	1998; Vol. 41, p. (915-934)

Original implementation is by Prof. Joze Korelc;
taken from AceFEM package (http://symech.fgg.uni-lj.si/), 

*)

ToQuadMesh::elmtype = "Only conversion of pure triangular meshes is supported."

distortion := distortion = Compile[
		{{n1, _Real, 1}, {n2, _Real, 1}, {n3, _Real, 1}, {n4, _Real, 1}},
		2/Pi * Max @ Map[
			Abs[Pi/2 - If[ # < 0, Pi - #, #] ]& @ 
				ArcTan[	#[[1, 1]] #[[2, 1]] + #[[1, 2]] #[[2, 2]],
						#[[1, 1]] #[[2, 2]] - #[[1, 2]] #[[2, 1]]
				]&,
			{{n2-n1, n4-n1}, {n3-n2, n1-n2}, {n4-n3, n2-n3}, {n1-n4, n3-n4}}
			]
	];

ToQuadMesh[meshIn_] /; ElementMeshQ[meshIn] && !BoundaryElementMeshQ[meshIn] &&
	meshIn["EmbeddingDimension"] === 2 :=
Module[
	{coor, econn, elem, dist, ncoor, taken, quad, triag, edge, nc, allquads, 
	alltriag, marker, nodes, markerQ, mesh, increaseOrderQ},

	mesh = meshIn;

	If[ mesh["MeshOrder"] =!= 1,
		mesh = MeshOrderAlteration[mesh, 1]];
		increaseOrderQ = True;

	If[ Union[Head /@ mesh["MeshElements"]] =!= {TriangleElement},
		Message[ToQuadMesh::elmtype]; Return[$Failed, Module]
	];

	coor = mesh["Coordinates"];
	econn = Join @@ mesh["ElementConnectivity"];
	(* zero is used as a default marker *)
	marker = Join @@ ElementMarkers[ mesh["MeshElements"]];
	markerQ = ElementMarkersQ[mesh["MeshElements"]];
	elem = Join @@ ElementIncidents[ mesh["MeshElements"]];	
	
	dist = MapThread[ (
		{
		If[ #2[[1]] == 0,
			Nothing,
			distortion[Sequence@@coor[[nodes = {#[[1]], #[[2]], Complement[elem[[#2[[1]]]], #1][[1]], #[[3]]}]]]->{#3, #2[[1]], nodes}
		],
		If[ #2[[2]] == 0,
			Nothing,
			distortion[Sequence@@coor[[nodes = {#[[1]], #[[2]], #[[3]], Complement[elem[[#2[[2]]]], #1][[1]]}]]]->{#3, #2[[2]], nodes}
		],
		If[ #2[[3]] == 0,
			Nothing,
			distortion[Sequence@@coor[[nodes = {#[[1]], Complement[elem[[#2[[3]]]], #1][[1]], #[[2]], #[[3]]}]]]->{#3, #2[[3]], nodes}
		]
		})&,
		{elem, econn, Range[elem//Length]}
	]//Flatten;

	taken = ConstantArray[False, elem//Length];
	quad = Map[
		If[ 
			Or@@taken[[#[[2, {1, 2}]]]]  || #[[1]]>0.8 || marker[[#[[2, 1]]]]  =!= marker[[#[[2, 2]]]],
			Nothing,
			taken[[#[[2, {1, 2}]]]] = True;{#[[2, 3]], marker[[#[[2, 1]]]]}
		]&,
		Sort[dist]
	]//Transpose;
	
	dist = Null;
	triag = {Extract[elem, #], Extract[marker, #]}& @ Position[taken, False];
	edge = DeleteDuplicates[
		Join[
			Flatten[Map[{#[[{1, 2}]]//Sort, #[[{2, 3}]]//Sort, #[[{3, 4}]]//Sort, #[[{4, 1}]]//Sort} &, quad[[1]]], 1],
			Flatten[Map[{#[[{1, 2}]]//Sort, #[[{2, 3}]]//Sort, #[[{3, 1}]]//Sort} &, triag[[1]]], 1]
		]
	];
	
	ncoor = coor//Length;
	nc = Dispatch[Flatten[Map[(ncoor++;{#->ncoor, #[[{2, 1}]]->ncoor})&, edge]]];
	
	allquads = MapThread[
		(ncoor++; 
		{Total[coor[[#]]]/4,
			{
			{#[[1]], #[[{1, 2}]]/.nc, ncoor, #[[{4, 1}]]/.nc},
			{#[[{1, 2}]]/.nc, #[[2]], #[[{2, 3}]]/.nc, ncoor},
			{ncoor, #[[{2, 3}]]/.nc, #[[3]], #[[{3, 4}]]/.nc},
			{#[[{4, 1}]]/.nc, ncoor, #[[{3, 4}]]/.nc, #[[4]]}
			},
		{#2, #2, #2, #2}
		})& ,
		quad
	];
	
	alltriag = MapThread[
		(ncoor++;
		{Total[coor[[#]]]/3,
			{
			{#[[1]], #[[{1, 2}]]/.nc, ncoor, #[[{3, 1}]]/.nc},
			{#[[{1, 2}]]/.nc, #[[2]], #[[{2, 3}]]/.nc, ncoor},
			{ncoor, #[[{2, 3}]]/.nc, #[[3]], #[[{3, 1}]]/.nc}
			},
		{#2, #2, #2}})&, 
		triag
	];

	mesh = ToElementMesh[
				"Coordinates" -> Join[
					coor,
					Map[(coor[[#[[1]]]]+coor[[#[[2]]]])/2&, edge],
					allquads[[All, 1]], alltriag[[All, 1]]
				], 
				"MeshElements" -> If[ markerQ,
					{QuadElement[Flatten[Join[
						allquads[[All, 2]],
						alltriag[[All, 2]]], 1],
						Flatten[{allquads[[All, 3]], alltriag[[All, 3]]}]]},
					{QuadElement[Flatten[Join[
						allquads[[All, 2]], alltriag[[All, 2]]], 1]]}
				],
				"RegionHoles" -> meshIn["RegionHoles"]
			];

	If[ !ElementMeshQ[mesh], Return[ $Failed, Module]];

	If[ increaseOrderQ, mesh = MeshOrderAlteration[mesh, meshIn["MeshOrder"]]];

	mesh
]


(* first seen at https://github.com/szhorvat/IGraphM#installation *)
UpdateFEMAddOns[] :=
Module[
	{json, download, target, msg},
	Check[
		json = Import[
			"https://api.github.com/repos/WolframResearch/FEMAddOns/releases"
			, "JSON"];
		download = Lookup[ First[ Lookup[ First[json], "assets"]],
			"browser_download_url"];
		msg = "Downloading FEMAddOns " <> Lookup[First[json], "tag_name"] <> " ...";
		target = FileNameJoin[{CreateDirectory[], "FEMAddOns.paclet"}];

		If[ $Notebooks,
			PrintTemporary@ Labeled[
				ProgressIndicator[Appearance -> "Necklace"], msg, Right]
		,
			Print[msg]
		];
		URLSave[download, target]
	,
		Return[$Failed]
	];

	If[ FileExistsQ[target],
		(* we overwrite if the paclet is already installed because
			testing for a version is tricky, we anyways have downloaded
			the paclet and trying to install the same version gives
			an error message that is not wanted in this case. *)
		PacletManager`PacletInstall[target, "IgnoreVersion" -> True]
	,
		$Failed
	]
]


(**)
(* Boundary ElementMesh Booleans *)
(**)

ClearAll[validInputQ]
validInputQ[bm1_, bm2_] := 
	BoundaryElementMeshQ[bm1] && BoundaryElementMeshQ[bm2] &&
		(bm1["EmbeddingDimension"] === bm2["EmbeddingDimension"]) &&
			(bm1["MeshOrder"] === bm2["MeshOrder"] === 1)

(*
For RegionDifference, Intersection... (but not ReginUnion) we need a full region representation to do the operation. A boundary representation like MeshRegion does not work as the operation means something different. For RegionUnion it matters less but it generates sometimes wanted sometimes not wanted internal boundaries. So if a boundary mesh region can not be generated and the operation is RegionUnion then we assume that we have internal boundaries and do the operation. The disadvantage of using BoundaryMeshRegion is that internal boundaries can not be dealt with. Alternatively, we could generate a full element mesh from the boundary, convert that to a MeshRegion do the operation and convert that to a boundary element mesh; that however, crashes.
*)
ClearAll[MakeBoundaryElementMeshBooleanOperation]
MakeBoundaryElementMeshBooleanOperation[op_, funName_] := 
Return[
	funName[bm1_, bm2_, opts:OptionsPattern[ToBoundaryMesh]] /; validInputQ[bm1, bm2] :=
	Module[
		{m1, m2, rop, r, successQ},
		m1 = BoundaryMeshRegion[bm1];
		m2 = BoundaryMeshRegion[bm2];
		If[ (!BoundaryMeshRegionQ[m1] || !BoundaryMeshRegionQ[m2]) &&
				op === RegionUnion,
			m1 = MeshRegion[bm1];
			m2 = MeshRegion[bm2];
		];

		rop = op[m1, m2];
		successQ = (BoundaryMeshRegionQ[rop] || MeshRegionQ[rop]);

		WarnIf[ !successQ,
			funName,
			"The `1` of the boundary element meshes `2` and `3` failed",
			op, bm1, bm2
		];
		If[ !successQ, Return[$Failed, Module]];

		r = ToBoundaryMesh[rop, opts];
		r
	]
]


(**)
(* Difference *)
(**)
MakeBoundaryElementMeshBooleanOperation[RegionDifference,
	BoundaryElementMeshDifference];

(**)
(* Intersection *)
(**)
MakeBoundaryElementMeshBooleanOperation[RegionIntersection,
	BoundaryElementMeshIntersection];



(**)
(* Join *)
(**)

BoundaryElementMeshJoin[bm1_, bm2_, opts:OptionsPattern[ToBoundaryMesh]] /; validInputQ[bm1, bm2] :=
Module[
	{c1, c2, nc1, newBCEle, newPEle, eleTypes, markers},

	c1 = bm1["Coordinates"];
	c2 = bm2["Coordinates"];
	nc1 = Length[c1];

	newBCEle = bm2["BoundaryElements"];
	eleTypes = Head /@ newBCEle;
	If[ ElementMarkersQ[newBCEle],
		markers = ElementMarkers[newBCEle]
	,
		markers = Sequence[]
	];
	newBCEle = MapThread[#1[##2] &,
		{eleTypes, ElementIncidents[newBCEle] + nc1, markers}];

	newPEle = bm2["PointElements"];
	eleTypes = Head /@ newPEle;
	If[ ElementMarkersQ[newPEle],
		markers = ElementMarkers[newPEle]
	,
		markers = Sequence[]
	];
	newPEle = MapThread[#1[##2] &,
		{eleTypes, ElementIncidents[newPEle] + nc1, markers}];


	ToBoundaryMesh[
		"Coordinates" -> Join[c1, c2], 
		"BoundaryElements" -> Flatten[{bm1["BoundaryElements"], newBCEle}],
		"PointElements" -> Flatten[{bm1["PointElements"], newPEle}],
		opts
	]
]

(*BoundaryElementMeshJoin[bm1_,bmRest__]/;validInputQ[bm1,bmRest[[1]]]\
:=Fold[BoundaryElementMeshJoin,bm1,{bmRest}]*)

BoundaryElementMeshJoin[r1_, r2_, r3__] := 
 BoundaryElementMeshJoin[BoundaryElementMeshJoin[r1, r2], r3];


(**)
(* Rotation *)
(**)
BoundaryElementMeshRotation[bm_, rm_, opts:OptionsPattern[ToBoundaryMesh]] /; validInputQ[bm, bm] && MatrixQ[rm] := 
Module[
	{c},
	c = bm["Coordinates"];
	c = rm.# & /@ c;

	ToBoundaryMesh[
		"Coordinates" -> c,
		"BoundaryElements" -> bm["BoundaryElements"],
		"PointElements" -> bm["PointElements"],
		opts
	]
]


(**)
(* SymmetricDifference *)
(**)
MakeBoundaryElementMeshBooleanOperation[RegionSymmetricDifference,
	BoundaryElementMeshSymmetricDifference];


(**)
(* Translate *)
(**)
BoundaryElementMeshTranslate[bm_, vector_, opts:OptionsPattern[ToBoundaryMesh]] /; validInputQ[bm, bm] &&
	ArrayQ[vector, 1, NumericQ] && Length[vector] === bm["EmbeddingDimension"] := 
Module[
	{c, tt},
	c = bm["Coordinates"];
	tt = TranslationTransform[vector];
	c = tt /@ c;

	ToBoundaryMesh[
		"Coordinates" -> c, 
		"BoundaryElements" -> bm["BoundaryElements"], 
		"PointElements" -> bm["PointElements"],
		opts
	]
]



(**)
(* Union *)
(**)
MakeBoundaryElementMeshBooleanOperation[RegionUnion, BoundaryElementMeshUnion];

BoundaryElementMeshUnion[r1_, r2_, r3__] :=
 BoundaryElementMeshUnion[BoundaryElementMeshUnion[r1, r2], r3]



End[];

EndPackage[];

