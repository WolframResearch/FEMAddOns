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


Begin["`Private`"];


(* ::Subsubsection:: *)
(*ElementMeshSmoothing*)


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


(* ::Subsubsection:: *)
(*Structured mesh*)


StructuredMesh::array = "Raster of input points must be a full array of numbers with depth of `1`.";
Options[StructuredMesh] = {InterpolationOrder->1};


getElementConnectivity[nx_, ny_] := Flatten[
	Table[{
		i+(j-1)(nx+1),
		i+(j-1)(nx+1)+1,
		i+j(nx+1)+1,
		i+j(nx+1)
		},
		{j, 1, ny},
        {i, 1, nx}
   ],
   1
]

getElementConnectivity[nx_, ny_, nz_] := Flatten[
	Table[{
		i+(j-1)(nx+1)+(k-1)(nx+1)(ny+1),
		i+(j-1)(nx+1)+(k-1)(nx+1)(ny+1)+1,
		i+j(nx+1)+(k-1)(nx+1)(ny+1)+1,
		i+j(nx+1)+(k-1)(nx+1)(ny+1),
		i+(j-1)(nx+1)+k(nx+1)(ny+1),
		i+(j-1)(nx+1)+k(nx+1)(ny+1)+1,
		i+j(nx+1)+k(nx+1)(ny+1)+1,
		i+j(nx+1)+k(nx+1)(ny+1)
        },
        {k, 1, nz},
        {j, 1, ny},
        {i, 1, nx}
    ],
    2
]


StructuredMesh[raster_, {nx_, ny_}, opts:OptionsPattern[]
	] /; Length[raster] > 1 :=
Module[
    {order, dim, restructured, xInt, yInt, zInt, nodes, connectivity},
    If[ Not @ ArrayQ[raster, 3, NumericQ],
		Message[StructuredMesh::array, 3+1];
		Return[$Failed, Module]
	];

    order = OptionValue[InterpolationOrder] /. Automatic->1;
    dim = Last @ Dimensions[raster];

    restructured = Transpose[raster, {3, 2, 1}];
    xInt = ListInterpolation[restructured[[1]], {{0, 1}, {0, 1}},
		InterpolationOrder->order];
    yInt = ListInterpolation[restructured[[2]], {{0, 1}, {0, 1}},
		InterpolationOrder->order];

    nodes = Flatten[#, 1]& @ If[ dim == 3,
        zInt = ListInterpolation[restructured[[3]], {{0, 1}, {0, 1}},
			InterpolationOrder->order];
        Table[{xInt[i, j], yInt[i, j], zInt[i, j]},
			{j, 0, 1, 1./ny}, {i, 0, 1, 1./nx}]
        , (* else *)
        Table[{xInt[i, j], yInt[i, j]},
			{j, 0, 1, 1./ny}, {i, 0, 1, 1./nx}]
    ];

    connectivity = getElementConnectivity[nx, ny];

    If[ dim == 3,
        ToBoundaryMesh["Coordinates"->nodes,
			"BoundaryElements"->{QuadElement[connectivity]}]
		, (* else *)
        ToElementMesh["Coordinates"->nodes,
			"MeshElements"->{QuadElement[connectivity]}]
    ]
]

StructuredMesh[raster_, {nx_, ny_, nz_}, opts:OptionsPattern[]
	] /; Length[raster] > 1 :=
Module[
    {order, restructured, xInt, yInt, zInt, nodes, connectivity},
    If[ Not @ ArrayQ[raster, 4, NumericQ],
		Message[StructuredMesh::array, 4+1];
		Return[$Failed, Module]
	];

    order = OptionValue[InterpolationOrder] /. Automatic->1;
       
    restructured = Transpose[raster, {4, 3, 2, 1}];
    xInt = ListInterpolation[restructured[[1]], {{0, 1}, {0, 1}, {0, 1}},
		InterpolationOrder->order];
    yInt = ListInterpolation[restructured[[2]], {{0, 1}, {0, 1}, {0, 1}},
		InterpolationOrder->order];
    zInt = ListInterpolation[restructured[[3]], {{0, 1}, {0, 1}, {0, 1}},
		InterpolationOrder->order];
    
    nodes = Flatten[
       Table[
          {xInt[i, j, k], yInt[i, j, k], zInt[i, j, k]},
          {k, 0, 1, 1./nz}, {j, 0, 1, 1./ny}, {i, 0, 1, 1./nx}
       ],
       2
    ];

    connectivity = getElementConnectivity[nx, ny, nz];
    
    ToElementMesh["Coordinates"->nodes,
		"MeshElements"->{HexahedronElement[connectivity]}]
]


(* ::Subsubsection:: *)
(*ToQuadMesh*)


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


(* ::Subsubsection:: *)
(*End Package*)


End[];

EndPackage[];

