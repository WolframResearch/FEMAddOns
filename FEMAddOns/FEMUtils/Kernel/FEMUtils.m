(* ::Package:: *)

(*
https://mathematica.stackexchange.com/questions/156611/improve-the-mesh-smoothing-procedure
*)


BeginPackage["FEMUtils`", {"NDSolve`FEM`"}];

ClearAll[ElementMeshSmoothing];
ElementMeshSmoothing::usage = "ElementMeshSmoothing[mesh] smoothes an ElementMesh.";
Options[ElementMeshSmoothing] = {Method -> Automatic};

ClearAll[StructuredMesh]
StructuredMesh::usage="StructuredMesh[raster,{nx,ny}] creates structured mesh of quadrilaterals.
StructuredMesh[raster,{nx,ny,nz}] creates structured mesh of hexahedra.";

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


getElementConnectivity[nx_,ny_]:=Flatten[
	Table[{
		i+(j-1)(nx+1),
		i+(j-1)(nx+1)+1,
		i+j(nx+1)+1,
		i+j(nx+1)
		},
		{j,1,ny},
        {i,1,nx}
   ],
   1
]

(* =====================================================================\[Equal] *)

getElementConnectivity[nx_,ny_,nz_]:=Flatten[
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
        {k,1,nz},
        {j,1,ny},
        {i,1,nx}
    ],
    2
]


StructuredMesh::array="Raster of input points must be full array of numbers with depth of `1`.";

StructuredMesh//Options={InterpolationOrder->1};
StructuredMesh[raster_,{nx_,ny_},opts:OptionsPattern[]]:=Module[
    {order,dim,restructured,xInt,yInt,zInt,nodes,connectivity},
    If[Not@ArrayQ[raster,3,NumericQ],Message[StructuredMesh::array,3];Return[$Failed]];

    order=OptionValue[InterpolationOrder]/.Automatic->1;
    dim=Last@Dimensions[raster];

    restructured=Transpose[raster,{3,2,1}];
    xInt=ListInterpolation[restructured[[1]],{{0,1},{0,1}},InterpolationOrder->order];
    yInt=ListInterpolation[restructured[[2]],{{0,1},{0,1}},InterpolationOrder->order];

    nodes=Flatten[#,1]&@If[dim==3,
        zInt=ListInterpolation[restructured[[3]],{{0,1},{0,1}},InterpolationOrder->order];
        Table[{xInt[i,j],yInt[i,j],zInt[i,j]},{j,0,1,1./ny},{i,0,1,1./nx}]
        ,
        Table[{xInt[i,j],yInt[i,j]},{j,0,1,1./ny},{i,0,1,1./nx}]
    ];

    connectivity=getElementConnectivity[nx,ny];

    If[dim==3,
        ToBoundaryMesh["Coordinates"->nodes,"BoundaryElements"->{QuadElement[connectivity]}],
        ToElementMesh["Coordinates"->nodes,"MeshElements"->{QuadElement[connectivity]}]
    ]
]

(* ===================================================================================== *)

StructuredMesh[raster_,{nx_,ny_,nz_},opts:OptionsPattern[]]:=Module[
    {order,restructured,xInt,yInt,zInt,nodes,connectivity},
    If[Not@ArrayQ[raster,4,NumericQ],Message[StructuredMesh::array,4];Return[$Failed]];

    order=OptionValue[InterpolationOrder]/.Automatic->1;
       
    restructured=Transpose[raster,{4, 3, 2, 1}];
    xInt=ListInterpolation[restructured[[1]],{{0,1},{0,1},{0,1}},InterpolationOrder->order];
    yInt=ListInterpolation[restructured[[2]],{{0,1},{0,1},{0,1}},InterpolationOrder->order];
    zInt=ListInterpolation[restructured[[3]],{{0,1},{0,1},{0,1}},InterpolationOrder->order];
    
    nodes=Flatten[
       Table[
          {xInt[i,j,k],yInt[i,j,k],zInt[i,j,k]},
          {k,0,1,1./nz},{j,0,1,1./ny},{i,0,1,1./nx}
       ],
       2
    ];

    connectivity=getElementConnectivity[nx,ny,nz];
    
    ToElementMesh["Coordinates"->nodes,"MeshElements"->{HexahedronElement[connectivity]}]
]


(* ::Subsubsection:: *)
(*End Package*)


End[];

EndPackage[];

