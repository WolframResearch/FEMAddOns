(*
https://mathematica.stackexchange.com/questions/156611/improve-the-mesh-smoothing-procedure
*)


BeginPackage["FEMUtils`", {"NDSolve`FEM`"}];

ClearAll[ElementMeshSmoothing];

ElementMeshSmoothing::usage = "ElementMeshSmoothing[mesh] smoothes an ElementMesh.";

Options[ElementMeshSmoothing] = {Method -> Automatic};

Begin["`Private`"];



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
	{n, vec, mat, adjacencymatrix2, mass2, laplacian2, 
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

	ToElementMesh["Coordinates" -> newCoords, 
		"MeshElements" -> mesh["MeshElements"], 
		"BoundaryElements" -> mesh["BoundaryElements"], 
		"PointElements" -> mesh["PointElements"], 
		"CheckIncidentsCompleteness" -> False, 
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

End[];

EndPackage[];

