(* ::Package:: *)

BeginPackage["FEMUtils`", {"NDSolve`FEM`"}];


ClearAll[ElementMeshSmoothing];
ElementMeshSmoothing::usage = "ElementMeshSmoothing[mesh] smoothes an ElementMesh.";
Options[ElementMeshSmoothing] = {Method -> Automatic};

ClearAll[StructuredMesh]
StructuredMesh::usage="StructuredMesh[raster,{nx,ny}] creates structured mesh of quadrilaterals.
StructuredMesh[raster,{nx,ny,nz}] creates structured mesh of hexahedra.";

ClearAll[TriangleToQuad]
TriangleToQuad::usage="TriangleToQuad[mesh] converts triangular mesh to quadrilateral mesh.";


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
    If[Not@ArrayQ[raster,3,NumericQ],Message[StructuredMesh::array,3+1];Return[$Failed]];

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
    If[Not@ArrayQ[raster,4,NumericQ],Message[StructuredMesh::array,4+1];Return[$Failed]];

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
(*TriangleToQuad*)


(*
Source of algorithm for mesh conversion is:
	Houman Borouchaki, Pascal J. Frey; Adaptive triangular\[Dash]quadrilateral mesh generation;
	International Journal for Numerical Methods in Engineering; 1998; Vol. 41, p. (915-934)

Implementation is taken from AceFEM package (http://symech.fgg.uni-lj.si/), 
written by prof. Joze Korelc.
 
*)

TriangleToQuad::elmtype="Only conversion of pure triangular meshes is supported."
TriangleToQuad::order="Only conversion of meshes with \"MeshOrder\"->1 is supported."

TriangleToQuad[mesh_ElementMesh]:=Module[
	{coor,econn,elem,distortion,dist,ncoor,taken,quad,triag,edge,nc,allquads,
	alltriag,region,nodes,regionout},
	
	If[mesh["MeshOrder"]=!=1,Message[TriangleToQuad::order];Return[$Failed]];
	
	coor=mesh["Coordinates"];
	econn=mesh["ElementConnectivity"][[1]];
	elem=Cases[mesh["MeshElements"],_TriangleElement];
	
	If[Length[elem]=!=1,
		Message[TriangleToQuad::elmtype];Return[$Failed],
		If[MatchQ[elem,{TriangleElement[_]}],
			elem=elem[[1,1]];
			region=ConstantArray[1,Length[elem]];
			regionout=False
			,
			region=elem[[1,2]];
			regionout=True;
			elem=elem[[1,1]]
		]
	];

	If[Length[econn]=!=Length[elem],
		Message[TriangleToQuad::elmtype];Return[$Failed]
	];

	distortion=Compile[{{n1,_Real,1},{n2,_Real,1},{n3,_Real,1},{n4,_Real,1}},
		2./\[Pi] * Max@Map[
			Abs[\[Pi]/2.-If[#<0,\[Pi]-#,#]]&@ArcTan[#[[1,1]]#[[2,1]]+#[[1,2]] #[[2,2]],#[[1,1]] #[[2,2]]-#[[1,2]] #[[2,1]]]&,
			{{n2-n1,n4-n1},{n3-n2,n1-n2},{n4-n3,n2-n3},{n1-n4,n3-n4}}
			]
	];

	dist=MapThread[ (
		{
		If[#2[[1]]==0,
			Nothing,
			distortion[Sequence@@coor[[nodes={#[[1]],#[[2]],Complement[elem[[#2[[1]]]],#1][[1]],#[[3]]}]]]->{#3,#2[[1]],nodes}
		],
		If[#2[[2]]==0,
			Nothing,
			distortion[Sequence@@coor[[nodes={#[[1]],#[[2]],#[[3]],Complement[elem[[#2[[2]]]],#1][[1]]}]]]->{#3,#2[[2]],nodes}
		],
		If[#2[[3]]==0,
			Nothing,
			distortion[Sequence@@coor[[nodes={#[[1]],Complement[elem[[#2[[3]]]],#1][[1]],#[[2]],#[[3]]}]]]->{#3,#2[[3]],nodes}
		]
		})&,
		{elem,econn,Range[elem//Length]}
	]//Flatten;

	taken=ConstantArray[False,elem//Length];
	quad=Map[
		If[
			Or@@taken[[#[[2,{1,2}]]]]  || #[[1]]>0.8 || region[[#[[2,1]]]] =!=region[[#[[2,2]]]],
			Nothing,
			taken[[#[[2,{1,2}]]]]=True;{#[[2,3]],region[[#[[2,1]]]]}
		]&,
		Sort[dist]
	]//Transpose;
	
	dist=.;
	triag={Extract[elem,#],Extract[region,#]}&@Position[taken,False];
	edge=DeleteDuplicates[
		Join[
			Flatten[Map[{#[[{1,2}]]//Sort,#[[{2,3}]]//Sort,#[[{3,4}]]//Sort,#[[{4,1}]]//Sort} &,quad[[1]]],1],
			Flatten[Map[{#[[{1,2}]]//Sort,#[[{2,3}]]//Sort,#[[{3,1}]]//Sort} &,triag[[1]]],1]
		]
	];
	
	ncoor=coor//Length;
	nc=Dispatch[Flatten[Map[(ncoor++;{#->ncoor,#[[{2,1}]]->ncoor})&,edge]]];
	
	allquads=MapThread[
		(ncoor++; 
		{Total[coor[[#]]]/4,
			{
			{#[[1]],#[[{1,2}]]/.nc,ncoor,#[[{4,1}]]/.nc},
			{#[[{1,2}]]/.nc,#[[2]],#[[{2,3}]]/.nc,ncoor},
			{ncoor,#[[{2,3}]]/.nc,#[[3]],#[[{3,4}]]/.nc},
			{#[[{4,1}]]/.nc,ncoor,#[[{3,4}]]/.nc,#[[4]]}
			},
		{#2,#2,#2,#2}
		})& ,
		quad
	];
	
	alltriag=MapThread[
		(ncoor++;
		{Total[coor[[#]]]/3,
			{
			{#[[1]],#[[{1,2}]]/.nc,ncoor,#[[{3,1}]]/.nc},
			{#[[{1,2}]]/.nc,#[[2]],#[[{2,3}]]/.nc,ncoor},
			{ncoor,#[[{2,3}]]/.nc,#[[3]],#[[{3,1}]]/.nc}
			},
		{#2,#2,#2}})& ,
		triag
	];

	ToElementMesh[
		"Coordinates"->Join[coor,Map[(coor[[#[[1]]]]+coor[[#[[2]]]])/2&,edge],allquads[[All,1]],alltriag[[All,1]]],
		"MeshElements"->If[regionout,
			{QuadElement[Flatten[Join[allquads[[All,2]],alltriag[[All,2]]],1],Flatten[{allquads[[All,3]],alltriag[[All,3]]}]]},
			{QuadElement[Flatten[Join[allquads[[All,2]],alltriag[[All,2]]],1]]}
		]
	]
]


(* ::Subsubsection:: *)
(*End Package*)


End[];

EndPackage[];

