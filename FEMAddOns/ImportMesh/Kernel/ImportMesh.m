(* ::Package:: *)

(* ::Section:: *)
(*Begin package*)


(* ::Subsection:: *)
(*Header comments*)


(* :Title: ImportMesh *)
(* :Context: ImportMesh` *)
(* :Author: Matevz Pintar *)
(* :Summary: Utilities for importing FEM meshes from other software. *)
(* :Copyright: C3M d.o.o., 2019 *)

(* :Acknowledgements: b3m2a1, Oliver Ruebenkoenig *)


(* ::Subsection:: *)
(*Begin package*)


If[
	$VersionNumber<11,
	Print["Package requires Mathematica 11 or higher."];Abort[]
];

(* Mathematica FEM functionality (context "NDSolve`FEM`") is needed. *)
BeginPackage["ImportMesh`",{"NDSolve`FEM`"}];

If[
	$VersionNumber <= 12.2,
	stringToDouble = Internal`StringToDouble
,
	stringToDouble = Internal`StringToMReal
];

ClearAll[ImportMesh];


(* ::Subsection:: *)
(*Messages*)


ImportMesh::usage=(
	"ImportMesh[\"file\"] imports data from mesh file, returning a ElementMesh object.\n"<>
	"ImportMesh[\"string\", fmt] imports \"string\" in the specified format.\n"<>
	"ImportMesh[stream, fmt] imports the InputStream stream in the specified format."
);
(* All error/warning messages are currently attached to the only public symbol. *)
ImportMesh::nosup="Mesh file format \".`1`\" is currently not supported.";
ImportMesh::eltype="Element type `1` is not supported.";
ImportMesh::abaqus="Incremental node or element generation (*NGEN and *ELGEN keywords) is not supported.";
ImportMesh::fail="Failed to extract mesh from ``";
ImportMesh::gmshfrmt="GMSH file format version `1` is not supported.";


(* ::Subsection:: *)
(*Package Level Functions*)


(* The purpose of this public "`Package`" subcontext is that it makes implementation functions more available. 
 It won't be loaded by default into the $ContextPath, but it can then be loaded by Needs["ImportMesh`Package`"] 
 for those who are interested in doing more development work. *)

BeginPackage["`Package`"];

importAbaqusMesh::usage="";
importComsolMesh::usage="";
importGmshMesh::usage="";
importElfenMesh::usage="";
convertToElementMesh::usage="";
$importMeshFormats::usage="";

EndPackage[];


(* ::Section:: *)
(*Code*)


(* Begin private context *)
Begin["`Private`"];

(* 
Implementation for each mesh file format has its own private subcontext (e.g. ImportMesh`Private`Gmsh`).
This is because low level helper functions (e.g. getNodes) are doing same things differently for different formats.
Some common private functions are implemented in ImportMesh`Package` context and inside other subcontext they 
to be called by their full name.
*)


(* ::Subsection:: *)
(*The main public function*)


(* This is declared first so that all package-level functions can inherit from it. *)
$importMeshFormats=
	<|
		"inp"->
			<|
				"Name"->
					"Abaqus",
				"Function"->
					importAbaqusMesh,
				"Elements"->
					{}
				|>,
		"mes"->
			<|
				"Name"->
					"Elfen",
				"Function"->
					importElfenMesh,
				"Elements"->
					{}
				|>,
		"mphtxt"->
			<|
				"Name"->
					"Comsol",
				"Function"->
					importComsolMesh,
				"Elements"->
					{}
				|>,
		"msh"->
			<|
				"Name"->
					"Gmsh",
				"Function"->
					importGmshMesh,
				"Elements"->
					{}
				|>
		|>;


(* In the case that importing mesh file is not successful the function ImportMesh should return $Failed, 
just like the system function Import. *)

ImportMesh//Options={
	"ScaleSize"->1,
	"SpatialDimension"->Automatic,
	"ReturnElement"->"Mesh"
};

ImportMesh[file:_String|_File, opts:OptionsPattern[]]/;(
	FileExistsQ[file]||Message[ImportMesh::noopen,file]
	):=
		Module[
			{ext, fn, res},
			(*PrintTemporary["Converting mesh..."];*)
			ext=ToLowerCase[FileExtension[file]];
			fn=
				If[KeyExistsQ[$importMeshFormats, ext],
					$importMeshFormats[ext, "Function"],
					Message[ImportMesh::nosup, ext]
					];
			res=
				Catch[
					fn[
						file,
						FilterRules[{opts}, 
							Options[fn]
							]
						]
					];
			If[res===$Failed, Message[ImportMesh::fail, file]];
			res
		];


(* ::Subsection:: *)
(*Common functions*)


updateIncidents[{}, map_] := {};

updateIncidents[eleIn_, map_] := Module[
{ele = eleIn, eleInci, eleMarker, eleType},
  If[AssociationQ[map],
   (* the /. map is inefficient since it unpacks *)
   
   eleInci = Developer`ToPackedArray /@ (ElementIncidents[ele] /. map);
   eleMarker = ElementMarkers[ele];
   eleType = Head /@ ele;
   ele = MapThread[#1[#2, #3] &, {eleType, eleInci, eleMarker}];
   eleInci = Null;
   eleMarker = Null;
   eleType = Null;
   ];
  ele
];


makeUniqueAndCompleteMeshElementIncidents[nodesIn_, eleIn_] := 
 Module[
  {nodes, ele, map,
   eleInci, givenUniqueInciIDs, neededUniqueInciIDs
   },
  nodes = nodesIn;
  ele = eleIn;
  map = {};
  
  eleInci = ElementIncidents[ele];
  givenUniqueInciIDs = Union[Join @@ (Flatten /@ eleInci)];
  neededUniqueInciIDs = Range[Length[givenUniqueInciIDs]];
  If[givenUniqueInciIDs =!= neededUniqueInciIDs,
   map = AssociationThread[
     Rule[givenUniqueInciIDs, neededUniqueInciIDs]];
   ele = updateIncidents[ele, map];
   nodes = nodes[[givenUniqueInciIDs]];
   ];
  {nodes, ele, map}
  ];


(* This function is common for all file formats. It accepts nodes and  the list with all elements 
specified in a file and creates ElementMesh. *)
convertToElementMesh[nodesIn_,allElements_]:=Module[
	{sDim,point,line,surface,solid,meshElements,boundaryElements,nodes,map},
	
	nodes=nodesIn;

	(* Spatial dimensions of the problem. *)
	sDim=Last@Dimensions[nodes];

	If[ !(sDim === 3 || sDim === 2 || sDim === 1), Return[ $Failed, Module]];
	
	point=Cases[allElements,PointElement[__],2];
	line=Cases[allElements,LineElement[__],2];
	surface=Cases[allElements,TriangleElement[__]|QuadElement[__],2];
	solid=Cases[allElements,TetrahedronElement[__]|HexahedronElement[__],2];

	(* If number of dimensions is 3 but no solid elements are specified, 
	then we use ToBoundaryMesh to create ElementMesh. And similarly for 2 dimensions. *)
	If[
		TrueQ[(sDim==3&&solid=={})||(sDim==2&&surface=={})]
		,
		meshElements=Automatic
		,
		meshElements=Switch[sDim,1,line,2,surface,3,solid]
	];

	boundaryElements = Switch[sDim,1,point,2,line,3,surface];

	(* make sure element incidents are unique and complete *)
	If[ meshElements =!= Automatic,
		{nodes, meshElements, map} =
			 makeUniqueAndCompleteMeshElementIncidents[nodes, meshElements];
		boundaryElements = updateIncidents[boundaryElements, map];
		point = updateIncidents[point, map];

	];

	(* "Coordinates" and "MeshElements" are the only required fields. 
	Order of rules given is important: first nodes and elements, then options  *)	
	ToElementMesh[
		"Coordinates"->nodes,
		"MeshElements"->meshElements,
		"BoundaryElements"->boundaryElements,
		"PointElements"->point,
		"CheckIntersections"->False
	]
];


(* ::Subsection:: *)
(*Abaqus (.inp)*)


(* Begin private context *)
Begin["`Abaqus`"];

(*
Abaqus has quite involved specification of .inp input file. Majority of code in this subcontext is needed just to 
extract element topology (e.g. HexahedronElement) and its order from element type string (e.g. "TYPE=C3D20R"). 
Another problem is that element connectivity specification for elements with a lot of nodes sometimes breaks to a 
new line in a text file and this has to be somehow taken into account.
*)


(* ::Subsubsection:: *)
(*Process elements*)


(* This is a quick ugly hack to pass around spatial dimension that relies on this global symbol. *)
$spatialDimension=3;


(* This "type" argument it the next few functions is used only to issue a message with full element type string. *)
processLine[type_,string_]:=Which[
	StringStartsQ[string,"1"],{LineElement,2},
	StringStartsQ[string,"2"],{LineElement,3},
	True,Message[ImportMesh::eltype,type];Throw[$Failed]
];


processSurface[type_,string_]:=Which[
	StringStartsQ[string,"3"],{TriangleElement,3},
	StringStartsQ[string,"6"],{TriangleElement,6},
	StringStartsQ[string,"4"],{QuadElement,4},
	StringStartsQ[string,"8"],{QuadElement,8},
	True,Message[ImportMesh::eltype,type];Throw[$Failed]
];


processVolume[type_,string_]:=Which[
	StringStartsQ[string,"4"],{TetrahedronElement,4},
	StringStartsQ[string,"10"],{TetrahedronElement,10},
	StringStartsQ[string,"8"],{HexahedronElement,8},
	StringStartsQ[string,"20"],{HexahedronElement,20},
	True,Message[ImportMesh::eltype,type];Throw[$Failed]
];


processContinuumType[type_,inString_]:=Module[
	{string=inString},
	Which[
		StringStartsQ[string,"PS"|"PE"|"AX"],
		string=StringTrim[string,"PS"|"PE"|"AX"];
		$spatialDimension=2;
		processSurface[type,string]
		,
		StringStartsQ[string,"3D"],
		string=StringTrim[string,"3D"];processVolume[type,string]
		,
		True,Message[ImportMesh::eltype,type];Throw[$Failed]
	]
];


processElementType[type_String]:=Module[
	{string=type,keysCont,keysStruct},
	(* Order of keys is important because two of them start with the same letter. *)
	keysCont=Alternatives@@{"C","DC","Q"};
	keysStruct=Alternatives@@{"M3D","DS","STRI","S"};
	
	$spatialDimension=3;
	
	Which[
		(* Continuum elements can be 2D or 3D and we have to process this separately. *)
		StringStartsQ[string,keysCont],
		string=StringTrim[string,keysCont];processContinuumType[type,string]
		,
		StringStartsQ[string,keysStruct],
		string=StringTrim[string,keysStruct];processSurface[type,string]
		,
		StringStartsQ[string,"B2"],
		string=StringTrim[string,"B2"];$spatialDimension=2;processLine[type,string]
		,
		True,Message[ImportMesh::eltype,type];Throw[$Failed]
	]
];


processElements[type_String,flattenedConnectivity_,marker_Integer]:=Block[
	{head,noNodes,connectivity},
	{head,noNodes}=processElementType[type];
	connectivity=Partition[flattenedConnectivity,noNodes+1][[All,2;;]];
	(* $nodeNumbering is a global symbol here. This is a quick hack and could be solved better. *)
	head[
		connectivity/.$nodeNumbering,
		ConstantArray[marker,Length[connectivity]]
	]
];


(* ::Subsubsection:: *)
(*Helper functions*)


(* Helper function to get position of keyword in a list of lists. It assumes Abaqus keywords (e.g. *NODE)
are always at the beginning of line. *)
getPosition[list_List,key_]:=Position[list,key][[All,1]];


(* Collect lines until another symbol "*" appears. *)
takeLines[list_,start_Integer]:=TakeWhile[Drop[list,start],(StringPart[First@#,1]=!="*")&];


(* Collect all node specifications (can be scattered around file) and return their numbering and coordinates (X,Y,Z). *)
getNodes[list_]:=Module[
	{startLines,allNodeData,numbering,crds},
	startLines=getPosition[list,"*NODE"];
	allNodeData=Join@@(takeLines[list,#]&/@startLines);
	numbering=ToExpression[allNodeData[[All,1]]];
	crds=Map[
		stringToDouble,
		allNodeData[[All,2;;]],
		{2}
	];
	{numbering,crds}
];


getElementType[list_,pos_Integer]:=StringDelete[list[[pos,2]],"TYPE="];


getElementSet[list_,pos_Integer]:=If[
	Length[list[[pos]]]>2,
	StringDelete[list[[pos,3]],"ELSET="],
	Missing[]
];


getElementMarkers[list_, pos_Integer, lookup_] :=
Module[
	{default, materialName, marker, elIDs},
	default = {{},0};
	If[ Length[ list[[pos]]] < 2, Return[default, Module]];
	materialName = StringDelete[list[[pos, 2]], "ELSET="];
	marker = lookup[materialName];
	If[ MissingQ[marker],
		default	
	,
		elIDs = ToExpression[Join @@ takeLines[list, pos]];
		{elIDs, marker}
	]
]


getMarkerLookup[list_, solidSectionPos_] :=
Module[
	{matString, matName},

	matString = (list[[#]] & /@ solidSectionPos)[[All, {2, -1}]];
	matString = Select[(StringSplit[#, {"="}] & /@ matString),
		(#[[1, 1]] == "ELSET" && #[[2, 1]] == "MATERIAL") &];

	matName = matString[[All, 2, 2]];

	Association[ Thread[Rule[matString[[All, 1, 2]], matName /.
		MapThread[Rule, {#, Range[Length[#]]} &[Union[matName]]]]]]
]

getElementConnectivity[list_,pos_]:=ToExpression[Join@@takeLines[list,pos]];


getElements[list_]:=Module[
	{startLines,elSetStrings,types,flattenedConnectivity,markers},
	startLines=getPosition[list,"*ELEMENT"];

	elSetStrings=getElementSet[list,#]&/@startLines;

	types=getElementType[list,#]&/@startLines;	
	flattenedConnectivity=getElementConnectivity[list,#]&/@startLines;
	
	(* Values of variable ELSET is are sorted and assigned consecutive integer. If ELSET is not set, then
	elements get marker value 0. This is save in global symbol $markerNumbering for further use.
	This is a hack, and could be solved more elegantly. *)
	$markerNumbering=MapIndexed[#1->First[#2]&,Union@DeleteMissing[elSetStrings]];

	markers=elSetStrings/.Prepend[$markerNumbering,Missing[]->0];

	ele = MapThread[ processElements[#1,#2,#3]&, {types,flattenedConnectivity,markers} ];

	If[ Union[markers] === {0},
		elSetPos = getPosition[list,"*ELSET"];
		solidSectionPos = getPosition[list, "*SOLIDSECTION"];

		If[ (Length[solidSectionPos] > 0) && (Length[ elSetPos] > 0),
			markerLookUp = getMarkerLookup[list, solidSectionPos];
			elSetMarker = getElementMarkers[list, #, markerLookUp]& /@ elSetPos;
			elSetMarker = DeleteCases[ elSetMarker , {{}, 0}];
			numEle = Total[Length[ElementIncidents[#]]&/@ele];
			markers = ConstantArray[0, {numEle}];
			(markers[[#[[1]]]] = #[[2]])& /@ elSetMarker;
			markers = Internal`PartitionRagged[ markers,
				Length /@ ElementMarkers[ele]];
			ele = MapThread[#1[#2,#3]&,
				{Head /@ ele, ElementIncidents[ele], markers}];
		];

	];	

	ele
];


(* ::Subsubsection:: *)
(*Main function*)


Options[importAbaqusMesh]=Options[ImportMesh];

$importMeshFormats["inp", "Elements"]=
	{
		"MeshNodes",
		"MeshElements"
		};

importAbaqusMesh[list_List, ops:OptionsPattern[]]:=Module[
	{nodes,numbering,allElements,dim, ret=OptionValue["ReturnElement"]},

	(* Currently incremental node and element generation is not supported.*)
	If[getPosition[list,"*NGEN"|"*ELGEN"]=!={},Message[ImportMesh::abaqus];Throw[$Failed]];
	
	{numbering,nodes}=getNodes[list];
	
	$nodeNumbering=MapIndexed[#1->First[#2]&, numbering];
	
	allElements=getElements[list];
	If[ret==="MeshElements", Return[allElements]];
	
	(* Here we use the ugly hack. Value of global symbol $spatialDimension set at processing the element type
		is used to determine if we have 2D or 3D space mesh. *)
	dim=
		Replace[OptionValue["SpatialDimension"], 
			{
				Except[_?IntegerQ]:>$spatialDimension
				}
			];
	nodes=nodes[[All,1;;dim]];
	If[ret==="MeshNodes", Return[nodes]];
				
	convertToElementMesh[nodes,allElements]
];


importAbaqusMesh[file:_String?FileExistsQ|_InputStream, opts:OptionsPattern[]]:=
	importAbaqusMesh[
		DeleteCases[
			StringDelete[Whitespace]/@
				ToUpperCase@
					ReadList[
						file,
						Word,
						RecordLists->True,
						WordSeparators->{","},
						RecordSeparators -> {"\n"}
						],
			{s_String/;StringStartsQ[s,"**"]}
			],
		opts
		];


importAbaqusMesh[str_String, opts:OptionsPattern[]]:=
	importAbaqusMesh[
		StringToStream[str],
		opts
		];


End[]; (* "`Abaqus`" *)


(* ::Subsection:: *)
(*Comsol (.mphtxt)*)


(* Begin private context *)
Begin["`Comsol`"];


(* ::Subsubsection:: *)
(*Helper functions*)


getNumber[list_List,key_]:=ToExpression@Flatten@StringCases[list,x__~~key:>x];


getPosition[list_List,key_]:=ToExpression@Flatten@Position[list,key];


getNodes[list_]:=Module[
	{noNodes,start,listOfStrings},
	noNodes=getNumber[list," # number of mesh points"|" # number of mesh vertices"];
	start=getPosition[list,"# Mesh point coordinates"|"# Mesh vertex coordinates"];
	listOfStrings=Join@@MapThread[
		Take[list,{#1+1,#1+#2}]&,
		{start,noNodes}
	];
	Map[
		stringToDouble,
		StringSplit@listOfStrings,
		{2}
	]
];


takeLines[list_List,start_Integer,noLines_Integer]:=ToExpression@StringSplit@Take[list,{start+1,start+noLines}];


switchType={
	"vtx"->PointElement,
	"edg"|"edg2"->LineElement,
	"tri"|"tri2"->TriangleElement,"quad"|"quad2"->QuadElement,
	"tet"|"tet2"->TetrahedronElement,"hex"|"hex2"->HexahedronElement
};


modification["tri2",nodes_]:=nodes[[{1,2,3,4,6,5}]];
modification["quad",nodes_]:=nodes[[{1,2,4,3}]];
modification["quad2",nodes_]:=nodes[[{1,2,4,3,5,8,9,6}]];
modification["hex",nodes_]:=nodes[[{1,2,4,3,5,6,8,7}]];
modification["hex2",nodes_]:=nodes[[{1,2,4,3,5,6,8,7,9,12,13,10,23,26,27,24,14,16,22,20}]];
modification["tet2",nodes_]:=nodes[[{1,2,3,4,5,7,6,9,10,8}]];
modification[type_,nodes_]:=nodes;


getElements[list_,type_,length_,startElement_,startDomain_]:=With[
	{head=type/.switchType},
	{
	If[ startDomain === 0,
		head[
			(* +1 because node counting starts from 0 *)
			(modification[type,#]&/@takeLines[list,startElement,length])+1
		]
	,
		head[
			(* +1 because node counting starts from 0 *)
			(modification[type,#]&/@takeLines[list,startElement,length])+1,
			Flatten@takeLines[list,startDomain,length]
		]
	]
	}
];


(* ::Subsubsection:: *)
(*Main function*)


Options[importComsolMesh]=Options[ImportMesh];

$importMeshFormats["mphtxt", "Elements"]=
	{
		"MeshElements",
		"MeshNodes",
		"MeshTypes",
		"MeshStartElements",
		"MeshMarkers"
		};
importComsolMesh[list:{__String}, opts:OptionsPattern[]]:=Module[
	{nodes,types,lengths,startElements,startMarkers,
		allElements,
		ret=OptionValue["ReturnElement"]
		},
	
	types=Flatten@StringCases[list,Whitespace~~x__~~" # type name":>x];
	If[ret==="MeshTypes", Return[types]];
	lengths=getNumber[list," # number of elements"];
	startElements=getPosition[list,"# Elements"];
	If[ret==="MeshStartElements", Return[startElements]];
	(* I think both "Domains"  and "Geometric entity indices" can be considered as markers. *)
	startMarkers=getPosition[list,"# Domains"|"# Geometric entity indices"];
	If[ret==="MeshMarkers", Return[startMarkers]];
	If[startMarkers === {}, startMarkers = ConstantArray[0, {Length[types]}]];
	
	nodes=getNodes[list];
	If[ret==="MeshNodes", Return[nodes]];
	
	allElements=MapThread[
		getElements[list,#1,#2,#3,#4]&,
		{types,lengths,startElements,startMarkers}
	];
	If[ret==="MeshElements", Return[allElements]];
	
	convertToElementMesh[nodes,allElements]
];


importComsolMesh[file:_String?FileExistsQ|_InputStream, opts:OptionsPattern[]]:=
	importComsolMesh[ReadList[file,String], opts];
importComsolMesh[str_String, opts:OptionsPattern[]]:=
	importComsolMesh[StringToStream[str], opts];


End[]; (* "`Comsol`" *)


(* ::Subsection:: *)
(*Gmsh (.msh)*)


(* Begin private context *)
Begin["`Gmsh`"];


(* ::Subsubsection:: *)
(*Helper functions*)


getStartPosition[list_List,key_]:=ToExpression@First[Flatten@Position[list,key],$Failed];


getNumber[list_List,key_]:=ToExpression@Part[list,getStartPosition[list,key]+1];


getNodes[list_]:=Module[
	{noNodes,start},
	start=getStartPosition[list,"$Nodes"];
	noNodes=ToExpression@Part[list,start+1];
	Map[
		stringToDouble,
		Rest/@StringSplit@Take[list,{start+2,start+1+noNodes}],
		{2}
	]
];


getMarkers[list_]:=Module[
	{start,noEntites},
	start=getStartPosition[list,"$PhysicalNames"];
	(* In case there are no markers specified return $Failed *)
	If[start===$Failed,Return[$Failed]];
	start=start+1;
	noEntites=ToExpression@Part[list,start];
	(* {dimension, integer name, string name} *)
	ToExpression@StringSplit@Take[list,{start+1,start+noEntites}]
];


(* {Head, order, noNodes, reordering} *)
elementData={
	1->{LineElement,1,2,{1,2}},
	2->{TriangleElement,1,3,{1,2,3}},
	3->{QuadElement,1,4,{1,2,3,4}},
	4->{TetrahedronElement,1,4,{1,2,3,4}},
	5->{HexahedronElement,1,8,{1,2,3,4,5,6,7,8}},
	8->{LineElement,2,3,{1,2,3}},
	9->{TriangleElement,2,6,{1,2,3,4,5,6}},
	11->{TetrahedronElement,2,10,{1,4,2,3,8,10,5,9,6,7}},
	15->{PointElement,1,1,{1}},
	16->{QuadElement,2,8,{1,2,3,4,5,6,7,8}},
	17->{HexahedronElement,2,20,{1,2,3,4,5,6,7,8,9,12,17,10,18,11,19,20,13,16,14,15}},
	_->$Failed
};


restructure[list_]:=Block[
	{head,order,noNodes,nodes,noMarkers,markers,reordering},
	{head,order,noNodes,reordering}=list[[1,1]]/.elementData;
	noMarkers=list[[1,2]];
	nodes=Part[list,All,(3+noMarkers);;];
	markers=Part[list,All,3];
	head[nodes[[All,reordering]],markers]
];


getElements[list_]:=Module[
	{start,noEntites,raw},
	start=getStartPosition[list,"$Elements"]+1;
	noEntites=ToExpression@Part[list,start];
	raw=Rest/@ToExpression@StringSplit@Take[list,{start+1,start+noEntites}];
	
	restructure/@(Values@GroupBy[raw,First])
];


(* ::Subsubsection:: *)
(*Main function*)


$importMeshFormats["msh", "Elements"]={"MeshElements", "MeshNodes", "MeshMarkers"};

importGmshMesh//Options=Options[ImportMesh];		

importGmshMesh[list_List, opts:OptionsPattern[]]:=Module[
	{version,nodes,markers,allElements,ret},
	
	version=First@StringSplit@Part[list,getStartPosition[list,"$MeshFormat"]+1];
	If[
		Not@MemberQ[{"2.2"},version],
		Message[ImportMesh::gmshfrmt,version];Return[$Failed]
	];
	
	ret=OptionValue["ReturnElement"];
	
	nodes=getNodes[list];
	If[ret==="MeshNodes", Return[nodes]];
	markers=getMarkers[list];
	If[ret==="MeshMarkers", Return[markers]];
	allElements=getElements[list];
	If[ret==="MeshElements", Return[allElements]];
	
	convertToElementMesh[nodes,allElements]
];


importGmshMesh[file:_String?FileExistsQ|_InputStream, opts:OptionsPattern[]]:= importGmshMesh[ReadList[file,String], opts];

importGmshMesh[str_String, opts:OptionsPattern[]]:= importGmshMesh[StringToStream[str], opts];


End[]; (* "`Gmsh`" *)


(* ::Subsection:: *)
(*Elfen (.mes)*)


(* Begin private context *)
Begin["`Elfen`"];


(* ::Subsubsection:: *)
(*Helper functions*)


getPosition[list_List,key_]:=ToExpression@Flatten@Position[list,key];


getTypes[list_]:=Module[
	{start,n},
	start=First@getPosition[list,"element_type_numbers"];
	n=ToExpression@list[[start+1]];
	
	ToExpression/@list[[start+2;;start+1+n]]
];


getNodes[list_]:=Module[
	{start,dim,n},
	(* We assume all nodes are listed only on one place in the file. *)
	start=First@getPosition[list,"coordinates"];
	dim=ToExpression@list[[start+1]];
	n=ToExpression@list[[start+2]];
	
	Partition[
		stringToDouble/@list[[start+3;;start+3+(n*dim)-1]],
		dim
	]
];


$elfenTypes={
	21->{TriangleElement,{1,2,3}},
	22->{QuadElement,{1,2,3,4}},
	23->{TriangleElement,{1,3,5,2,4,6}},
	24|25->{QuadElement,{1,3,5,7,2,4,6,8}},
	31->{TetrahedronElement,{1,2,3,4}},
	33->{HexahedronElement,{1,2,3,4,5,6,7,8}},
	34->{TetrahedronElement,{1,3,5,10,2,4,6,8,9,7}},
	36->{HexahedronElement,{1,3,5,7,13,15,17,19,2,4,6,8,14,16,18,20,9,10,11,12}}
};


processElements[list_,startLine_Integer,type_Integer]:=Module[
	{nNodes,nElms,connectivity,head,reordering},
	
	nNodes=ToExpression@list[[startLine+1]];
	nElms=ToExpression@list[[startLine+2]];
	
	connectivity=Partition[
		ToExpression/@list[[startLine+3;;startLine+3+(nNodes*nElms)-1]],
		nNodes
	];
	{head,reordering}=type/.$elfenTypes;
	
	head[connectivity[[All,reordering]] ]
];


getElements[list_]:=Module[
	{startPositions,types},
	
	startPositions=getPosition[list,"element_topology"];
	types=getTypes[list];
	
	MapThread[
		processElements[list,#1,#2]&,
		{startPositions,types}
	]
];


(* ::Subsubsection:: *)
(*Main function*)


Options[importElfenMesh]=Options[ImportMesh];

$importMeshFormats["mes", "Elements"]=
	{
		"MeshElements",
		"MeshNodes"
		};
		
importElfenMesh[list_List, opts:OptionsPattern[]]:=Module[
	{nodes,allElements,ret=OptionValue["ReturnElement"]},
	
	nodes=getNodes[list];
	If[ret==="MeshNodes", Return[nodes]];
	allElements=getElements[list];
	If[ret==="MeshElements", Return[allElements]];
	
	convertToElementMesh[nodes,allElements]
];


importElfenMesh[file:_String?FileExistsQ|_InputStream, opts:OptionsPattern[]]:=
	importElfenMesh[
		ReadList[
			file,
			Word,
			RecordLists->True,
			RecordSeparators -> {"#","\"","{","}","*","\n"}
			]//Flatten, 
		opts
		];
importElfenMesh[str_String, opts:OptionsPattern[]]:=
	importElfenMesh[StringToStream[str], opts];


End[]; (* "`Elfen`" *)


(* ::Subsection:: *)
(*Register converters*)


(* Provide converters for the different mesh types defined. *)
If[!TrueQ@$importRegistered,
	KeyValueMap[
		Function[
			With[
				{
					func=#2["Function"], 
					names=ToUpperCase@{#, #2["Name"]<>"Mesh"}, 
					els=#2["Elements"]
					},
				Map[
					Function[
						ImportExport`RegisterImport[#,
							Join[
								Map[
									With[{el=#},
										el:>
											Function[{el->func[##, "ReturnElement"->el]}]
										]&,
									els
									],
								{
									"Mesh":>Function[{"Mesh"->func[##]}],
									"Elements":>Function[{"Elements"->Prepend[els, "Mesh"]}],
									func
									}
								],
							"FunctionChannels"->{"Streams"},
							"AvailableElements"->Prepend[els, "Mesh"]
							]
						],
					names
					]
				]
			],
		$importMeshFormats
		];
	ImportExport`RegisterImport["ElementMesh", ImportMesh]
	];
	
$importRegistered=True;


(* ::Section:: *)
(*End package*)


End[]; (* "`Private`" *)


EndPackage[];
