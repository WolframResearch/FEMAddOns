(*
   DistMesh - Wolfram Language Implementation
   Copyright (c) 2015 Wolfram Research, Inc.  All rights reserved.
   Author: Oliver Ruebenkoenig

   This software is based on the work of Per-Olof Persson and Gilbert
   Strang: http://persson.berkeley.edu/distmesh/

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*)


BeginPackage["DistMesh`", {"NDSolve`FEM`", "TriangleLink`", "TetGenLink`"}];

ClearAll[DistMesh];

DistMesh::usage = "DistMesh[r] creates an ElementMesh from a region r.";
DistMesh::"bcond" = "`1` can not be used as predicate.";
DistMesh::"bqal" = "The mesh quality goal of `1` could not be satisfied. The current mesh has a quality of `2`. Increasing the number of \"MaxIterations\" or lowering the \"MeshQualityGoal\" may help.";
DistMesh::"bgrad" = "A close to zero derivative has been found during the computation. Setting \"ScaleDerivative\" -> False may help.";
DistMesh::"noele" = "The requested \"MaxCellMeasure\" seems to coarse. Try setting \"MaxCellMeasure\" to a smaler value.";

Options[DistMesh] = Sort[Join[{
  MaxIterations -> 150,
  "MeshQualityFunction" -> Min,
  "ScaleDerivative" -> True,
  "DistMeshRefinementFunction" -> None},
  Options[ToElementMesh]]
];

Begin["`Private`"];

myDelaunay[pts_] :=
    Block[{tmp, dim, res},
      dim = Last[Dimensions[pts]];
      Switch[dim,
        2, res = TriangleDelaunay[pts],
        3, res = TetGenDelaunay[pts],
        True, res = {$Failed, $Failed}
      ]
    ];

(* this should be called lessThenZeroPosition *)
greaterEqualZeroPosition[th_] :=
    Flatten[SparseArray[1 - UnitStep[th]]["NonzeroPositions"]];

mkCompileCode[vars_, code_, idx_] :=
    With[{fVars = Flatten[vars]},
      Compile[{{in, _Real, idx}},
        Block[fVars,
          fVars = Flatten[in];
          code
        ],
        RuntimeAttributes -> Listable
      ]
    ];

ClearAll[mkcForces];
mkcForces[chCode_, dim_Integer] :=
    Compile[{{p, _Real, 2}, {bars, _Integer, 2}, {fScale, _Real, 0}},
      Block[{barVec, len, hBars, l0, force, forceVec},
        barVec = p[[ bars[[ All, 1]] ]] - p[[ bars[[ All, 2 ]] ]];
        len = Sqrt[Plus @@@ Power[barVec, 2]];
        hBars = chCode[((Plus @@ p[[#]]) & /@ bars) / 2];
        l0 = hBars * fScale * Power[(Plus @@ (len^dim)) / (Plus @@ (hBars^dim)), 1 / dim];
        force = Max[#, 0.] & /@ (l0 - len);
        forceVec = Transpose[{#, -#}] &[(force / len) * barVec];
        forceVec
      ]
    ];

Clear[moveMeshPoints];
moveMeshPoints[p_, bars_, cForces_, fScale_] :=
    Module[{totForce, forceVec},
      forceVec = cForces[p, bars, fScale];
      totForce = NDSolve`FEM`AssembleMatrix[
        {bars, ConstantArray[Range[Last[Dimensions[p]]], {Length[bars]}]},
        Flatten[forceVec, {{1}, {2, 3}}], Dimensions[p]
      ];
      totForce
    ];

initialMeshPoints[{{x1_, x2_}, {y1_, y2_}}, h0_] :=
    Module[{xRow, yColumn, tmp},
      xRow = Range[x1, x2, h0];
      yColumn = Range[ y1, y2, h0 * Sqrt[3] / 2 ];
      tmp = Transpose[{
        Flatten[Transpose[Table[If[OddQ[i], xRow, xRow + h0 / 2], {i, Length[yColumn ]}]]],
        Flatten[Table[yColumn, {Length[xRow]}]]
      }];
      Developer`ToPackedArray[tmp]
    ];

initialMeshPoints[{{x1_, x2_}, {y1_, y2_}, {z1_, z2_}}, h0_] :=
    Module[{xP, yP, zP},
      xP = Range[x1, x2, h0];
      yP = Range[y1, y2, h0];
      zP = Range[z1, z2, h0];
      Flatten[Transpose[Outer[List, xP, yP, zP], {3, 2, 1}], 2]
    ];

validRDFQ[a_] := NumericQ[a];
validRDFQ[___] := False;

DistMesh[em_?ElementMeshQ, opts : OptionsPattern[DistMesh]] :=
    Module[{nr},
      nr = ToNumericalRegion[em];
      DistMesh[ nr, opts]
    ];

DistMesh[sr_?RegionQ, opts : OptionsPattern[DistMesh]] :=
    Module[{nr},
      nr = ToNumericalRegion[sr];
      DistMesh[ nr, opts]
    ];

DistMesh[nr_NumericalRegion, opts : OptionsPattern[DistMesh]] :=
    Module[
      {
        maxIterations, iterationNumber, meshQualityFactor,
        meshQualityFunction, pfix, dim,
        p, ix, deps, cdfs, d, df, cond, vars, bbox,
        dgradDim, dgradN, depsIM, dgrad,
        pMid, t, bars, geps, rdf,
        crdf, cForces, mq, mqMax,
        barSelectors, saveBestMeshQ, mcm, mcmd, etype,
        dptol, ttol, fScale, deltat, r0, n, pOld, tmp, scaledR0, totForce,
        dgradx, dgrady, dgrad2, bmesh,
        pBest = $Failed, tBest = $Failed, maxf, mr, h0
      },

      maxIterations = OptionValue[MaxIterations];

      meshQualityFactor = OptionValue[MeshQualityGoal];
      If[! NumericQ[meshQualityFactor], meshQualityFactor = 0.5];
      meshQualityFactor = Min[1, Max[0, meshQualityFactor]];

      meshQualityFunction = OptionValue["MeshQualityFunction"];

      saveBestMeshQ = True;

      pfix = N[OptionValue["IncludePoints"]];
      If[! MatrixQ[pfix, NumericQ], pfix = {}];

      bbox = nr["Bounds"];
      dim = nr["EmbeddingDimension"];

      mcm = OptionValue["MaxCellMeasure"];
      If[etype === Automatic,
        Switch[dim,
          1, etype = NDSolve`FEM`LineElement,
          2, etype = NDSolve`FEM`TriangleElement,
          3, etype = NDSolve`FEM`TetrahedronElement
        ];
      ];
      mcmd = Region`Mesh`Utilities`ProcessMaxCellMeasure[
        bbox,
        MaxCellMeasure -> mcm,
        "ElementType" -> etype
      ];
      (* we always want the edge length *)
      h0 = mcmd[1];

      iterationNumber = 0;
      mqMax = 0;

      dptol = 0.001;
      ttol = 0.1;
      fScale = 1 + 0.4 / 2^(dim - 1);
      If[dim === 2,
        geps = 0.001 * h0;
        deltat = 0.2;
        ,
        geps = 0.1 * h0;
        deltat = 0.1;
      ];
      deps = Sqrt[$MachineEpsilon] * h0;
      cond = nr["Predicates"];
      vars = nr["PredicateVariables"];

      If[cond =!= Automatic,
      (* compile the distance set fuction *)

        df = Region`LevelFunction[cond];
        If[df === $Failed,
          Message[DistMesh::"bcond", cond];
          Return[$Failed];
        ];
        cdfs = mkCompileCode[vars, df, 1];
        ,
        bmesh = nr["BoundaryMesh"];
        If[!BoundaryElementMeshQ[bmesh],
          Message[DistMesh::"bcond", cond];
          Return[$Failed];
        ];
        mr = BoundaryMeshRegion[bmesh];
        If[BoundaryMeshRegionQ[mr],
          vars = Table[Unique["x", {Temporary}], {dim}];
          cdfs = SignedRegionDistance[mr];
          pfix = Join[pfix, bmesh["Coordinates"]];
        ];
      ];

      rdf = OptionValue["DistMeshRefinementFunction"];
      rdf = If[!validRDFQ[rdf @@ bbox[[-1]]],
        Function[1.],
        rdf
      ];

      (* compile relative mesh coarsness function *)

      crdf = mkCompileCode[vars, rdf @@ vars, 1];
      cForces = mkcForces[crdf, dim];

      barSelectors = Subsets[Range[dim + 1], {2}];

      p = initialMeshPoints[bbox, h0];
      p = p[[greaterEqualZeroPosition[cdfs[p] - geps]]];

      r0 = 1 / crdf[p]^dim;
      (* distmeshnd uses Min and a different scheme *)

      scaledR0 = r0 / Max[r0];

      p = p[[greaterEqualZeroPosition[RandomReal[{0, 1}, {Length[p]}] - scaledR0]]];
      If[pfix =!= {}, p = Sort[Complement[p, pfix]]];
      p = Join[pfix, p];
      n = Length[p];

      pOld = Infinity;
      While[True,
        iterationNumber++;

        If[Max[Sqrt[Total[(p - pOld)^2, {2}]] / h0] > ttol,
          pOld = p;
          If[ Length[p] < 3,
            Message[DistMesh::"noele"];
            Return[$Failed, Module]
          ];
          {p, t} = myDelaunay[p];
          pMid =
              Total[NDSolve`FEM`GetElementCoordinates[p, t], {2}] / (dim + 1);
          t = t[[greaterEqualZeroPosition[cdfs[pMid] + geps]]];
          bars = Union[Sort /@ Flatten[t[[ All, #]] & /@ barSelectors, 1]];
        ];

        totForce = moveMeshPoints[p, bars, cForces, fScale];
        totForce[[Range[Length[pfix]]]] = ConstantArray[0., {Length[pfix], dim}];
        p = p + deltat * totForce;

        d = cdfs[p];
        ix = greaterEqualZeroPosition[-(d + 0.)];

        dgradDim = Dimensions[p[[ix]]];
        dgradN = ConstantArray[0., dgradDim];
        depsIM = deps * IdentityMatrix[dim];

        Do[
          dgradN[[All, i]] = (cdfs[p[[ix]] +
              ConstantArray[depsIM[[i]], {Length[p[[ix]]]}]] - d[[ix]]) / deps,
          {i, Last[ dgradDim]}
        ];
        If[TrueQ[OptionValue["ScaleDerivative"]],
          dgrad = Total[dgradN^2, {2}];
          If[ Min[Abs[dgrad]] <= $MachineEpsilon,
            Message[DistMesh::"bgrad"];
            Break[]
		];
          (* dgrad could have zeros *)
          p[[ix]] = p[[ix]] - (d[[ix]] * dgradN) / dgrad;
          ,
          p[[ix]] = p[[ix]] - (d[[ix]] * dgradN);
        ];

        (* due to wrong ordering incident ordering may be negative *)
        mq = meshQualityFunction[NDSolve`FEM`MeshElementQuality[NDSolve`FEM`GetElementCoordinates[p, t]]];

        If[(mq >= meshQualityFactor),
          Break[]
        ];

        If[saveBestMeshQ && mq > mqMax,
          mqMax = mq;
          pBest = p;
          tBest = t
        ];

        maxf = Max[deltat * Sqrt[Total[totForce[[greaterEqualZeroPosition[d + geps]]]^2, {2}]] / h0];

        If[maxf < dptol * h0,
          Break[]
        ];

        If[iterationNumber >= maxIterations,
          If[ saveBestMeshQ && mq <= mqMax && tBest =!= $Failed,
            mq = mqMax;
            p = pBest;
            t = tBest
          ];
          Break[]
        ];
      ];

      If[ mq < meshQualityFactor,
        Message[DistMesh::"bqal", meshQualityFactor, mq]
      ];

      (*Print["it num: ", iterationNumber, "  mq: ", mq];*)

      ToElementMesh[
        "Coordinates" -> p,
        "MeshElements" -> {MeshElementType[dim, Last[Dimensions[t]], 1][t]},
        FilterRules[{opts}, Options[ToElementMesh]]
      ]
    ];

End[];
EndPackage[];

