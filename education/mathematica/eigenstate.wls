#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["global`*"]


rmin=-8;rmax=8;


(* ::Subsubsection:: *)
(*1D Potential*)


\[ScriptCapitalH]=SchrodingerPDEComponent[{\[CapitalPsi][x],{x}},<|"ReducedPlanckConstant"->1,"SchrodingerPotential"->0.5*x^2,"Mass"->1|>];


{vals,funs}=NDEigensystem[\[ScriptCapitalH],\[CapitalPsi][x],{x,rmin,rmax},11,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{MaxCellMeasure->0.05}}}}];


Show[Plot[Evaluate[funs+vals],{x,rmin,rmax}],Plot[x^2,{x,rmin,rmax}],PlotRange->{{rmin,rmax},{0,11}}]


(* ::Subsubsection:: *)
(*2D Potential*)


\[ScriptCapitalH]=SchrodingerPDEComponent[{\[CapitalPsi][x,y],{x,y}},<|"ReducedPlanckConstant"->1,"SchrodingerPotential"->0.5*(x^2+y^2),"Mass"->1|>];


{vals,funs}=NDEigensystem[\[ScriptCapitalH],\[CapitalPsi][x,y],{x,rmin,rmax},{y,rmin,rmax},25,Method->{"PDEDiscretization"->{"FiniteElement",{"MeshOptions"->{MaxCellMeasure->0.05}}}}];


GraphicsGrid[ArrayReshape[{DensityPlot[Evaluate[Abs@#^2],{x,rmin,rmax},{y,rmin,rmax},PlotPoints->100,PlotRange->Full]&/@funs},{5,5}]]
