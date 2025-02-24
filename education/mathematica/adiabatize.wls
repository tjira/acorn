#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"];


PlotPotentials[Ud_,xmin_,xmax_]:=GraphicsRow[{Plot[Evaluate@Diagonal@Ud,{x,xmin,xmax}],Plot[Evaluate@Sort@Eigenvalues@Ud,{x,xmin,xmax}]}]


(* ::Subsubsection:: *)
(*Harmonic*)


Subscript[Us,hm]:={
{0.5*x^2}
};Subscript[Us,hm]//MatrixForm


PlotPotentials[Subscript[Us,hm],-8,8]


(* ::Subsubsection:: *)
(*Double Well*)


Subscript[Us,dw]:={
{0.05*(0.6*x^2-16)^2}
};Subscript[Us,dw]//MatrixForm


PlotPotentials[Subscript[Us,dw],-8,8]


(* ::Subsubsection:: *)
(*Tully #1*)


Subscript[Ud,tully,1]:={
{Sign[x]*0.01*(1-Exp[-Sign[x]*1.6*x]),0.005*Exp[-x^2]},
{0.005*Exp[-x^2],-Sign[x]*0.01*(1-Exp[-Sign[x]*1.6*x])}
};Subscript[Ud,tully,1]//MatrixForm


PlotPotentials[Subscript[Ud,tully,1],-8,8]


(* ::Subsubsection:: *)
(*Tully #2*)


Subscript[Ud,tully,2]:={
{0,0.015*Exp[-0.06*x^2]},
{0.015*Exp[-0.06*x^2],-0.1*Exp[-0.28*x^2]+0.05}
};Subscript[Ud,tully,2]//MatrixForm


PlotPotentials[Subscript[Ud,tully,2],-8,8]


(* ::Subsubsection:: *)
(*Tully #3*)


Subscript[Ud,tully,3]:={
{6*^-4,If[x>0,0.1*(2-Exp[-0.9*x]),0.1*Exp[0.9*x]]},
{If[x>0,0.1*(2-Exp[-0.9*x]),0.1*Exp[0.9*x]],-6*^-4}
};Subscript[Ud,tully,3]//MatrixForm


PlotPotentials[Subscript[Ud,tully,3],-8,8]


(* ::Subsubsection:: *)
(*Double State #1*)


Subscript[Ud,2,1]:={
{0.001*x,0.001*Exp[-0.05*x^2]},
{0.001*Exp[-0.05*x^2],-0.001*x}
};Subscript[Ud,2,1]//MatrixForm


PlotPotentials[Subscript[Ud,2,1],-8,8]


(* ::Subsubsection:: *)
(*Double State #2*)


Subscript[Ud,2,2]:={
{0.01*Tanh[0.6*x],0.001*Exp[-x^2]},
{0.001*Exp[-x^2],-0.01*Tanh[0.6*x]}
};Subscript[Ud,2,2]//MatrixForm


PlotPotentials[Subscript[Ud,2,2],-8,8]


(* ::Subsubsection:: *)
(*Triple State #1*)


Subscript[Ud,3,3]:={
{0.001*x,0.001*Exp[-0.01*x^2],0},
{0.001*Exp[-0.01*x^2],0,0.001*Exp[-0.01*x^2]},
{0,0.001*Exp[-0.01*x^2],-0.001*x}
};Subscript[Ud,3,3]//MatrixForm


PlotPotentials[Subscript[Ud,3,3],-8,8]


(* ::Subsubsection:: *)
(*Triple State #2*)


Subscript[Ud,3,2]:={
{0.01*Tanh[0.5*x],0.001*Exp[-x^2],0},
{0.001*Exp[-x^2],0,0.001*Exp[-x^2]},
{0,0.001*Exp[-x^2],-0.01*Tanh[0.5*x]}
};Subscript[Ud,3,2]//MatrixForm


PlotPotentials[Subscript[Ud,3,2],-8,8]


(* ::Subsubsection:: *)
(*Triple State #3*)


Subscript[Ud,3,3]:={
{0.03*(Tanh[1.6*x]+Tanh[1.6*(x+7)]),0.005*Exp[-x^2],0.005*Exp[-(x+7)^2]},
{0.005*Exp[-x^2],-0.03*(Tanh[1.6*x]+Tanh[1.6*(x-7)]),0.005*Exp[-(x-7)^2]},
{0.005*Exp[-(x+7)^2],0.005*Exp[-(x-7)^2],-0.03*(Tanh[1.6*(x+7)]-Tanh[1.6*(x-7)])}
};Subscript[Ud,3,3]//MatrixForm


PlotPotentials[Subscript[Ud,3,3],-12,12]


(* ::Subsubsection:: *)
(*Quadruple State #1*)


Subscript[Ud,4,1]:={
{0.03*(Tanh[1.6*(x+7)]-Tanh[1.6*(x-7)])+0.03,0.005*Exp[-(x-7)^2],0.005*Exp[-(x+7)^2],0},
{0.005*Exp[-(x-7)^2],0.03*(Tanh[1.6*(x-7)]+Tanh[1.6*x]+Tanh[1.6*(x+7)]),0.005*Exp[-x^2],0.005*Exp[-(x+7)^2]},
{0.005*Exp[-(x+7)^2],0.005*Exp[-x^2],-0.03*(Tanh[1.6*x]+Tanh[1.6*(x-7)]+Tanh[1.6*(x+7)]),0.005*Exp[-(x-7)^2]},
{0,0.005*Exp[-(x+7)^2],0.005*Exp[-(x-7)^2],-0.03*(Tanh[1.6*(x+7)]-Tanh[1.6*(x-7)])-0.03}
};Subscript[Ud,4,1]//MatrixForm


PlotPotentials[Subscript[Ud,4,1],-12,12]


(* ::Subsubsection:: *)
(*Quadruple State #2*)


Subscript[Ud,4,2]:={
{0.01*Tanh[0.6*x],0.001*Exp[-x^2],0,0},
{0.001*Exp[-x^2],0.004*Tanh[0.6*x],0.001*Exp[-x^2],0},
{0,0.001*Exp[-x^2],-0.004*Tanh[0.6*x],0.001*Exp[-x^2]},
{0,0,0.001*Exp[-x^2],-0.01*Tanh[0.6*x]}
};Subscript[Ud,4,2]//MatrixForm


PlotPotentials[Subscript[Ud,4,2],-8,8]
