#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"];


PlotPotentials[Ud_,xmin_,xmax_]:=GraphicsRow[{Plot[Evaluate@Diagonal@Ud,{x,xmin,xmax}],Plot[Evaluate@Sort@Eigenvalues@Ud,{x,xmin,xmax}]}]


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


Subscript[Ud,3,1]:={
{0.001*(x-5),0.001*Exp[-0.01*x^2],0.002*Exp[-0.01*x^2]},
{0.001*Exp[-0.01*x^2],0,0.001*Exp[-0.01*x^2]},
{0.002*Exp[-0.01*x^2],0.001*Exp[-0.01*x^2],-0.001*(x+5)}
};Subscript[Ud,3,1]//MatrixForm


PlotPotentials[Subscript[Ud,3,1],-12,12]


(* ::Subsubsection:: *)
(*Triple State #2*)


Subscript[Ud,3,2]:={
{0.03*(Tanh[1.6*x]+Tanh[1.6*(x+7)]),0.005*Exp[-x^2],0.005*Exp[-(x+7)^2]},
{0.005*Exp[-x^2],-0.03*(Tanh[1.6*x]+Tanh[1.6*(x-7)]),0.005*Exp[-(x-7)^2]},
{0.005*Exp[-(x+7)^2],0.005*Exp[-(x-7)^2],-0.03*(Tanh[1.6*(x+7)]-Tanh[1.6*(x-7)])}
};Subscript[Ud,3,2]//MatrixForm


PlotPotentials[Subscript[Ud,3,2],-12,12]


(* ::Subsubsection:: *)
(*Triple State #3*)


Subscript[Ud,3,3]:={
{0.001*x,0.001*Exp[-0.01*x^2],0.002*Exp[-0.01*x^2]},
{0.001*Exp[-0.01*x^2],0,0.001*Exp[-0.01*x^2]},
{0.002*Exp[-0.01*x^2],0.001*Exp[-0.01*x^2],-0.001*x}
};Subscript[Ud,3,3]//MatrixForm


PlotPotentials[Subscript[Ud,3,3],-8,8]


(* ::Subsubsection:: *)
(*Triple State #4*)


Subscript[Ud,3,4]:={
{0.01*Tanh[0.6*x],0.002*Exp[-x^2],0.001*Exp[-x^2]},
{0.002*Exp[-x^2],0,0.002*Exp[-x^2]},
{0.001*Exp[-x^2],0.002*Exp[-x^2],-0.01*Tanh[0.6*x]}
};Subscript[Ud,3,4]//MatrixForm


PlotPotentials[Subscript[Ud,3,4],-8,8]


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
{0.01*Tanh[0.6*x],0.001*Exp[-x^2],0.001*Exp[-x^2],0.001*Exp[-x^2]},
{0.001*Exp[-x^2],0.004*Tanh[0.6*x],0.001*Exp[-x^2],0.001*Exp[-x^2]},
{0.001*Exp[-x^2],0.001*Exp[-x^2],-0.004*Tanh[0.6*x],0.001*Exp[-x^2]},
{0.001*Exp[-x^2],0.001*Exp[-x^2],0.001*Exp[-x^2],-0.01*Tanh[0.6*x]}
};Subscript[Ud,4,2]//MatrixForm


PlotPotentials[Subscript[Ud,4,2],-8,8]
