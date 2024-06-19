#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"];


Ud[x_]:={
{0.01*Tanh[0.6*x],0.001*Exp[-x^2],0.001*Exp[-x^2]},
{0.001*Exp[-x^2],0,0.001*Exp[-x^2]},
{0.001*Exp[-x^2],0.001*Exp[-x^2],-0.01*Tanh[0.6*x]}
};Ud[x]//MatrixForm


GraphicsRow[{Plot[Evaluate@Diagonal@Ud[x],{x,-8,8}],Plot[Evaluate@Sort@Eigenvalues[Ud[x]],{x,-8,8}]}]
