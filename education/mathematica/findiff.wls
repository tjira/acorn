#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


points={0,-1,-2};order=2;


A=Join[{ConstantArray[1,Length@points]},Table[points^k,{k,1,Length@points-1}]];A//MatrixForm
b=D[Table[x^k,{k,0,Length@points-1}],{x,order}]/.x->points[[1]];b//MatrixForm


LinearSolve[A,b]
