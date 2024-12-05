#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


\[Psi][x_]=Exp[-\[Gamma](x-x0)^2/2]*Exp[+I*p0*x]


W[x_,p_]=Factor@Simplify[Integrate[Simplify[Refine[\[Psi][x+y]\[Conjugate],{p0\[Element]Reals,x\[Element]Reals,x0\[Element]Reals,y\[Element]Reals,\[Gamma]\[Element]Reals}]*\[Psi][x-y]]*Exp[2*I*p*y],{y,-Infinity,Infinity}],Assumptions->{\[Gamma]\[Element]Reals,\[Gamma]>0}]


Subscript[\[Gamma], 0]=1/2;
Subscript[\[Sigma], p]=Sqrt[Subscript[\[Gamma], 0]/2]
Subscript[\[Sigma], q]=Sqrt[1/Subscript[\[Gamma], 0]/2]
