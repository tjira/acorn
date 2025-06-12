#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


\[Psi][x_]=Exp[-\[Gamma](x-x0)^2/2]*Exp[I*p0*x];\[Gamma]=2;\[Psi][x]


W[x_,p_]=Factor@Simplify[Integrate[Simplify[Refine[\[Psi][x+y]\[Conjugate],{p0\[Element]Reals,x\[Element]Reals,x0\[Element]Reals,y\[Element]Reals,\[Gamma]\[Element]Reals}]*\[Psi][x-y]]*Exp[2*I*p*y],{y,-Infinity,Infinity}],Assumptions->{\[Gamma]\[Element]Reals,\[Gamma]>0}]


Solve[1/2/Subscript[\[Sigma],q]^2==\[Gamma],Subscript[\[Sigma],q],Assumptions->{Subscript[\[Sigma],q]>0,\[Gamma]>0}][[1]]
Solve[1/2/Subscript[\[Sigma],p]^2==1/\[Gamma],Subscript[\[Sigma],p],Assumptions->{Subscript[\[Sigma],p]>0,\[Gamma]>0}][[1]]



