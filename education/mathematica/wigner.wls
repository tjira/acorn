#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


\[Psi][x_]=Exp[-(x-x0)^2]*Exp[+I*p0*x]


W[x_,p_]=Factor@Simplify@Integrate[Simplify[Refine[\[Psi][x+y]\[Conjugate],{p0\[Element]Reals,x\[Element]Reals,x0\[Element]Reals,y\[Element]Reals}]*\[Psi][x-y]]*Exp[2*I*p*y],{y,-Infinity,Infinity}]
