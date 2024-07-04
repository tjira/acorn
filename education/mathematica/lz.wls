#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"];


(* ::Subsubsection:: *)
(*Define the Model Variables and Hamiltonian*)


\[Alpha]=-7.59233*^-5;H12=0.000999428;start=-1*^6;end=1*^6;


H={{0.5*\[Alpha]*t,H12},{H12,-0.5*\[Alpha]*t}};ics={A[start]==1,B[start]==0};


(* ::Subsubsection:: *)
(*Exact Probability*)


Exp[-2*Pi*H12^2/-\[Alpha]]


(* ::Subsubsection:: *)
(*Define the Wavefunction and Schr\[ODoubleDot]dinger Equation*)


psi={A[t]*Exp@Integrate[-I*H[[1,1]],t],B[t]*Exp@Integrate[-I*H[[2,2]],t]};se=Simplify@Thread[I*D[psi,t]==H . psi];


(* ::Subsubsection:: *)
(*Solve the ODE and Calculate the Asymptotic Value*)


sol=DSolveValue[Flatten@{se,ics},{A[t],B[t]},{t,start,end}];Abs@#&/@sol^2/.t->end
