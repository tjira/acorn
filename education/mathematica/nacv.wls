#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


(* ::Subsection:: *)
(*Predefined Functions*)


SortedEigensystem[H_,I_]:=Module[{S=Eigensystem@H},Transpose@Table[{S[[1,i]],S[[2,i]]/Sqrt[S[[2,i]] . S[[2,i]]]},{i,I}]]


NonadiabaticCouplingMatrix[U_]:=Table[FullSimplify[U[[i]] . D[U[[j]],x]],{i,1,Length@U},{j,1,Length@U}]


(* ::Subsection:: *)
(*Symmetric Two-State Crossing*)


(* ::Subsubsection:: *)
(*Hamiltonian*)


H={{H11[x],H12[x]},{H12[x],-H11[x]}};H//MatrixForm


(* ::Subsubsection:: *)
(*Eigenvalues and Eigenvectors*)


{A,U}=SortedEigensystem[H,{1,2}];A//MatrixForm


(* ::Subsubsection:: *)
(*Nonadiabatic Coupling Matrix*)


If[#===0,0,Nonzero]&/@#&/@NonadiabaticCouplingMatrix@U//MatrixForm


(* ::Subsection:: *)
(*Symmetric Three-State Crossing*)


(* ::Subsubsection:: *)
(*Hamiltonian*)


H={{H11[x],H12[x],0},{H12[x],0,H12[x]},{0,H12[x],-H11[x]}};H//MatrixForm


(* ::Subsubsection:: *)
(*Eigenvalues and Eigenvectors*)


{A,U}=SortedEigensystem[H,{2,1,3}];A//MatrixForm


(* ::Subsubsection:: *)
(*Nonadiabatic Coupling Matrix*)


If[#===0,0,Nonzero]&/@#&/@NonadiabaticCouplingMatrix@U//MatrixForm
