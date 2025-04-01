#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


(* ::Subsection:: *)
(*Assumptions*)


$Assumptions={a1\[Element]Reals,a1>0,a2\[Element]Reals,a2>0};


(* ::Subsection:: *)
(*Coordinates*)


(* ::Input:: *)
(*r={x,y,z};r1={x1,y1,z1};r2={x2,y2,z2};*)


(* ::Subsection:: *)
(*Gaussians and Norm*)


(* ::Input:: *)
(*G1=Exp[-a1*(r-r1) . (r-r1)]*)
(*G2=Exp[-a2*(r-r2) . (r-r2)]*)


Simplify[1/Sqrt@Integrate[Simplify[G1^2],{x,-Infinity,Infinity},{y,-Infinity,Infinity},{z,-Infinity,Infinity}]]


(* ::Subsection:: *)
(*Overlap*)


Simplify@Integrate[Simplify[G1*G2],{x,-Infinity,Infinity},{y,-Infinity,Infinity},{z,-Infinity,Infinity}]
