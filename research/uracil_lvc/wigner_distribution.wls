#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["global`*"]


\[Omega]={734,771,1193,1383,1406,1673,1761,1794}/219474.63068;


T=273.15;m=1;\[Beta]=1/T/3.1668114*^-6


(* ::Subsection:: *)
(*Standard Deviations of Position and Momentum*)


q=Sqrt[0.5*Tanh[\[Beta]*\[Omega]/2]/\[Omega]/m];q//MatrixForm


p=Sqrt[0.5*Tanh[\[Beta]*\[Omega]/2]*\[Omega]*m];p//MatrixForm
