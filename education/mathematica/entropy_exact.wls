#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


kB=3.167*^-6;
m=2000;
T=100000;


beta=1/kB/T;
omega=1/Sqrt[m];
a=beta omega;


SkB=a/(Exp[a]-1)-Log[1-Exp[-a]]
