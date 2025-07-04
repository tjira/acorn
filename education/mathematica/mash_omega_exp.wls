#!/usr/bin/env wolframscript
(* ::Package:: *)

A={{0,-dV,2*T},{dV,0,0},{-2*T,0,0}};A//MatrixForm


FullSimplify@MatrixExp[A]//MatrixForm
