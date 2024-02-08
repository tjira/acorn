#!/usr/bin/env wolframscript

"define mass, force constant and other parameters"
m = 1; k = 1; range = {-4, 4}; states = 3; time = 20;

"kinetic and potential energy operator"
T[x_] := -Laplacian[#, {x}] / (2 * m) &
V[x_] := 0.5 * k * x^2 # &

"define hamiltonian operator"
H[x_] := T[x][#] + V[x][#] &

"solve the time independent schrodinger equation"
{eigv, eigf} = NDEigensystem[{H[x][psi[x]], DirichletCondition[psi[x] == 0, True]}, psi[x], {x, range[[1]], range[[2]]}, states,
    Method -> {"SpatialDiscretization" -> {"FiniteElement", "MeshOptions" -> {"MaxCellMeasure" -> 0.001}}}
]

"solve the time dependent schrodinger equation"
tdsol = psi[x, t] /. First@# & /@ NDSolve[
    {I * D[psi[x, t], {t, 1}] == H[x][psi[x, t]], psi[x, 0] == #, DirichletCondition[psi[x, t] == 0, True]}, psi[x, t], {x, range[[1]], range[[2]]}, {t, 0, time}
] & /@ eigf;

"get the limits of the time dependent plot"
ymin = Min@Table[{V[x][1], #[x, t] & /@ Flatten[{Function[{x, t}, Re@#], Function[{x, t}, Im@#]} & /@ (tdsol + eigv + I*eigv)]}, {x, range[[1]], range[[2]], 0.1}, {t, 0, time, 0.1}];
ymax = Max@Table[{#[x, t] & /@ Flatten[{Function[{x, t}, Re@#], Function[{x, t}, Im@#]} & /@ (tdsol + eigv + I*eigv)]}, {x, range[[1]], range[[2]], 0.1}, {t, 0, time, 0.1}];

"print state energies"
Print[Transpose@{Range[states] - 1, eigv}]

"save the plot of optimized states"
Export["adiabatic_quantum_dynamics-1.png", Plot[Evaluate@Flatten@{V[x][1], eigf + eigv}, {x, range[[1]], range[[2]]}, PlotRange -> All]]

"save the time dependent animation"
Export["adiabatic_quantum_dynamics-2.gif", Table[Plot[Evaluate@Flatten@{V[x][1], #[x, t] & /@ Flatten[{Function[{x, t}, Re@#], Function[{x, t}, Im@#]} & /@ (tdsol + eigv + I*eigv)]}, {x, range[[1]], range[[2]]},
    AxesOrigin -> {0, 0}, PlotRange -> {ymin - 0.1, ymax + 0.1}
], {t, 0, time, 0.1}]]
