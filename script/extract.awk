#!/bin/awk -f

{
    # set the flag to 1 if the corresponding block is found
    if ($0 ~ "STATE 1 ACF") acf = 1
    if ($0 ~ "STATE POPULATION") population = 1
    if ($0 ~ "POTENTIAL POINTS") potential = 1
    if ($0 ~ "QND WFN PROPAGATION") qndpropag = 1
    if ($0 ~ "STATE 1 SPECTRUM") spectrum = 1

    # reset the flag if the block is finished
    if ($0 ~ "^$") {
        acf = 0; population = 0; potential = 0; qndpropag = 0; spectrum = 0
    }

    # print the block if the flag is set
    if (acf > 2) print > "acf1.mat"
    if (population > 2) print > "population.mat"
    if (potential > 2) print > "potential.mat"
    if (qndpropag > 2) printf "%9.4f %6.4f %6.4f %6.4f\n", $2, $6, $7, $8 > "rho.mat"
    if (spectrum > 2) print > "spectrum1.mat"

    # increment the flag for each line in the block
    if (acf > 0) acf++
    if (population > 0) population++
    if (potential > 0) potential++
    if (qndpropag > 0) qndpropag++
    if (spectrum > 0) spectrum++
}
