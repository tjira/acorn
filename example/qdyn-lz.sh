#!/bin/bash

rm -rf result && mkdir result

TULLY_1="if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x))) 0.005*exp(-x^2)
         0.005*exp(-x^2)                                   (-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))"
TULLY_1_DIS=1; TULLY_1_AIS=1; TULLY_1_IX=-7; TULLY_1_IP=9

TULLY_2="0                    0.015*exp(-0.06*x^2)
         0.015*exp(-0.06*x^2) (-0.1)*exp(-0.28*x^2)+0.05"
TULLY_2_DIS=1; TULLY_2_AIS=1; TULLY_2_IX=-10; TULLY_2_IP=12

DS_1="0.01*tanh(0.6*x)  0.001*exp(-x^2)
      0.001*exp(-x^2) (-0.01)*tanh(0.6*x)"
DS_1_DIS=1; DS_1_AIS=1; DS_1_IX=-7; DS_1_IP=9

DS_2="0.001*x              0.001*exp(-0.05*x^2)
      0.001*exp(-0.05*x^2) (-0.001)*x"
DS_2_DIS=1; DS_2_AIS=1; DS_2_IX=-8; DS_2_IP=9

TS_1="0.03*(tanh(1.6*x)+tanh(1.6*(x+7))) 0.005*exp(-x^2)                       0.005*exp(-(x+7)^2)
      0.005*exp(-x^2)                    (-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7))) 0.005*exp(-(x-7)^2)
      0.005*exp(-(x+7)^2)                0.005*exp(-(x-7)^2)                   (-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))"
TS_1_DIS=1; TS_1_AIS=2; TS_1_IX=-10; TS_1_IP=15

TS_2="0.001*x         0.001*exp(-0.01*x^2) 0.002*exp(-0.01*x^2)
      0.001*exp(-0.01*x^2) 0                0.001*exp(-0.01*x^2)
      0.002*exp(-0.01*x^2) 0.001*exp(-0.01*x^2) (-0.001)*x"
TS_2_DIS=2; TS_2_AIS=2; TS_2_IX=-15; TS_2_IP=15

TS_3="0.01*tanh(0.6*x) 0.001*exp(-x^2) 0.001*exp(-x^2)
      0.002*exp(-x^2)  0               0.002*exp(-x^2)
      0.001*exp(-x^2)  0.002*exp(-x^2) (-0.01)*tanh(0.6*x)"
TS_3_DIS=2; TS_3_AIS=2; TS_3_IX=-10; TS_3_IP=15

QS_1="0.03*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))+0.03 0.005*exp(-(x-7)^2)                                0.005*exp(-(x+7)^2)                                   0
      0.005*exp(-(x-7)^2)                         0.03*(tanh(1.6*(x-7))+tanh(1.6*x)+tanh(1.6*(x+7))) 0.005*exp(-x^2)                                       0.005*exp(-(x+7)^2)
      0.005*exp(-(x+7)^2)                         0.005*exp(-x^2)                                    (-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7))+tanh(1.6*(x+7))) 0.005*exp(-(x-7)^2)
      0                                           0.005*exp(-(x+7)^2)                                0.005*exp(-(x-7)^2)                                   (-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))-0.03"
QS_1_DIS=2; QS_1_AIS=3; QS_1_IX=-10; QS_1_IP=15

QS_2="0.01*tanh(0.6*x) 0.001*exp(-x^2)   0.001*exp(-x^2)      0.001*exp(-x^2)
      0.001*exp(-x^2)  0.003*tanh(0.6*x) 0.001*exp(-x^2)      0.001*exp(-x^2)
      0.001*exp(-x^2)  0.001*exp(-x^2)   (-0.003)*tanh(0.6*x) 0.001*exp(-x^2)
      0.001*exp(-x^2)  0.001*exp(-x^2)   0.001*exp(-x^2)      (-0.01)*tanh(0.6*x)"
QS_2_DIS=3; QS_2_AIS=3; QS_2_IX=-10; QS_2_IP=15

A=0.05; C=tab:blue; I=3500; L=500; M=2000; R=1; T=100

for POT in TULLY_1 TULLY_2 DS_1 DS_2 TS_1 TS_2 TS_3 QS_1 QS_2; do
# for POT in QS_2; do
    P=${POT}_IP; X=${POT}_IX; DIS=${POT}_DIS; AIS=${POT}_AIS; STATES=$(echo "sqrt($(echo ${!POT} | wc -w))" | bc)

    PSI=$(for i in $(seq 1 $((2*${!DIS}))); do echo " 0"; done)" exp(-(x-(${!X}))^2) 0"$(for i in $(seq 1 $((2*($STATES-${!DIS}-1)))); do echo " 0"; done)

    acorn_expression -g -24 40 -o U_DIA.mat -p 1025 -e ${!POT}
    acorn_expression -g -24 40 -o PSI_DIA_GUESS.mat -p 1025 -e $PSI
    acorn_qdyn -i 350 -m $M -p ${!P} -s 10 --adiabatic --savewfn --align --factor 0.01

    for MODE in "" "--adiabatic"; do
        acorn_cdyn -c ${!X} -e $([ -z $MODE ] && echo ${!DIS} || echo ${!AIS}) -i $I -l $L -m $M -p ${!P} -r $R -s 1 -t $T -u ${!POT} $MODE --savetraj
    done

    plot.py PSI_ADIA_0.mat -c $(($STATES*2)) --png --title "Adiabatic Potential & Initial Wavefunction" && mv output.png U_PSI_ADIA_$POT.png
    plot.py PSI_DIA_0.mat -c $(($STATES*2)) --png --title "Adiabatic Potential & Initial Wavefunction" && mv output.png U_PSI_DIA_$POT.png
    plot.py P_ADIA_0.mat P_LZ-ADIA_ADIA.mat -c $(($STATES*$STATES)) -e $(seq 1 $(($STATES+1)) $(($STATES*$STATES))) --png --title "EXACT vs LZ-ADIA: Adiabatic Populations" && mv output.png P_EXACT_v_LZ-ADIA_ADIA_$POT.png
    plot.py P_DIA_0.mat P_LZ-ADIA_DIA.mat -c $(($STATES*$STATES)) -e $(seq 1 $(($STATES+1)) $(($STATES*$STATES))) --png --title "EXACT vs LZ-ADIA: Diabatic Populations" && mv output.png P_EXACT_v_LZ-ADIA_DIA_$POT.png
    plot.py P_ADIA_0.mat P_LZ-DIA_ADIA.mat -c $(($STATES*$STATES)) -e $(seq 1 $(($STATES+1)) $(($STATES*$STATES))) --png --title "EXACT vs LZ-DIA: Adiabatic Populations" && mv output.png P_EXACT_v_LZ-DIA_ADIA_$POT.png
    plot.py P_DIA_0.mat P_LZ-DIA_DIA.mat -c $(($STATES*$STATES)) -e $(seq 1 $(($STATES+1)) $(($STATES*$STATES))) --png --title "EXACT vs LZ-DIA: Diabatic Populations" && mv output.png P_EXACT_v_LZ-DIA_DIA_$POT.png
    plot.py P_LZ-ADIA_DIA.mat P_LZ-ADIA_ADIA.mat -c $(($STATES*$STATES)) -e $(seq 1 $(($STATES+1)) $(($STATES*$STATES))) --png --title "LZ-ADIA: Both Populations" && mv output.png P_LZ-ADIA_BOTH_$POT.png
    plot.py P_LZ-DIA_DIA.mat P_LZ-DIA_ADIA.mat -c $(($STATES*$STATES)) -e $(seq 1 $(($STATES+1)) $(($STATES*$STATES))) --png --title "LZ-DIA: Both Populations" && mv output.png P_LZ-DIA_BOTH_$POT.png
    montage -geometry "2400x1800" -title "Dynamics Summary (P=${!P}, TRAJS=$T)" U_PSI_DIA_$POT.png P_EXACT_v_LZ-DIA_DIA_$POT.png P_EXACT_v_LZ-DIA_ADIA_$POT.png P_LZ-DIA_BOTH_$POT.png U_PSI_ADIA_$POT.png P_EXACT_v_LZ-ADIA_DIA_$POT.png P_EXACT_v_LZ-ADIA_ADIA_$POT.png P_LZ-ADIA_BOTH_$POT.png 1_$POT.png

    # plot.py ETOT_LZ-ADIA_*.mat -p $C -a $A --png --title "Adiabatic Trajectory: Total Energy" && mv output.png ETOT_LZ-ADIA_$POT.png
    # plot.py ETOT_LZ-DIA_*.mat -p $C -a $A --png --title "Diabatic Trajectory: Total Energy" && mv output.png ETOT_LZ-DIA_$POT.png
    # plot.py R_LZ-ADIA_*.mat -p $C -a $A --png --title "Adiabatic Trajectory: Position" && mv output.png R_LZ-ADIA_$POT.png
    # plot.py R_LZ-DIA_*.mat -p $C -a $A --png --title "Diabatic Trajectory: Position" && mv output.png R_LZ-DIA_$POT.png
    # plot.py V_LZ-ADIA_*.mat -p $C -a $A --png --title "Adiabatic Trajectory: Velocity" && mv output.png V_LZ-ADIA_$POT.png
    # plot.py V_LZ-DIA_*.mat -p $C -a $A --png --title "Diabatic Trajectory: Velocity" && mv output.png V_LZ-DIA_$POT.png
    # plot.py A_LZ-ADIA_*.mat -p $C -a $A --png --title "Adiabatic Trajectory: Acceleration" && mv output.png A_LZ-ADIA_$POT.png
    # plot.py A_LZ-DIA_*.mat -p $C -a $A --png --title "Diabatic Trajectory: Acceleration" && mv output.png A_LZ-DIA_$POT.png
    # montage -geometry "2400x1800" ETOT_LZ-DIA_$POT.png R_LZ-DIA_$POT.png V_LZ-DIA_$POT.png A_LZ-DIA_$POT.png ETOT_LZ-ADIA_$POT.png R_LZ-ADIA_$POT.png V_LZ-ADIA_$POT.png A_LZ-ADIA_$POT.png 2_$POT.png

    # plot.py E_LZ-ADIA_*.mat -p $C -a $A --png --title "Adiabatic Simulation: TD Particle Energy" && mv output.png E_LZ-ADIA_$POT.png
    # plot.py E_LZ-DIA_*.mat -p $C -a $A --png --title "Diabatic Simulation: TD Particle Energy" && mv output.png E_LZ-DIA_$POT.png
    # plot.py ED_LZ-ADIA_*.mat -p $C -a $A -c 2 -e 1 --png --title "Adiabatic Simulation: TD Energy Differences" && mv output.png ED_LZ-ADIA_$POT.png
    # plot.py ED_LZ-DIA_*.mat -p $C -a $A -c 2 -e 1 --png --title "Diabatic Simulation: TD Energy Differences" && mv output.png ED_LZ-DIA_$POT.png
    # plot.py DED_LZ-ADIA_*.mat -p $C -a $A -c 2 -e 1 --png --title "Adiabatic Simulation: TD Energy Differences - 1st Derivative" && mv output.png DED_LZ-ADIA_$POT.png
    # plot.py DED_LZ-DIA_*.mat -p $C -a $A -c 2 -e 1 --png --title "Diabatic Simulation: TD Energy Differences - 1st Derivative" && mv output.png DED_LZ-DIA_$POT.png
    # plot.py DDED_LZ-ADIA_*.mat -p $C -a $A -c 2 -e 1 --png --title "Adiabatic Simulation: TD Energy Differences - 2nd Derivative" && mv output.png DDED_LZ-ADIA_$POT.png
    # plot.py DDED_LZ-DIA_*.mat -p $C -a $A -c 2 -e 1 --png --title "Diabatic Simulation: TD Energy Differences - 2nd Derivative" && mv output.png DDED_LZ-DIA_$POT.png
    # montage -geometry "2400x1800" E_LZ-DIA_$POT.png ED_LZ-DIA_$POT.png DED_LZ-DIA_$POT.png DDED_LZ-DIA_$POT.png E_LZ-ADIA_$POT.png ED_LZ-ADIA_$POT.png DED_LZ-ADIA_$POT.png DDED_LZ-ADIA_$POT.png 3_$POT.png

    convert -append [1-9]_$POT.png result/$POT.png && rm -rf *.mat
done

rm -rf *.png
