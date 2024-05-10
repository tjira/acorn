$0 ~ "\\|dD\\|     S1" {flagRHO=1; lineRHO=0} $0 ~ "^$" {flagRHO=0} ++lineRHO > 1 && flagRHO {printf "%9.4f %6.4f %6.4f %6.4f\n", $2, $6, $7, $8 > "rho.mat"}
$0 ~ "STATE 1 SPECTRUM" {flagSPECTRUM=1; lineSPECTRUM=0} $0 ~ "^$" {flagSPECTRUM=0} ++lineSPECTRUM > 2 && flagSPECTRUM {print > "spectrum1.mat"}
$0 ~ "POTENTIAL POINTS" {flagPOT=1; linePOT=0} $0 ~ "^$" {flagPOT=0} ++linePOT > 2 && flagPOT {print > "potential.mat"}
$0 ~ "STATE POPULATION" {flagPOP=1; linePOP=0} $0 ~ "^$" {flagPOP=0} ++linePOP > 2 && flagPOP {print > "pop.mat"}
$0 ~ "STATE 1 ACF" {flagACF=1; lineACF=0} $0 ~ "^$" {flagACF=0} ++lineACF > 2 && flagACF {print > "acf1.mat"}
