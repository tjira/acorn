BEGIN {if (i == "") i = 1;}

$0 ~ "STATE " i " SPECTRUM" {flag=1; next} $0 ~ "^$" {flag=0} flag && NF == 2 {print}
