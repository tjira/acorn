$0 ~ "POTENTIAL POINTS" {flag=1; next} $0 ~ "^$" {flag=0} flag && (NF == 2 || NF == 3) {print}
