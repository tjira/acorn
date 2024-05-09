$0 ~ "STATE POPULATION" {flag=1} $0 ~ "^$" {flag=0} flag && NF == 3 {printf("%9.4f %6.4f %6.4f\n"), $1, $2, $3}
