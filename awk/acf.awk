$0 ~ "STATE " i " ACF" {flag=1; next} $0 ~ "^$" {flag=0} flag && NF == 3 {print}
