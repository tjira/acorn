# set the terminal and output name
set term eps; set output "output.eps"

# ignore the header
set key autotitle columnhead

# set the Y range
set yrange [-0.05:1.05]

# plot the files
plot for [i=1:words(paths)] for [c=2:*] word(paths, i) u 1:c w li not
