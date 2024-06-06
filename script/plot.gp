#!/bin/env -S gnuplot -c

# set the text variables
if (!exists("output")) {output="output"} # output file name
if (!exists("format")) {format="webp"} # output file format
if (!exists("alpha")) {alpha="00"} # transparency of the lines in hexadecimals
if (!exists("scale")) {scale="1"} # scaling factor for the data that should have the same number of elements as there are columns
if (!exists("shift")) {shift="0"} # shifting factor for the data that should have the same number of elements as there are columns
if (!exists("only")) {only=""} # chosen columns to plot from the interval [1, columns]
if (!exists("titles")) {titles=""} # titles for the plots

# fill the titles if not specified
if (titles eq "") do for [i=1:ARGC*columns] {titles = titles . i . " "}

# fill the columns to plot if not specified
if (only eq "") do for [i=1:columns] {only = only . i . " "}

# set the integer variables
if (!exists("height")) {height=480} # height of the output
if (!exists("columns")) {columns=1} # columns to plot each frame
if (!exists("frames")) {frames=1} # number of frames to plot
if (!exists("width")) {width=640} # width of the output

# set the terminal and output name
term = format; set output output . "." . format;

# fix the size of the output
if (format eq "eps") {height = height/160; width = width/160}
if (format eq "pdf") {height = height/160; width = width/160}

# set the format specific terminals
if (format eq "eps") {term = "epscairo"}
if (format eq "pdf") {term = "pdfcairo"}
if (format eq "png") {term = "pngcairo"}

# set the output type and format
if (frames == 1) {set term term lw 2 size width,height}
else {set term term animate lw 2 size width,height}

# ignore the header and emitted warnings
set key autotitle columnhead; unset warnings

# set the location of the legend
set key inside right top opaque spacing 1.1 outside

# set the default linetypes
set linetype  1 lc rgb "#" . alpha . "1f78b4"; set linetype  2 lc rgb "#" . alpha . "33a02c"
set linetype  3 lc rgb "#" . alpha . "e31a1c"; set linetype  4 lc rgb "#" . alpha . "ff7f00"
set linetype  5 lc rgb "#" . alpha . "6a3d9a"; set linetype  6 lc rgb "#" . alpha . "a6cee3"
set linetype  7 lc rgb "#" . alpha . "b2df8a"; set linetype  8 lc rgb "#" . alpha . "fb9a99"
set linetype  9 lc rgb "#" . alpha . "fdbf6f"; set linetype 10 lc rgb "#" . alpha . "cab2d6"

# set the domain and range
if (exists("xlabel")) set xlabel xlabel
if (exists("ylabel")) set ylabel ylabel
if (exists("xmin")) set xrange [xmin:]
if (exists("xmax")) set xrange [:xmax]
if (exists("ymin")) set yrange [ymin:]
if (exists("ymax")) set yrange [:ymax]

# plot the files
do for [i=2:columns*frames+1:columns] {
    plot for [j=1:ARGC] for [k=1:words(only)] ARGV[j] u 1:(word(scale, (word(only, k)-1)%words(scale)+1)*column(i+word(only, k)-1)+word(shift, (word(only, k)-1)%words(shift)+1)) w li s u t word(titles, (j-1)*words(only)+k)
}
