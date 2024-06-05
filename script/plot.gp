#!/bin/env -S gnuplot -c

# set the text variables
if (!exists("background")) {background=""} # file to plot in the background
if (!exists("output")) {output="output"} # output file name
if (!exists("format")) {format="webp"} # output file format
if (!exists("alpha")) {alpha="00"} # transparency of the lines in hexadecimals
if (!exists("scale")) {scale="1"} # scaling factor for the data that should have the same number of elements as there are columns
if (!exists("shift")) {shift="0"} # shifting factor for the data that should have the same number of elements as there are columns

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

# set the default linetypes
set linetype  1 lc rgb "#" . alpha . "1f77b4"; set linetype  2 lc rgb "#" . alpha . "ff7f0e"
set linetype  3 lc rgb "#" . alpha . "2ca02c"; set linetype  4 lc rgb "#" . alpha . "d62728"
set linetype  5 lc rgb "#" . alpha . "9467bd"; set linetype  6 lc rgb "#" . alpha . "8c564b"
set linetype  7 lc rgb "#" . alpha . "e377c2"; set linetype  8 lc rgb "#" . alpha . "7f7f7f"
set linetype  9 lc rgb "#" . alpha . "bcbd22"; set linetype 10 lc rgb "#" . alpha . "17becf"

# set the domain and range
if (exists("xlabel")) set xlabel xlabel
if (exists("ylabel")) set ylabel ylabel
if (exists("xmin")) set xrange [xmin:]
if (exists("xmax")) set xrange [:xmax]
if (exists("ymin")) set yrange [ymin:]
if (exists("ymax")) set yrange [:ymax]

# plot the files
do for [i=2:columns*frames+1:columns] {
    plot for [j=1:ARGC] for [k=0:columns-1] ARGV[j] u 1:(word(scale, k%words(scale)+1)*column(i+k)+word(shift, k%words(shift)+1)) w li s u not, for [i=2:*] background u 1:i w li s u not
}
