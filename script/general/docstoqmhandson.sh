#!/bin/bash

# compile the pdf
./script/general/docstopdf.sh

# create the folder and copy the library
mkdir -p docs/tex/qmhandson && cp docs/tex/library.bib docs/tex/qmhandson/

# extract the electronic structure
sed -n '/tableofcontents/,/appendices/{/tableofcontents/b;/appendices/b;p}' docs/tex/main.tex > docs/tex/qmhandson/elstructure.tex

# extract the code solutions
sed -n '/\\label{sec:hf_code_solution}/,/appendices/{/appendices/b;p}' docs/tex/main.tex > docs/tex/qmhandson/elstructure_solutions.tex

# remove first ald last line of the files
sed -i '1d;$d' docs/tex/qmhandson/elstructure.tex && sed -i '$d' docs/tex/qmhandson/elstructure_solutions.tex
