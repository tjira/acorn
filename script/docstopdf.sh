#!/bin/bash

PAGES=("hartreefockmethod" "mollerplessetperturbationtheory" "configurationinteraction" "coupledcluster")

ACRONYMS=(
    "Restricted Hartree--Fock/RHF"
    "post-Hartree--Fock/post-HF"
    "Hartree--Fock/HF"
    "Møller--Plesset Perturbation Theory/MPPT"
    "Configuration Interaction Singles, Doubles and Triples/CISDT"
    "Configuration Interaction Singles and Doubles/CISD"
    "Full Configuration Interaction/FCI"
    "Configuration Interaction/CI"
    "Coupled Cluster Singles and Doubles/CCSD"
    "Coupled Cluster Doubles/CCD"
    "Coupled Cluster/CC"
    "Self-Consistent Field/SCF"
    "Molecular Spinorbital/MS"
)

# create the main LaTeX document
cat > docs/tex/main.tex << EOL
\documentclass[open=any,parskip=half,11pt]{scrbook}

\usepackage{amsmath}
\usepackage{braket}
\usepackage[left=2cm,top=2.5cm,right=2cm,bottom=2.5cm]{geometry}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\usepackage{mathrsfs}

\usepackage[backend=biber,style=chem-acs]{biblatex}
\addbibresource{library.bib}

\usepackage[acronym,automake,nogroupskip,toc]{glossaries}
\makeglossaries

\title{Algorithms of Quantum Chemistry}
\author{Tom\'a\v s J\'ira}
EOL

# add the acronym definitions
echo "" >> docs/tex/main.tex && for ACRONYM in "${ACRONYMS[@]}"; do
    IFS='/'; AS=($ACRONYM); unset IFS; echo "\\newacronym{${AS[1],,}}{${AS[1]}}{${AS[0]}}" >> docs/tex/main.tex
done && echo "" >> docs/tex/main.tex

# begin the document
cat >> docs/tex/main.tex << EOL
\begin{document}

\maketitle
\tableofcontents
EOL

# loop over all pages
echo "" >> docs/tex/main.tex && for PAGE in ${PAGES[@]}; do

    # replace some MD quirks with LaTeX quirks
    cat docs/pages/$PAGE.md | sed 's/\\\\\\\\\\/\\\\/g ; s/\\_/_/g ; s/\\|/|/g ; s/<!--//g ; s/-->//g ; /mathjax/d ; /{:.cite}/d ; /{:.note}/d ; /> /d' > temp.md

    # convert MD to LaTeX
    pandoc temp.md --mathjax --wrap=none -o temp.tex && cat temp.tex >> docs/tex/main.tex && rm temp.md temp.tex
done && echo "" >> docs/tex/main.tex

# end the document
cat >> docs/tex/main.tex << EOL
\printglossary[type=\acronymtype,title=List of Acronyms,toctitle=List of Acronyms]

\printbibliography

\end{document}
EOL

# replace section references
sed -i 's/\\href{hartreefockmethod.html\\#the-integral-transforms}{here}/in section \\ref{integral-transforms-to-the-basis-of-molecular-spinorbitals}/g' docs/tex/main.tex

# loop over all acronyms
for ACRONYM in "${ACRONYMS[@]}"; do

    # split the full and short acronyms
    IFS='/'; AS=($ACRONYM); unset IFS;

    # replace the full acronyms with the short acronyms and save the result to a temporary file
    awk -v full="${AS[0]}" -v short="${AS[1],,}" '(/\\section/ || /\\subsection/ || /\\subsubsection/ || /acronym/) {print; next} {gsub(full, "\\acrshort{"short"}"); print}' docs/tex/main.tex > temp.tex

    # replace the first occurrence of the short acronyms with the full acronyms and save the result to the main file
    sed -i '0,/\\acrshort{'"${AS[1],,}"'}/s//\\acrfull{'"${AS[1],,}"'}/' temp.tex && mv temp.tex docs/tex/main.tex
done

# replace section with chapters and subsections with sections
sed -i 's/\\section/\\chapter/g ; s/\\subsection/\\section/g ; s/\\subsubsection/\\subsection/g' docs/tex/main.tex

# compile the main LaTeX file to PDF
cd docs/tex && pdflatex main && biber main && makeglossaries main && pdflatex main && pdflatex main && cd ../..

# create the library LaTeX document
cat > docs/tex/library.tex << EOL
\documentclass{scrbook}

\usepackage[left=2cm,top=2.5cm,right=2cm,bottom=2.5cm]{geometry}

\usepackage[backend=biber,doi=true,style=chem-acs]{biblatex}
\addbibresource{library.bib}

\begin{document}

\nocite{*}

\printbibliography[heading=none]

\end{document}
EOL

# compile the library LaTeX file to PDF
cd docs/tex && pdflatex library && biber library && pdflatex library && pdflatex library && cd ../..

# create the MD file from the library LaTeX file
cd docs/tex && pandoc -t markdown_strict --wrap=none --citeproc library.tex -o library.md --bibliography library.bib && cd ../..

# loop over all pages
for PAGE in ${PAGES[@]}; do

    # remove old citations
    awk '/{:.cite}/ {exit} {print}' docs/pages/$PAGE.md > temp.md && mv temp.md docs/pages/$PAGE.md && echo "{:.cite}" >> docs/pages/$PAGE.md

    # get the citations
    CITATIONS=($(grep -o "cite{.*}" docs/pages/$PAGE.md | sed 's/cite{// ; s/}//'))

    # remove duplicates from the citations
    CITATIONS=($(echo "${CITATIONS[@]}" | tr " " "\n" | sort -u | tr "\n" " "))

    # find and append the citation
    for CITATION in ${CITATIONS[@]}; do
        grep $CITATION docs/tex/library.md | awk '{print "> " $0 "\n>"}' >> docs/pages/$PAGE.md
    done
done

# remove MD library
rm docs/tex/library.md
