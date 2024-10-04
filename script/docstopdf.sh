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
\documentclass[headsepline=true,parskip=half,open=any,11pt]{scrbook}

\usepackage{amsmath} % all the math environments and symbols
\usepackage{braket} % braket notation
\usepackage[left=2cm,top=2.5cm,right=2cm,bottom=2.5cm]{geometry} % page layout
\usepackage[colorlinks=true,linkcolor=blue]{hyperref} % hyperlinks
\usepackage{mathrsfs} % mathscr environment
\usepackage{xcolor} % colors

\usepackage[backend=biber,style=chem-acs]{biblatex} % bibliography
\addbibresource{library.bib}

\usepackage[acronym,automake,nogroupskip,toc]{glossaries} % acronyms
\makeglossaries

\usepackage{listings}
\lstdefinestyle{code}{
    backgroundcolor=\color[rgb]{0.95,0.95,0.92},   
    commentstyle=\color[rgb]{0,0.6,0},
    keywordstyle=\color{magenta},
    numberstyle=\tiny\color[rgb]{0.5,0.5,0.5},
    stringstyle=\color[rgb]{0.58,0,0.82},
    basicstyle=\ttfamily\footnotesize,
    breakatwhitespace=false,         
    breaklines=true,                 
    captionpos=b,                    
    keepspaces=true,                 
    numbers=left,                    
    numbersep=5pt,                  
    showspaces=false,                
    showstringspaces=false,
    showtabs=false,                  
    tabsize=2
}
\lstset{style=code}

\usepackage[automark]{scrlayer-scrpage} % page numbering in header
\clearpairofpagestyles
\ohead{\headmark}
\ihead{\pagemark}

\usepackage{tcolorbox} % colored boxes
\newtcolorbox[auto counter,number within=section]{example}[2][]{
	colback=blue!5!white,
	colframe=blue!75!black,
	fonttitle=\bfseries,
	title=Example~\thetcbcounter: #2,#1
}

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
    CITATION_GROUPS=($(grep -o "cite{.*}--" docs/pages/$PAGE.md | sed 's/cite{// ; s/}--//')); CITATIONS=()

    # split the citation groups
    for CITATION_GROUP in ${CITATION_GROUPS[@]}; do
        IFS=','; CITATIONS+=($CITATION_GROUP); unset IFS;
    done

    # remove duplicates from the citations
    CITATIONS=($(echo "${CITATIONS[@]}" | tr " " "\n" | sort -u | tr "\n" " "))

    # find and append the citation
    for CITATION in ${CITATIONS[@]}; do
        grep $(echo "$CITATION" | sed 's/!/(/ ; s/!/)/') docs/tex/library.md | awk '{print "> " $0 "\n>"}' >> docs/pages/$PAGE.md
    done
done

# remove MD temp
rm docs/tex/*.md
