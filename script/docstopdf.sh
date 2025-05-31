#!/bin/bash

BIB_COMPILER="biber"; GLOSSARY_COMPILER="makeglossaries"; TEX_COMPILER="lualatex"

PAGES=(
    "mathematicalbackground"
    "hilbertspaces"
    "electronicstructuremethods"
    "hartreefock"
    "mollerplessetperturbationtheory"
    "configurationinteraction"
    "coupledcluster"
    "timeevolution"
    "realtimepropagation"
    "imaginarytimepropagation"
    "mixedquantumclassical"
    "trajectorysurfacehopping"
    "fewestswitches"
    "landauzener"
    "mappingapproach"
    "mathematicalmethods"
    "matrixexponential"
    "generalizedeigenvalueproblem"
)

CHAPTERS=(
    "Mathematical Background"
    "Mathematical Methods"
    "Electronic Structure Methods"
    "Time Evolution in Quantum Mechanics"
    "Mixed Quantum-Classical Dynamics"
)

APPENDIX="Mathematical Methods"

ACRONYMS=(
    "Restricted Hartree--Fock/RHF"
    "post-Hartree--Fock/post-HF"
    "Hartree--Fock/HF"
    "Møller--Plesset Perturbation Theory of 2nd Order/MP2"
    "Møller--Plesset Perturbation Theory of 3rd Order/MP3"
    "Møller--Plesset Perturbation Theory/MPPT"
    "Configuration Interaction Singles, Doubles and Triples/CISDT"
    "Configuration Interaction Singles and Doubles/CISD"
    "Full Configuration Interaction/FCI"
    "Configuration Interaction/CI"
    "Coupled Cluster Singles and Doubles/CCSD"
    "Coupled Cluster Doubles/CCD"
    "Coupled Cluster/CC"
    "Time-Dependent Schrödinger Equation/TDSE"
    "Direct Inversion in the Iterative Subspace/DIIS"
    "Molecular Spinorbital/MS"
    "Fast Fourier Transform/FFT"
    "Inverse Fourier Transform/IFT"
    "Discrete Fourier Transform/DFT"
    "Fourier Transform/FT"
    "Self-Consistent Field/SCF"
)

# create the document class in the tex file
rm -rf docs/tex && mkdir docs/tex && cat > docs/tex/main.tex << EOL
% This file was transpiled using a script from the online version in the tjira.github.io/acorn repository. Do not edit it directly.

% Compile with the "$TEX_COMPILER main && $BIB_COMPILER main && $GLOSSARY_COMPILER main && $TEX_COMPILER main && $TEX_COMPILER main" command.

\documentclass[headsepline=true,parskip=half,open=any,11pt]{scrbook}\title{Algorithms of Quantum Chemistry}\author{Tomáš \textsc{Jíra}}
EOL

# add the rest of the header
cat >> docs/tex/main.tex << EOL

\begin{filecontents*}{library.bib}
@article{10.1002/wcms.58,
    author = {Dieter Cremer},
    title = {Møller--Plesset perturbation theory: from small molecule methods to methods for thousands of atoms},
    journal = {WIREs Computational Molecular Science},
    volume = {1},
    pages = {509-530},
    year = {2011},
    doi = {10.1002/wcms.58}
}

@article{10.1063/1.460620,
    author = {John F. Stanton and Jürgen Gauss and John D. Watts and Rodney J. Bartlett},
    title = {A direct product decomposition approach for symmetry exploitation in many-body methods. I. Energy calculations},
    journal = {The Journal of Chemical Physics},
    volume = {94},
    pages = {4334-4345},
    year = {1991},
    doi = {10.1063/1.460620}
}

@article{10.1021/acs.jpclett.9b01641,
    author = {Luke W. Bertels and Joonho Lee and Martin Head-Gordon},
    title = {Third-Order Møller–Plesset Perturbation Theory Made Useful? Choice of Orbitals and Scaling Greatly Improves Accuracy for Thermochemistry, Kinetics, and Intermolecular Interactions},
    journal = {The Journal of Physical Chemistry Letters},
    volume = {10},
    pages = {4170-4176},
    year = {2019},
    doi = {10.1021/acs.jpclett.9b01641}
}

@book{10.1016/S0065-3276!08!60019-2,
    author = {Peter M.W. Gill},
    title = {Molecular integrals Over Gaussian Basis Functions},
    publisher = {Academic Press},
    year = {1994},
    doi = {10.1016/S0065-3276(08)60019-2}
}

@book{1014569052,
    author = {Attila Szabo and Neil S. Ostlund},
    title = {Modern quantum chemistry: introduction to advanced electronic structure theory},
    publisher = {Courier Corporation},
    year = {1996},
    url = {https://www.worldcat.org/title/1014569052}
}

@article{10.1016/0010-4655!73!90016-7,
    author = {Josef Paldus and Heymans H.C. Wong},
    title = {Computer generation of Feynman diagrams for perturbation theory I. General algorithm},
    journal = {Computer Physics Communications},
    volume = {6},
    pages = {1-7},
    year = {1973},
    doi = {10.1016/0010-4655(73)90016-7}
}

@article{10.1016/0010-4655!73!90017-9,
    author = {Heymans H.C. Wong and Josef Paldus},
    title = {Computer generation of Feynman diagrams for perturbation theory II. Program description},
    journal = {Computer Physics Communications},
    volume = {6},
    pages = {9-16},
    year = {1973},
    doi = {10.1016/0010-4655(73)90017-9}
}

@article{10.48550/arXiv.1903.11240,
    author = {Benyamin Ghojogh and Fakhri Karray and Mark Crowley},
    title = {Eigenvalue and Generalized Eigenvalue Problems: Tutorial},
    year = {2019},
    doi = {10.48550/arXiv.1903.11240}
}

@book{10.1002/9780470749593.hrs006,
    author = {Yukio Yamaguchi and Henry F. Schaefer},
    title = {Analytic Derivative Methods in Molecular Electronic Structure Theory: A New Dimension to Quantum Chemistry and its Applications to Spectroscopy},
    publisher = {John Wiley \& Sons, Ltd},
    year = {2011},
    doi = {10.1002/9780470749593.hrs006}
}
\end{filecontents*}

\usepackage[colorlinks=true,linkcolor=blue,pdfa]{hyperref} % hyperlink package should be on top

\usepackage{amsmath} % all the math environments and symbols
\usepackage{amssymb} % additional math symbols
\usepackage[toc,page]{appendix} % appendices
\usepackage[backend=biber,style=chem-acs]{biblatex} % bibliography package
\usepackage{bm} % bold math symbols
\usepackage{braket} % braket notation
\usepackage[format=plain,labelfont=bf]{caption} % bold captions
\usepackage[left=1.5cm,top=2cm,right=1.5cm,bottom=2cm]{geometry} % page layout
\usepackage[acronym,nogroupskip,nomain,toc]{glossaries} % glossary package
\usepackage{mathrsfs} % mathscr environment
\usepackage[automark]{scrlayer-scrpage} % page styles
\usepackage{subcaption} % subfigures

% pdf metadata setup
\makeatletter\hypersetup{
    pdfauthor={\@author},
    pdftitle={\@title},
    pdfsubject={Quantum Chemistry},
    pdfkeywords={quantum,chemistry,algorithm}
}\makeatother

\addbibresource{library.bib} % set the library file
\clearpairofpagestyles % clear the default page styles
\ihead{\pagemark} % add page number to the inner header
\makeglossaries % make the glossary
\ohead{\headmark} % add chapter name to the outer header
EOL

# add the acronym definitions
echo "" >> docs/tex/main.tex && for ACRONYM in "${ACRONYMS[@]}"; do
    IFS='/'; AS=($ACRONYM); unset IFS; echo "\\newacronym{${AS[1],,}}{${AS[1]}}{${AS[0]}}" >> docs/tex/main.tex
done && echo "" >> docs/tex/main.tex

# begin the document
cat >> docs/tex/main.tex << EOL
\begin{document}

\makeatletter\begin{titlepage}
    \center
    \textsc{\LARGE University of Chemistry and Technology, Prague}\\\\[1.5cm]
    \rule{\linewidth}{0.5mm}\\\\[0.4cm]
    {\huge\bfseries\@title}\\\\[0.4cm]
    \rule{\linewidth}{0.5mm}\\\\[1.5cm]
    {\large\textit{Author}}\\\\\@author
    \vfill\vfill\vfill
    {\large\today}
    \vfill
\end{titlepage}\makeatother

\tableofcontents
EOL

# loop over all pages
for PAGE in ${PAGES[@]}; do

    # replace kramdown math wrappers and remove markdown links, comments and hints
    sed 's/\$\$/\$/g ; s/^\$$//g ; s/\[.*\](.*)//g ; s/<!--\|-->//g ; /^[{>]/d' "docs/pages/$PAGE.md" > temp.md

    # convert MD to LaTeX
    pandoc --from=markdown-auto_identifiers --mathjax --to latex --wrap=none --output temp.tex temp.md

    # insert new line after the page and remove temporary
    echo "" >> docs/tex/main.tex && cat temp.tex >> docs/tex/main.tex && rm temp.md temp.tex
done

# replace inline math and bold varepsilon
sed -i 's/\\(\|\\)/$/g ; s/\\mathbf{\\varepsilon}/\\bm{\\varepsilon}/g' docs/tex/main.tex

# set the chapters
for CHAPTER in "${CHAPTERS[@]}"; do
    sed -i "s/\\\\section{\\\\texorpdfstring{$CHAPTER/\\\\chapter{\\\\texorpdfstring{$CHAPTER/g" docs/tex/main.tex
done

# loop over all acronyms
for ACRONYM in "${ACRONYMS[@]}"; do

    # split the full and short acronyms
    IFS='/'; AS=($ACRONYM); unset IFS;

    # replace the full acronyms with the short acronyms and save the result to a temporary file
    awk -v full="${AS[0]}" -v short="${AS[1],,}" 'BEGIN {code = 0} {
        if (/\\section/ || /\\subsection/ || /\\subsubsection/ || /acronym/) {print; next}
        if (/^"""/) {code = !code}; if (!code) {gsub(full, "\\acrshort{"short"}")}; print;
    }' docs/tex/main.tex > temp.tex

    # replace the first occurrence of the short acronyms with the full acronyms and save the result to the main file
    sed -i '0,/\\acrshort{'"${AS[1],,}"'}/s//\\acrfull{'"${AS[1],,}"'}/' temp.tex && mv temp.tex docs/tex/main.tex
done

# set the start of the appendix
sed -i "s/\\\\chapter{\\\\texorpdfstring{Mathematical Methods/\\n\\\\begin{appendices}\\\\chapter{\\\\texorpdfstring{Mathematical Methods/" docs/tex/main.tex

# end the document
cat >> docs/tex/main.tex << EOL
\end{appendices}

\printglossary[type=\acronymtype]

\printbibliography

\end{document}
EOL

# first pass of the compiler
cd docs/tex && $TEX_COMPILER --interaction=batchmode --no-shell-escape main | grep -v "entering extended mode" && cd ../..

# generate the glossary and library
cd docs/tex && $GLOSSARY_COMPILER -q main && $BIB_COMPILER --quiet main && cd ../..

# second and third pass of the compiler
cd docs/tex && $TEX_COMPILER --interaction=batchmode --no-shell-escape main | grep -v "entering extended mode" && cd ../..
cd docs/tex && $TEX_COMPILER --interaction=batchmode --no-shell-escape main | grep -v "entering extended mode" && cd ../..
