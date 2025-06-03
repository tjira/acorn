#!/bin/bash

BIB_COMPILER="biber"; GLOSSARY_COMPILER="makeglossaries"; TEX_COMPILER="lualatex"

PAGES=(
    "mathematical_background"
    "hilbert_space"
    "electronic_structure_methods"
    "hartree_fock"
    "moller_plesset_perturbation_theory"
    "configuration_interaction"
    "coupled_cluster"
    "time_evolution_in_quantum_mechanics"
    "real_time_propagation"
    "imaginary_time_propagation"
    "mixed_quantum_classical_dynamics"
    "electronic_amplitude_propagation"
    "time_derivative_coupling"
    "fewest_switches"
    "landau_zener"
    "mapping_approach"
    "mathematical_methods"
    "matrix_exponential"
    "generalized_eigenvalue_problem"
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
    "Time-Independent Schrödinger Equation/TISE"
    "Direct Inversion in the Iterative Subspace/DIIS"
    "Molecular Spinorbital/MS"
    "Fast Fourier Transform/FFT"
    "Inverse Fourier Transform/IFT"
    "Discrete Fourier Transform/DFT"
    "Fourier Transform/FT"
    "Self-Consistent Field/SCF"
    "Time Derivative Coupling/TDC"
    "Nonadiabatic Coupling Vector/NACV"
    "Trajectory Surface Hopping/TSH"
    "Born--Oppenheimer Approximation/BOA"
    "Potential Energy Surface/PES"
    "Hammes-Schiffer Tully/HST"
    "Norm-Preserving Interpolation/NPI"
)

# create the document class and define the header in the tex file
mkdir -p docs/tex && rm -f docs/tex/* && cat > docs/tex/main.tex << EOL
% This file was transpiled using a script from the online version in the tjira.github.io/acorn repository. Do not edit it directly.

% Compile with the "$TEX_COMPILER main && $BIB_COMPILER main && $GLOSSARY_COMPILER main && $TEX_COMPILER main && $TEX_COMPILER main" command.

\DocumentMetadata{
    lang = en,
    pdfversion = 1.7,
    pdfstandard = A-3b,
}

\documentclass[headsepline=true,parskip=half,open=any,11pt]{scrbook}\title{Algorithms of Quantum Chemistry}\author{Tomáš \textsc{Jíra}}

\usepackage[hidelinks]{hyperref}\makeatletter\hypersetup{pdfauthor={\@author},pdftitle={\@title}}\makeatother % hyperlinks and metadata

\usepackage{amsmath} % all the math environments and symbols
\usepackage{amssymb} % additional math symbols
\usepackage[toc,page]{appendix} % appendices
\usepackage[backend=biber,style=chem-acs]{biblatex} % bibliography package
\usepackage{fontspec} % font selection
\usepackage{braket} % braket notation
\usepackage[format=plain,labelfont=bf]{caption} % bold captions
\usepackage[left=1.5cm,top=2cm,right=1.5cm,bottom=2cm]{geometry} % page layout
\usepackage[acronym,nogroupskip,nomain,toc]{glossaries} % glossary package
\usepackage{mathrsfs} % mathscr environment
\usepackage[automark]{scrlayer-scrpage} % page styles
\usepackage{subcaption} % subfigures
\usepackage{unicode-math} % unicode math support

\addbibresource{library.bib} % set the library file
\clearpairofpagestyles % clear the default page styles
\ihead{\pagemark} % add page number to the inner header
\makeglossaries % make the glossary
\ohead{\headmark} % add chapter name to the outer header
\setmainfont{Libertinus Serif} % set the main font
\setmathfont{Libertinus Math} % set the math font
\setsansfont{Libertinus Sans} % set the sans-serif font
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

    # replace kramdown math wrappers, remove escapes in front of | symbols and remove markdown links, comments and hints
    sed 's/\$\$/\$/g ; s/^\$$//g ; s/\\\\|/|/g ; s/\[.*\](.*)//g ; s/<!--\|-->//g ; /^[{>]/d' "docs/pages/$PAGE.md" > temp.md

    # convert MD to LaTeX
    pandoc --from=markdown-auto_identifiers --mathjax --to latex --wrap=none --output temp.tex temp.md

    # insert new line after the page and remove temporary
    echo "" >> docs/tex/main.tex && cat temp.tex >> docs/tex/main.tex && rm temp.md temp.tex
done

# some post-processing
sed -i 's/ \\eqref/~\\eqref/g' docs/tex/main.tex

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
sed -i "s/\\\\chapter{\\\\texorpdfstring{Mathematical Methods/\\\\begin{appendices}\\\\chapter{\\\\texorpdfstring{Mathematical Methods/" docs/tex/main.tex

# additional sections
cat >> docs/tex/main.tex << EOL
\end{appendices}

\printglossary[type=\acronymtype]

\printbibliography
EOL

# save the bibliography and end the document
echo -e "\n\\\begin{filecontents*}{library.bib}" >> docs/tex/main.tex && cat library.bib >> docs/tex/main.tex && echo -e "\\\end{filecontents*}\n\n\\\end{document}" >> docs/tex/main.tex

# first pass of the compiler
cd docs/tex && $TEX_COMPILER main && cd ../..

# generate the glossary and library
cd docs/tex && $GLOSSARY_COMPILER -q main && $BIB_COMPILER --quiet main && cd ../..

# second and third pass of the compiler
cd docs/tex && $TEX_COMPILER main && cd ../..
cd docs/tex && $TEX_COMPILER main && cd ../..
