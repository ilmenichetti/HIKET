# Overleaf compiles from the project root, while the paper lives in manuscript/.
# Add manuscript/ to the search paths so figures/, the table input and library.bib
# resolve from the root. No effect on local builds run inside manuscript/.
$ENV{'TEXINPUTS'} = './manuscript//:' . ($ENV{'TEXINPUTS'} // '');
$ENV{'BIBINPUTS'} = './manuscript//:' . ($ENV{'BIBINPUTS'} // '');
