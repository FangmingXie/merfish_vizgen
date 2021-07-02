#!/bin/bash
# use in the conda env routine

notebook=$1 # .ipynb file
texfile=${notebook/.ipynb/.tex}

outdir="../notes"

# convert to latex
jupyter nbconvert --to latex $notebook
# convert to pdf
xelatex -output-directory $outdir $texfile

# clean up
rm $texfile