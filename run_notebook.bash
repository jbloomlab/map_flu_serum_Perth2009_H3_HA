#!/bin/bash

RESULTSDIR="results"
NB="analysis_notebook.ipynb"

mkdir -p $RESULTSDIR

# run notebook
jupyter nbconvert \
    --to notebook \
    --execute \
    --inplace \
    --ExecutePreprocessor.timeout=-1 \
    $NB

# make Markdown rendering of notebook
jupyter nbconvert \
    --output-dir $RESULTSDIR \
    --to markdown \
    $NB
