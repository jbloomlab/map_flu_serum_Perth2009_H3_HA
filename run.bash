#!/bin/bash

set -e

RESULTSDIR="results/notebooks"

mkdir -p $RESULTSDIR

declare -a nbs=(
                "analyze_map.ipynb"
                "analyze_neut.ipynb"
                "analyze_natseqs.ipynb"
                )

for nb in "${nbs[@]}"
do
    echo "Running $nb"

    jupyter nbconvert \
        --to notebook \
        --execute \
        --inplace \
        --ExecutePreprocessor.timeout=-1 \
        $nb

    echo "Converting $nb to Markdown"
    jupyter nbconvert \
        --output-dir $RESULTSDIR \
        --to markdown \
        $nb
done

map_on_struct="parameterize_map_on_struct.py"
echo "Running $map_on_struct to parameterize notebooks to map on structure."
python $map_on_struct
echo "Completed running $map_on_struct"
echo "Now run the parameterized structure-mapping notebooks interactively."
