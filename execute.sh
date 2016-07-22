#!/bin/bash

# Exit on error
set -o errexit

# Create or overwrite the cognoma-cancer-data conda environment
conda env create --quiet --force --file environment.yml

# Activate the conda environment
source activate cognoma-cancer-data

# Execute notebooks in order
jupyter nbconvert --inplace --execute --ExecutePreprocessor.timeout=-1 *.ipynb
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
