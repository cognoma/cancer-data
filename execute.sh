#!/bin/bash

# Exit on error
set -o errexit

# Execute notebooks in order
jupyter nbconvert --inplace --execute \
  --ExecutePreprocessor.timeout=-1 \
  --ExecutePreprocessor.kernel_name=python \
  *.ipynb

# Export to .py files
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
