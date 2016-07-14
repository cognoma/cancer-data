# Create or overwrite the cognoma-cancer-data conda environment
conda env create --quiet --force --file environment.yml

# Activate the conda environment
source activate cognoma-cancer-data

# Execute notebooks in order
find -maxdepth 1 -iname "*.ipynb" | \
  sort | \
  xargs \
  jupyter nbconvert --execute --inplace --ExecutePreprocessor.timeout=-1
