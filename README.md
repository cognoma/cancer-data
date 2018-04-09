# Cancer data acquisition and processing for Project Cognoma

This is a mixed notebook and data repository for retrieving cancer data for [Project Cognoma](https://github.com/cognoma/cognoma).
Currently, all data is from the [TCGA Pan-Cancer collection](https://xenabrowser.net/datapages/?cohort=TCGA%20PanCanAtlas "UCSC Xena Browser cohort: TCGA PanCanAtlas") of the UCSC Xena Browser.

## Workflow

The data acquisition and analysis is executing by running Jupyter notebooks in the following order:

+ [`1.TCGA-download.ipynb`](1.TCGA-download.ipynb) — download and compress TCGA datasets.
+ [`2.TCGA-process.ipynb`](2.TCGA-process.ipynb) — convert downloaded TCGA datasets into sample × gene matrixes.

The [`execute.sh`](execute.sh) script executes the notebooks in order.
After installing and activating the [environment](#environment), run with the command `bash execute.sh` from the repository's root directory.

## Directories

The repository contains the following directories:

+ [`download`](download) — contains files retrieved from an external location whose content is unmodified.
Large downloaded files are tracked using Git LFS.
Associated metadata files are also retained for versioning.
+ [`data`](data) — contains generated datasets.
The complete matrix files are not currently tracked due to file size, but randomly-subsetted versions are available for development use (see [`data/subset`](data/subset)).

## Download

[![DOI: 10.6084/m9.figshare.3487685](https://img.shields.io/badge/DOI-10.6084/m9.figshare.3487685-blue.svg)](https://doi.org/10.6084/m9.figshare.3487685 "Complete datasets on figshare")

The complete datasets created by this repository (`data/expression-matrix.tsv.bz2` and `data/mutation-matrix.tsv.bz2`) are uploaded to [figshare](https://doi.org/10.6084/m9.figshare.3487685).
Since this is a manual process, check the figshare REFERENCES section to see which commit these datasets derive from.
In other words, the latest version on figshare may lag behind this repository.

## Environment

This repository uses [conda](https://conda.io/docs/) to manage its environment, which is named `cognoma-cancer-data`.
The required packages and versions are listed in [`environment.yml`](environment.yml).
If as a developer, you require an additional package, add it to `environment.yml`.

The following commands install and activate the environment:

```sh
# Create or overwrite the cognoma-cancer-data conda environment
conda env create --file=environment.yml

# Activate the conda environment (assumes conda >= 4.4)
conda activate cognoma-cancer-data
```

## License

This repository is dual licensed as [BSD 3-Clause](LICENSE-BSD.md) and [CC0 1.0](LICENSE-CC0.md), meaning any repository content can be used under either license.
This licensing arrangement ensures source code is available under an [OSI-approved License](https://opensource.org/licenses/alphabetical), while non-code content — such as figures, data, and documentation — is maximally reusable under a public domain dedication.
