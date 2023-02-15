# Rhea Pipeline Shiny App
The Rhea pipeline is a powerful tool for the downstream analysis of 16S rRNA gene amplicon profiles. To make it easier for users to use the Rhea pipeline, we have developed an R-Shiny app that provides a user-friendly interface. This app allows users to upload their OTU tables, metadata, and phylogenetic tree files, select the desired analysis options, and download the results in standard output format.

## Requirements
In order to run the scripts, it is important to have the R language and environment installed first (available here). We strongly recommend using the R-studio to simplify implementation and enhance productivity (available here). Required packages would be automatically installed the first time the scripts are run so internet connection at least for the first time is expected. The scripts are platform independent and should work in all systems that support R.

## Installation
To use the Rhea Pipeline Shiny App, you can follow these steps:

1. Clone or download the repository.
2. Open RStudio and navigate to the downloaded repository.
3. Open the app.R file.
4. Run the app by clicking on the "Run App" button in the top-right corner of the code editor.
5. The necessary packages will be automatically installed and the app will be launched in your browser.

## Usage
To use the Rhea Pipeline Shiny App, you can follow these steps:

1. Upload your OTU table, metadata, and phylogenetic tree files when asked.

2. Run the analysis per tab page.

3. Wait for the analysis to finish.

4. Download the results in standard output format. (A dowload tab is available in the top-right corner of the app. This lets you download the results in a zip file.)

## Features
The Rhea Pipeline Shiny App provides the following features:

- A user-friendly interface for the Rhea pipeline.

- Upload and download data in various formats.

- Various analysis options, including normalization steps, alpha- and beta-diversity analysis, taxonomic composition, statistical comparisons, and calculation of correlations.

- Interactive visualizations of the results.

- Automatic installation of necessary packages.

## Contributions
Contributions to the Rhea Pipeline Shiny App are welcome. If you find a bug, please open an issue, and if you have any suggestions or improvements, please feel free to create a pull request.

## Citation
Rhea was developped by Ilias Lagkouvardos. You can find the paper describing the Rhea pipeline here:
Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. PeerJ 5:e2836 https://doi.org/10.7717/peerj.2836
The github repository for the original Rhea pipeline is available here:
https://github.com/Lagkouvardos/Rhea

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.