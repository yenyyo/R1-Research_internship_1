# Research Internship 1

This repository contains the code and documentation created during my first research internship. The project focuses on immunosequencing to identify signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire, based on the work by Emerson et al.

**Note**: Folders preceded by '00 -' are auxiliary and not essential for understanding the main project.

## Table of Contents
1. [Installation](#installation)
2. [Code Overview](#code-overview)
   - [Data Preparation](#data-preparation)
   - [Clustering](#clustering)
   - [Analysis](#analysis)

## Installation

To get started, follow these steps:

1. Clone the repository.
2. Download the dataset used in the Emerson paper:
   
   Emerson, R., DeWitt, W., Vignali, M. et al. Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire. Nat Genet 49, 659–665 (2017). [https://doi.org/10.1038/ng.3822](https://doi.org/10.1038/ng.3822)
   
3. Follow the installation steps inside `00-Docu/01-Setup/`.

## Code Overview

### Data Preparation

Preprocessing techniques applied to the original dataset can be found in `01-Code/data_preparation/`:

1. `cugraph_clustering`: Original notebook.
2. `data_analysis`: Initial analysis of the dataset.
3. `read_emerson_v0` / `aux_funcs_v0`: Generated the dataset used in later steps. TCRs were defined as a CDR3 sequence.
4. `read_emerson_v2` / `aux_funcs_v2`: Generated the dataset used in later steps. TCRs were defined as a V_gene + CDR3 sequence.

### Clustering

1. `Create_clusters`: Creates the clusters of TCRs based on co-occurrence.
2. `Analyse_clusters` / `aux_clustering`: Loads data to compare and filter different clusters. Generates the file `gmm_results.pkl`.

### Analysis

1. `hla_analysis`: Analyzes each of the groups that divide each cluster and extracts HLA differences between both.
