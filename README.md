# Nature Immunology Review

This repository contains code for Figure 2 in the manuscript "Antigen-presenting cells: arbiters of mucosal tolerance and immunity."

## Overview

The repo contains code for analyzing and integrating 9 publicly available scRNA-seq datasets of ROR𝛾t+ APCs.

## Setup

Code in this repo requires R 4.2.3 and Python 3.8 as well as the following dependencies:
- Seurat v4.4.0
- ArchR
- dplyr
- ggplot2
- Scanpy v1.9.3
- Palantir v1.2

## Data availability

The datasets used in this repo are available at the following locations:

| Dataset Name             | Data accession                        | Publication                                                           |
|------------------------|------------------------|------------------------|
| MLN_RORgt_MHCII_multiome | GSE174405                             | <https://www.nature.com/articles/s41586-022-05309-5>                  |
| TC_all_LN                | GSE294005                             | <https://www.science.org/doi/10.1126/science.adp0535>                 |
| Colonna                  | GSE289268                             | <https://www.sciencedirect.com/science/article/pii/S0092867425002934> |
| Kedmi                    | GSE200148                             | <https://www.nature.com/articles/s41586-022-05089-y>                  |
| Wang                     | GSE176282                             | <https://www.science.org/doi/10.1126/sciimmunol.abl5053>              |
| Littman                  | <https://zenodo.org/records/15212000> | <https://www.nature.com/articles/s41586-025-08982-4>                  |
| Lyu                      | GSE184175                             | <https://www.nature.com/articles/s41586-022-05141-x>                  |
| Gardner.A                | GSE273746                             | <https://doi.org/10.1084/jem.20250573>                                |
| Gardner.E                | GSE285182                             | <https://doi.org/10.1084/jem.20250573>                                |

The results of integration including Seurat R object are available at [Box](https://mskcc.box.com/s/0z5mdy8e1zkmo4m5prj56hjcgse9l0ya).
