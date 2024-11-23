# GLQ
A Package for Comparative Phylogeography

This page serves as a repository for the GLQ R package and the dataset used in the study:

Dapporto et al. (2024) The genetic legacy of the Quaternary ice ages for West Palearctic butterflies. Science Advances, 10(38), eadm8596.
Read the full article here.
https://www.science.org/doi/full/10.1126/sciadv.adm8596

Data and Workspace for Replication
To enable replication of the analyses presented in the paper, the complete dataset and accompanying files can be downloaded here:
XXXXXXXXXX

A workspace containing all analyses (excluding FEEMS runs) is also included.
Please note the following considerations:

The workspace was re-run and saved in November 2024, after updates to several R packages.
Some analyses involve random selection of specimens and random ordination of datasets.
As a result, minor differences compared to the results published in the original paper may be observed.

To install GLQ you also need to install recluster and iodatabase. Use:
```
remotes::install_github("leondap/recluster")
remotes::install_github("leondap/GLQ")
```


