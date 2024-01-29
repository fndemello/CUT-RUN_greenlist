# CUT-RUN greenlist

This repository houses the supporting data for the CUT&RUN/Cut&Tag Greenlist, by de Mello et al. 2024:

*The CUT&RUN greenlist: genomic regions of consistent noise are effective normalizing factors for quantitative epigenome mapping
Fabio N. de Mello, Ana C. Tahira, Maria Gabriela Berzoti-Coelho, Sergio Verjovski-Almeida
Briefings in Bioinformatics, Volume 25, Issue 2, March 2024, bbad538, doi: [10.1101/2023.10.26.564165](https://doi.org/10.1101/2023.10.26.564165)*

Greenlist bedfiles are available on the main folder, and also grouped together as an Excel file in 
Supplementary Material S1. R scripts used for greenlist construction and validatation are available 
under /original_scripts.

Previous preprint version available as doi:[10.1101/2023.10.26.564165](https://doi.org/10.1101/2023.10.26.564165), bioRxiv 2023.10.26.564165.

## Instructions of usage

The greenlist is devised to be flexible — we focused on reporting that these regions feature
consistent noise profiles, consistent enough to be used as normalizers. However, the exact
approach on how to use them is up to the reader's preference.

For our manuscript, we chose to quantify the read pile-up on these regions using [deeptools multiBamSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html), followed by calculating [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) size factors for each sample. These 
size factors could then be used along the DESeq2 analysis pipeline.

We offer a brief tutorial for this pipeline as a Jupyter notebook (CUTandRUNAnalysis.ipynb), using the
dataset [GSE157095](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157095) by [Singh et al. (2022)](https://doi.org/10.1126/scitranslmed.abq2096) as an example.
