# CUT-RUN greenlist

This repository houses the supporting data for the CUT&RUN/Cut&Tag Greenlist, by Mello et al.:

*The CUT&RUN Greenlist: genomic regions of consistent noise are effective normalizing factors for quantitative epigenome mapping
Fabio N. de Mello, Ana C. Tahira, Maria Gabriela Berzoti-Coelho, Sergio Verjovski-Almeida
bioRxiv 2023.10.26.564165; doi: [10.1101/2023.10.26.564165](https://doi.org/10.1101/2023.10.26.564165)*

Greenlist bedfiles are available on the main folder (and also as Supplementary Material S1 to 
the bioRxiv preprint). R scripts used for greenlist construction and validatation are available 
under /original_scripts.

## Instructions of usage

The greenlist is devised to be flexible – we focused on reporting that these regions feature
consistent noise profiles, consistent enough to be used as normalizers. However, the exact
approach on how to use them is up to the reader's preference.

For our manuscript, we chose to quantify the read pile-up on these regions using [deeptools multiBamSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html), followed by calculating [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) size factors for each sample. These 
size factors could then be used along the DESeq2 analysis pipeline.
