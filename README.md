# code-for-masters-thesis
This repository contains the code for my Master’s Thesis and my Master's Thesis paper in the respective uploaded files:
Final_script_for_masters_project.R
and
Master's Thesis.pdf

The initial file that is being loaded in the script is the unzipped GSE217215_201218_ATAC_subsample.rds file containing (meta)data of hESCs that were transduced with a library of 198 transcription factor ORFs and cultured for 7 days before single-cell ATAC- and RNA-seq joint profiling. The other file needed in order to run the script is the GSE217215_201218_ATAC_fragments.tsv.gz compressed file that contains the raw scATAC-seq fragments file.
These files are too large for this git repository, but they are available at the Gene Expression Omnibus database (GEO) with the GEO accession number: GSE217215 and it is from the paper by Joung J, Ma S, Tay T, Geiger-Schuller KR et al. A transcription factor atlas of directed differentiation. Cell 2023 Jan 5;186(1):209-229.e26. PMID: 36608654. It is recommended to download the https versions of them.

Make sure that the paths to these files are appropriately changed in the script in order to run it and that the appropriate packages/libraries in the script are loaded.
